use std::path::PathBuf;

use anyhow::{anyhow, bail};
use molar::prelude::*;
use numpy::{
    nalgebra::{self, Const, Dyn, VectorView},
    PyArrayLike1, PyArrayMethods, PyReadonlyArray2, PyUntypedArrayMethods, ToPyArray,
};
use pyo3::{prelude::*, types::PyTuple, IntoPyObjectExt};

mod utils;
use triomphe::Arc;
use utils::*;

mod atom;
use atom::Atom;

mod particle;
use particle::Particle;

mod periodic_box;
use periodic_box::PeriodicBox;

mod membrane;
use membrane::*;

//-------------------------------------------

#[pyclass(unsendable)]
struct Topology(triomphe::Arc<molar::core::Topology>);

#[pyclass(unsendable)]
struct State(triomphe::Arc<molar::core::State>);

#[pymethods]
impl State {
    fn __len__(&self) -> usize {
        self.0.len()
    }

    #[getter]
    fn get_time(&self) -> f32 {
        self.0.get_time()
    }

    #[setter]
    fn set_time(&self, t: f32) {
        self.0.set_time(t);
    }

    #[getter]
    fn get_box(&self) -> anyhow::Result<PeriodicBox> {
        Ok(PeriodicBox(
            self.0
                .get_box()
                .ok_or_else(|| anyhow!("No periodic box"))?
                .clone(),
        ))
    }

    #[setter]
    fn set_box(&mut self, val: Bound<'_, PeriodicBox>) -> anyhow::Result<()> {
        let b = self
            .0
            .get_box_mut()
            .ok_or_else(|| anyhow!("No periodic box"))?;
        *b = val.borrow().0.clone();
        Ok(())
    }
}

#[pyclass(unsendable)]
struct FileHandler(
    Option<molar::io::FileHandler>,
    Option<molar::io::IoStateIterator>,
);

const ALREADY_TRANDFORMED: &str = "file handler is already transformed to state iterator";

#[pymethods]
impl FileHandler {
    #[new]
    fn new(fname: &str, mode: &str) -> anyhow::Result<Self> {
        match mode {
            "r" => Ok(FileHandler(
                Some(molar::io::FileHandler::open(fname)?),
                None,
            )),
            "w" => Ok(FileHandler(
                Some(molar::io::FileHandler::create(fname)?),
                None,
            )),
            _ => Err(anyhow!("Wrong file open mode")),
        }
    }

    fn read(&mut self) -> anyhow::Result<(Topology, State)> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        let (top, st) = h.read()?;
        Ok((Topology(top.into()), State(st.into())))
    }

    fn read_topology(&mut self) -> anyhow::Result<Topology> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        let top = h.read_topology()?;
        Ok(Topology(top.into()))
    }

    fn read_state(&mut self) -> anyhow::Result<State> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        let st = h.read_state()?;
        Ok(State(st.into()))
    }

    fn write(&mut self, data: Bound<'_, PyAny>) -> anyhow::Result<()> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        if let Ok(s) = data.extract::<PyRef<'_, System>>() {
            h.write(&s.0)?;
        } else if let Ok(s) = data.extract::<PyRef<'_, Sel>>() {
            h.write(&s.0)?;
        } else if let Ok(s) = data.downcast::<PyTuple>() {
            if s.len() != 2 {
                return Err(anyhow!("Tuple must have two elements"));
            }
            let top = s
                .iter()
                .next()
                .unwrap()
                .extract::<PyRefMut<'_, Topology>>()?;
            let st = s.iter().next().unwrap().extract::<PyRefMut<'_, State>>()?;
            h.write(&(Arc::clone(&top.0), Arc::clone(&st.0)))?;
        } else {
            return Err(anyhow!(
                "Invalid data type {} when writing to file",
                data.get_type()
            ));
        }
        Ok(())
    }

    fn write_topology(&mut self, data: Bound<'_, PyAny>) -> anyhow::Result<()> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        if let Ok(s) = data.extract::<PyRef<'_, System>>() {
            h.write_topology(&s.0)?;
        } else if let Ok(s) = data.extract::<PyRef<'_, Sel>>() {
            h.write_topology(&s.0)?;
        } else if let Ok(s) = data.extract::<PyRefMut<'_, Topology>>() {
            h.write_topology(&s.0)?;
        } else {
            return Err(anyhow!(
                "Invalid data type {} when writing to file",
                data.get_type()
            ));
        }
        Ok(())
    }

    fn write_state(&mut self, data: Bound<'_, PyAny>) -> anyhow::Result<()> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        if let Ok(s) = data.extract::<PyRef<'_, System>>() {
            h.write_state(&s.0)?;
        } else if let Ok(s) = data.extract::<PyRef<'_, Sel>>() {
            h.write_state(&s.0)?;
        } else if let Ok(s) = data.extract::<PyRefMut<'_, State>>() {
            h.write_state(&s.0)?;
        } else {
            return Err(anyhow!(
                "Invalid data type {} when writing to file",
                data.get_type()
            ));
        }
        Ok(())
    }

    fn __iter__(mut slf: PyRefMut<'_, Self>) -> PyRefMut<'_, Self> {
        if slf.1.is_none() {
            let h = slf.0.take().unwrap();
            slf.1 = Some(h.into_iter());
        }
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PyObject> {
        let st = slf.1.as_mut().unwrap().next().map(|st| State(st.into()));
        if st.is_some() {
            Python::with_gil(|py| Some(st.unwrap().into_py_any(py)))
                .unwrap()
                .ok()
        } else {
            None
        }
    }

    fn skip_to_frame(&mut self, fr: usize) -> PyResult<()> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        h.skip_to_frame(fr).map_err(|e| anyhow!(e))?;
        Ok(())
    }

    fn skip_to_time(&mut self, t: f32) -> anyhow::Result<()> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        h.skip_to_time(t)?;
        Ok(())
    }

    fn tell_first(&self) -> anyhow::Result<(usize, f32)> {
        let h = self
            .0
            .as_ref()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        Ok(h.tell_first()?)
    }

    fn tell_current(&self) -> anyhow::Result<(usize, f32)> {
        let h = self
            .0
            .as_ref()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        Ok(h.tell_current()?)
    }

    fn tell_last(&self) -> anyhow::Result<(usize, f32)> {
        let h = self
            .0
            .as_ref()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        Ok(h.tell_last()?)
    }

    #[getter]
    fn stats(&self) -> anyhow::Result<FileStats> {
        let h = self
            .0
            .as_ref()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        Ok(FileStats(h.stats.clone()))
    }

    #[getter]
    fn file_name(&self) -> anyhow::Result<PathBuf> {
        let h = self
            .0
            .as_ref()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        Ok(h.file_path.clone())
    }
}

#[pyclass]
struct FileStats(molar::io::FileStats);

#[pymethods]
impl FileStats {
    #[getter]
    fn elapsed_time(&self) -> std::time::Duration {
        self.0.elapsed_time
    }

    #[getter]
    fn frames_processed(&self) -> usize {
        self.0.frames_processed
    }

    #[getter]
    fn cur_t(&self) -> f32 {
        self.0.cur_t
    }

    fn __repr__(&self) -> String {
        format!("{}", self.0)
    }

    fn __str__(&self) -> String {
        format!("{}", self.0)
    }
}

#[pyclass(unsendable, sequence)]
struct System(molar::core::System);

#[pymethods]
impl System {
    #[new]
    #[pyo3(signature = (*py_args))]
    fn new<'py>(py_args: &Bound<'py, PyTuple>) -> PyResult<Self> {
        if py_args.len() == 1 {
            // From file
            Ok(System(
                molar::core::System::from_file(&py_args.get_item(0)?.extract::<String>()?)
                    .map_err(|e| anyhow!(e))?,
            ))
        } else if py_args.len() == 2 {
            let top = py_args
                .get_item(0)?
                .downcast::<Topology>()?
                .try_borrow_mut()?;
            let st = py_args.get_item(1)?.downcast::<State>()?.try_borrow_mut()?;
            Ok(System(
                molar::core::System::new(Arc::clone(&top.0), Arc::clone(&st.0))
                    .map_err(|e| anyhow!(e))?,
            ))
        } else {
            // Empty builder
            Ok(System(molar::core::System::new_empty()))
        }
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn select_all(&mut self) -> anyhow::Result<Sel> {
        Ok(Sel::new_owned(self.0.select_all()?))
    }

    fn select(&mut self, sel_str: &str) -> anyhow::Result<Sel> {
        Ok(Sel::new_owned(self.0.select(sel_str)?))
    }

    #[pyo3(signature = (arg=None))]
    fn __call__(&self, arg: Option<&Bound<'_, PyAny>>) -> anyhow::Result<Sel> {
        if let Some(arg) = arg {
            if let Ok(val) = arg.extract::<String>() {
                if val.is_empty() {
                    Ok(Sel::new_owned(self.0.select_all()?))
                } else {
                    Ok(Sel::new_owned(self.0.select(val)?))
                }
            } else if let Ok(val) = arg.extract::<(usize, usize)>() {
                Ok(Sel::new_owned(self.0.select(val.0..val.1)?))
            } else if let Ok(val) = arg.extract::<Vec<usize>>() {
                Ok(Sel::new_owned(self.0.select(val)?))
            } else {
                Err(anyhow!(
                    "Invalid argument type {} when creating selection",
                    arg.get_type()
                )
                .into())
            }
        } else {
            Ok(Sel::new_owned(self.0.select_all()?))
        }
    }

    fn set_state(&mut self, st: &State) -> anyhow::Result<State> {
        let old_state = self.0.set_state(Arc::clone(&st.0))?;
        Ok(State(old_state))
    }

    fn set_topology(&mut self, top: &Topology) -> anyhow::Result<Topology> {
        let old_top = self.0.set_topology(Arc::clone(&top.0))?;
        Ok(Topology(old_top))
    }

    // fn get_state(&self) -> State {
    //     State(self.0.get_state())
    // }

    // fn get_topology(&self) -> Topology {
    //     Topology(self.0.get_topology())
    // }

    fn save(&self, fname: &str) -> anyhow::Result<()> {
        Ok(self.0.save(fname)?)
    }

    fn remove(&self, arg: &Bound<'_, PyAny>) -> anyhow::Result<()> {
        // In the future other types can be used as well
        if let Ok(sel) = arg.downcast::<Sel>() {
            Ok(self.0.remove(&sel.borrow().0)?)
        } else if let Ok(sel_str) = arg.extract::<String>() {
            let sel = self.0.select(sel_str)?;
            Ok(self.0.remove(&sel)?)
        } else if let Ok(list) = arg.extract::<Vec<usize>>() {
            Ok(self.0.remove(&list)?)
        } else {
            unreachable!()
        }
    }

    fn append(&self, arg: &Bound<'_, PyAny>) -> anyhow::Result<()> {
        // In the future other types can be used as well
        if let Ok(sel) = arg.downcast::<Sel>() {
            self.0.append(&sel.borrow().0);
        } else if let Ok(sel) = arg.downcast::<System>() {
            self.0.append(&sel.borrow().0);
        } else {
            anyhow::bail!("Unsupported type to append a Source")
        }
        Ok(())
    }

    #[getter]
    fn get_time(&self) -> f32 {
        self.0.get_time()
    }

    // #[setter]
    // fn set_time(&self, t: f32) {
    //     self.0.set_time(t);
    // }

    fn set_box_from(&self, sys: &System) {
        self.0.set_box_from(&sys.0);
    }
}

//====================================

#[pyclass(sequence, unsendable)]
struct Sel(molar::core::Sel);

impl Sel {
    fn new_owned(sel: molar::core::Sel) -> Self {
        Self(sel)
    }

    fn new_ref(sel: &molar::core::Sel) -> Self {
        Self(sel.new_view())
    }
}

#[pymethods]
impl Sel {
    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn __call__(&self, arg: &Bound<'_, PyAny>) -> PyResult<Sel> {
        if let Ok(val) = arg.extract::<String>() {
            Ok(Sel::new_owned(self.0.select(val).map_err(|e| anyhow!(e))?))
        } else if let Ok(val) = arg.extract::<(usize, usize)>() {
            Ok(Sel::new_owned(
                self.0.select(val.0..=val.1).map_err(|e| anyhow!(e))?,
            ))
        } else if let Ok(val) = arg.extract::<Vec<usize>>() {
            Ok(Sel::new_owned(self.0.select(val).map_err(|e| anyhow!(e))?))
        } else {
            Err(anyhow!(
                "Invalid argument type {} when creating selection",
                arg.get_type()
            )
            .into())
        }
    }

    // Indexing
    fn __getitem__(slf: Bound<Self>, i: isize) -> PyResult<Py<PyAny>> {
        let s = slf.borrow();
        let ind = if i < 0 {
            if i.abs() > s.__len__() as isize {
                return Err(anyhow!(
                    "Negative index {i} is out of bounds {}:-1",
                    -(s.__len__() as isize)
                )
                .into());
            }
            s.__len__() - i.unsigned_abs()
        } else if i >= s.__len__() as isize {
            return Err(anyhow!("Index {} is out of bounds 0:{}", i, s.__len__()).into());
        } else {
            i as usize
        };

        // Call Rust function
        let p = s.0.get_particle_mut(ind).unwrap();
        Ok(Particle {
            atom: unsafe { &mut *(p.atom as *mut molar::core::Atom) },
            pos: map_pyarray_to_pos(slf.py(), p.pos, &slf),
            id: p.id,
        }
        .into_py_any(slf.py())?)
    }

    // Iteration protocol
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, ParticleIterator> {
        Bound::new(
            slf.py(),
            ParticleIterator {
                sel: slf.into(),
                cur: 0,
            },
        )
        .unwrap()
        .borrow()
    }

    fn get_index<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray1<usize>> {
        numpy::PyArray1::from_iter(py, self.0.iter_index())
    }

    fn get_coord<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f32>> {
        // We allocate an uninitialized PyArray manually and fill it with data.
        // By doing this we save on unnecessary initiallization and extra allocation
        unsafe {
            let arr = numpy::PyArray2::<f32>::new(py, [3, self.0.len()], true);
            let arr_ptr = arr.data();
            for i in 0..self.0.len() {
                let pos_ptr = self.0.get_pos_unchecked(i).coords.as_ptr();
                // This is faster than copying by element with uget_raw()
                std::ptr::copy_nonoverlapping(pos_ptr, arr_ptr.offset(i as isize * 3), 3);
            }
            arr
        }
    }

    fn set_coord(&mut self, arr: PyReadonlyArray2<f32>) -> PyResult<()> {
        // Check if the shape is correct
        if arr.shape() != [3, self.__len__()] {
            return Err(anyhow!(
                "Array shape must be [3, {}], not {:?}",
                self.__len__(),
                arr.shape()
            ))?;
        }
        let ptr = arr.data();

        unsafe {
            for i in 0..self.__len__() {
                let pos_ptr = self.0.get_pos_mut_unchecked(i).coords.as_mut_ptr();
                std::ptr::copy_nonoverlapping(ptr.offset(i as isize * 3), pos_ptr, 3);
            }
        }

        Ok(())
    }

    fn set_state(&mut self, st: &State) -> anyhow::Result<State> {
        let old_state = self.0.set_state(Arc::clone(&st.0))?;
        Ok(State(old_state))
    }

    fn set_state_from(&mut self, arg: &Bound<'_, PyAny>) -> anyhow::Result<State> {
        if let Ok(val) = arg.downcast::<System>() {
            Ok(State(self.0.set_state_from(&val.borrow().0)?))
        } else if let Ok(val) = arg.downcast::<Sel>() {
            Ok(State(self.0.set_state_from(&val.borrow().0)?))
        } else {
            Err(anyhow!(
                "Invalid argument type {} in set_state_from()",
                arg.get_type()
            )
            .into())
        }
    }

    fn set_topology(&mut self, top: &Topology) -> anyhow::Result<Topology> {
        let old_top = self.0.set_topology(Arc::clone(&top.0))?;
        Ok(Topology(old_top))
    }

    pub fn set_same_chain(&self, val: char) {
        self.0.set_same_chain(val)
    }

    pub fn set_same_resname(&mut self, val: &str) {
        self.0.set_same_resname(val)
    }

    pub fn set_same_resid(&mut self, val: i32) {
        self.0.set_same_resid(val)
    }

    pub fn set_same_name(&mut self, val: &str) {
        self.0.set_same_name(val)
    }

    pub fn set_same_mass(&mut self, val: f32) {
        self.0.set_same_mass(val)
    }

    pub fn set_same_bfactor(&mut self, val: f32) {
        self.0.set_same_bfactor(val)
    }

    #[getter]
    fn get_time(&self) -> f32 {
        self.0.get_time()
    }

    #[pyo3(signature = (dims=[false,false,false]))]
    fn com<'py>(
        &self,
        py: Python<'py>,
        dims: [bool; 3],
    ) -> PyResult<Bound<'py, numpy::PyArray1<f32>>> {
        let pbc_dims = PbcDims::new(dims[0], dims[1], dims[2]);
        Ok(clone_vec_to_pyarray1(
            &self
                .0
                .center_of_mass_pbc_dims(pbc_dims)
                .map_err(|e| anyhow!(e))?
                .coords,
            py,
        ))
    }

    #[pyo3(signature = (dims=[false,false,false]))]
    fn cog<'py>(
        &self,
        py: Python<'py>,
        dims: [bool; 3],
    ) -> PyResult<Bound<'py, numpy::PyArray1<f32>>> {
        let pbc_dims = PbcDims::new(dims[0], dims[1], dims[2]);
        Ok(clone_vec_to_pyarray1(
            &self
                .0
                .center_of_geometry_pbc_dims(pbc_dims)
                .map_err(|e| anyhow!(e))?
                .coords,
            py,
        ))
    }

    fn principal_transform(&self) -> anyhow::Result<IsometryTransform> {
        let tr = self.0.principal_transform()?;
        Ok(IsometryTransform(tr))
    }

    fn principal_transform_pbc(&self) -> anyhow::Result<IsometryTransform> {
        let tr = self.0.principal_transform_pbc()?;
        Ok(IsometryTransform(tr))
    }

    fn apply_transform(&self, tr: &IsometryTransform) {
        self.0.apply_transform(&tr.0);
    }

    fn gyration(&self) -> anyhow::Result<f32> {
        Ok(self.0.gyration()?)
    }

    fn gyration_pbc(&self) -> anyhow::Result<f32> {
        Ok(self.0.gyration_pbc()?)
    }

    fn inertia<'py>(
        &self,
        py: Python<'py>,
    ) -> anyhow::Result<(
        Bound<'py, numpy::PyArray1<f32>>,
        Bound<'py, numpy::PyArray2<f32>>,
    )> {
        let (moments, axes) = self.0.inertia()?;
        let mom = clone_vec_to_pyarray1(&moments, py);
        let ax = axes.to_pyarray(py);
        Ok((mom, ax))
    }

    fn inertia_pbc<'py>(
        &self,
        py: Python<'py>,
    ) -> anyhow::Result<(
        Bound<'py, numpy::PyArray1<f32>>,
        Bound<'py, numpy::PyArray2<f32>>,
    )> {
        let (moments, axes) = self.0.inertia_pbc()?;
        let mom = clone_vec_to_pyarray1(&moments, py);
        let ax = axes.to_pyarray(py);
        Ok((mom, ax))
    }

    fn save(&self, fname: &str) -> anyhow::Result<()> {
        Ok(self.0.save(fname)?)
    }

    fn translate<'py>(&self, arg: PyArrayLike1<'py, f32>) -> anyhow::Result<()> {
        let vec: VectorView<f32, Const<3>, Dyn> = arg
            .try_as_matrix()
            .ok_or_else(|| anyhow!("conversion to Vector3 has failed"))?;
        self.0.translate(&vec);
        Ok(())
    }

    fn split_resindex(&self) -> Vec<Sel> {
        self.0.split_resindex_iter().map(|s| Sel(s)).collect()
    }

    fn split_chain(&self) -> Vec<Sel> {
        self.0
            .split_iter(|p| Some(p.atom.chain))
            .map(|s| Sel(s))
            .collect()
    }

    fn split_molecule(&self) -> Vec<Sel> {
        self.0
            .split_mol_iter()
            .map(|sel| Sel::new_owned(sel))
            .collect()
    }

    fn to_gromacs_ndx(&self, name: &str) -> String {
        self.0.as_gromacs_ndx_str(name)
    }

    /// operator |
    fn __or__(&self, rhs: &Sel) -> Sel {
        Sel::new_owned(&self.0 | &rhs.0)
    }

    /// operator &
    fn __and__(&self, rhs: &Sel) -> Sel {
        Sel::new_owned(&self.0 & &rhs.0)
    }

    /// -= (remove other from self)
    fn __sub__(&self, rhs: &Sel) -> Sel {
        Sel::new_owned(&self.0 - &rhs.0)
    }

    /// ~ operator
    fn __invert__(&self) -> Sel {
        Sel::new_owned(!&self.0)
    }

    fn sasa(&self) -> SasaResults {
        SasaResults(self.0.sasa())
    }
}

#[pyclass(unsendable)]
struct SasaResults(molar::core::SasaResults);

#[pymethods]
impl SasaResults {
    #[getter]
    fn areas(&self) -> &[f32] {
        self.0.areas()
    }

    #[getter]
    fn volumes(&self) -> &[f32] {
        self.0.volumes()
    }

    #[getter]
    fn total_area(&self) -> f32 {
        self.0.total_area()
    }

    #[getter]
    fn total_volume(&self) -> f32 {
        self.0.total_volume()
    }
}

#[pyclass]
struct IsometryTransform(nalgebra::IsometryMatrix3<f32>);

// Free functions

#[pyfunction(name = "fit_transform")]
fn fit_transform_py(sel1: &Sel, sel2: &Sel) -> anyhow::Result<IsometryTransform> {
    let tr = molar::prelude::fit_transform(&sel1.0, &sel2.0)?;
    Ok(IsometryTransform(tr))
}

#[pyfunction(name = "fit_transform_matching")]
fn fit_transform_matching_py(sel1: &Sel, sel2: &Sel) -> anyhow::Result<IsometryTransform> {
    let tr = molar::core::fit_transform_matching(&sel1.0, &sel2.0)?;
    Ok(IsometryTransform(tr))
}

#[pyfunction]
fn rmsd(sel1: &Sel, sel2: &Sel) -> anyhow::Result<f32> {
    Ok(molar::core::Sel::rmsd(&sel1.0, &sel2.0)?)
}

#[pyfunction(name = "rmsd_mw")]
fn rmsd_mw_py(sel1: &Sel, sel2: &Sel) -> anyhow::Result<f32> {
    Ok(molar::core::rmsd_mw(&sel1.0, &sel2.0)?)
}

#[pyclass]
struct ParticleIterator {
    sel: Py<Sel>,
    cur: isize,
}

#[pymethods]
impl ParticleIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PyObject> {
        let ret = Python::with_gil(|py| {
            let s = slf.sel.bind(py);
            Sel::__getitem__(s.clone(), slf.cur)
        })
        .ok();
        slf.cur += 1;
        ret
    }
}

#[pyfunction]
#[pyo3(signature = (cutoff,data1,data2=None,dims=[false,false,false]))]
fn distance_search<'py>(
    py: Python<'py>,
    cutoff: &Bound<'py, PyAny>,
    data1: &Bound<'py, Sel>,
    data2: Option<&Bound<'py, Sel>>,
    dims: [bool; 3],
) -> anyhow::Result<Bound<'py, PyAny>> {
    let mut res: Vec<(usize, usize, f32)>;
    let pbc_dims = PbcDims::new(dims[0], dims[1], dims[2]);
    let sel1 = data1.borrow();

    if let Ok(d) = cutoff.extract::<f32>() {
        // Distance cutoff
        if let Some(d2) = data2 {
            let sel2 = d2.borrow();
            if pbc_dims.any() {
                res = molar::core::distance_search_double_pbc(
                    d,
                    sel1.0.iter_pos(),
                    sel2.0.iter_pos(),
                    sel1.0.iter_index(),
                    sel2.0.iter_index(),
                    sel1.0.get_box().ok_or_else(|| anyhow!("no periodic box"))?,
                    pbc_dims,
                );
            } else {
                res = molar::core::distance_search_double(
                    d,
                    sel1.0.iter_pos(),
                    sel2.0.iter_pos(),
                    sel1.0.iter_index(),
                    sel2.0.iter_index(),
                );
            }
        } else {
            if pbc_dims.any() {
                res = molar::core::distance_search_single_pbc(
                    d,
                    sel1.0.iter_pos(),
                    sel1.0.iter_index(),
                    sel1.0.get_box().ok_or_else(|| anyhow!("no periodic box"))?,
                    pbc_dims,
                );
            } else {
                res =
                    molar::core::distance_search_single(d, sel1.0.iter_pos(), sel1.0.iter_index());
            }
        }
    } else if let Ok(s) = cutoff.extract::<String>() {
        if s != "vdw" {
            bail!("Unknown cutoff type {s}");
        }

        // VdW cutof
        let vdw1: Vec<f32> = sel1.0.iter_atoms().map(|a| a.vdw()).collect();

        if sel1.0.len() != vdw1.len() {
            bail!("Size mismatch 1: {} {}", sel1.0.len(), vdw1.len());
        }

        if let Some(d2) = data2 {
            let sel2 = d2.borrow();
            let vdw2: Vec<f32> = sel2.0.iter_atoms().map(|a| a.vdw()).collect();

            if sel2.0.len() != vdw2.len() {
                bail!("Size mismatch 2: {} {}", sel2.0.len(), vdw2.len());
            }

            if pbc_dims.any() {
                res = molar::core::distance_search_double_vdw(
                    sel1.0.iter_pos(),
                    sel2.0.iter_pos(),
                    &vdw1,
                    &vdw2,
                );
            } else {
                res = molar::core::distance_search_double_vdw_pbc(
                    sel1.0.iter_pos(),
                    sel2.0.iter_pos(),
                    &vdw1,
                    &vdw2,
                    sel1.0.get_box().ok_or_else(|| anyhow!("no periodic box"))?,
                    pbc_dims,
                );
            }

            // Convert local indices to global
            unsafe {
                for el in &mut res {
                    el.0 = sel1.0.get_index_unchecked(el.0);
                    el.1 = sel2.0.get_index_unchecked(el.1);
                }
            }
        } else {
            bail!("VdW distance search is not yet supported for single selection");
        }
    } else {
        unreachable!()
    };

    // Subdivide the result into two arrays
    unsafe {
        // Pairs array
        let pairs_arr = numpy::PyArray2::<usize>::new(py, [res.len(), 2], true);
        for i in 0..res.len() {
            pairs_arr.uget_raw([i, 0]).write(res[i].0);
            pairs_arr.uget_raw([i, 1]).write(res[i].1);
        }

        // Distances array
        let dist_arr = numpy::PyArray1::<f32>::new(py, [res.len()], true);
        for i in 0..res.len() {
            dist_arr.uget_raw(i).write(res[i].2);
        }

        Ok((pairs_arr, dist_arr).into_bound_py_any(py)?)
    }
}

#[pyclass]
struct NdxFile(molar::core::NdxFile);

#[pymethods]
impl NdxFile {
    #[new]
    fn new(fname: &str) -> anyhow::Result<Self> {
        Ok(NdxFile(molar::core::NdxFile::new(fname)?))
    }

    fn get_group_as_sel(&self, gr_name: &str, src: &System) -> anyhow::Result<Sel> {
        Ok(Sel::new_owned(self.0.get_group_as_sel(gr_name, &src.0)?))
    }
}

//====================================
#[pyfunction]
fn greeting() {
    molar::greeting("molar_python");
}

/// A Python module implemented in Rust.
#[pymodule(name = "molar")]
//#[pymodule]
fn molar_python(m: &Bound<'_, PyModule>) -> PyResult<()> {
    pyo3_log::init();
    m.add_class::<Atom>()?;
    m.add_class::<Particle>()?;
    m.add_class::<Topology>()?;
    m.add_class::<State>()?;
    m.add_class::<PeriodicBox>()?;
    m.add_class::<FileHandler>()?;
    m.add_class::<System>()?;
    m.add_class::<Sel>()?;
    m.add_class::<SasaResults>()?;
    m.add_class::<NdxFile>()?;
    m.add_class::<Histogram1D>()?;
    m.add_function(wrap_pyfunction!(greeting, m)?)?;
    m.add_function(wrap_pyfunction!(fit_transform_py, m)?)?;
    m.add_function(wrap_pyfunction!(fit_transform_matching_py, m)?)?;
    m.add_function(wrap_pyfunction!(rmsd, m)?)?;
    m.add_function(wrap_pyfunction!(rmsd_mw_py, m)?)?;
    m.add_function(wrap_pyfunction!(distance_search, m)?)?;
    m.add_class::<LipidMolecule>()?;
    m.add_class::<Membrane>()?;
    Ok(())
}
