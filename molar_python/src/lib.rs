use std::path::Path;

use anyhow::{anyhow, bail};
use molar::prelude::*;
use numpy::{
    nalgebra::{self, Const, Dyn, VectorView},
    PyArrayLike1, PyArrayMethods, PyReadonlyArray2, PyUntypedArrayMethods, ToPyArray,
};
use pyo3::{
    prelude::*,
    types::PyTuple,
    IntoPyObjectExt,
};

mod utils;
use utils::*;

mod atom;
use atom::Atom;

mod particle;
use particle::Particle;

mod periodic_box;
use periodic_box::PeriodicBox;

//-------------------------------------------

#[pyclass(unsendable)]
struct Topology(molar::core::Holder<molar::core::Topology, BuilderSerial>);

#[pyclass(unsendable)]
struct State(molar::core::Holder<molar::core::State, BuilderSerial>);

#[pymethods]
impl State {
    fn __len__(&self) -> usize {
        self.0.len()
    }

    #[getter]
    fn time(&self) -> f32 {
        self.0.get_time()
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
struct FileHandler(molar::io::FileHandler);

#[pymethods]
impl FileHandler {
    #[new]
    fn new(fname: &str) -> anyhow::Result<Self> {
        Ok(FileHandler(molar::io::FileHandler::open(fname)?))
    }

    fn read(&mut self) -> anyhow::Result<(Topology, State)> {
        let (top, st) = self.0.read()?;
        Ok((Topology(top.into()), State(st.into())))
    }

    fn read_topology(&mut self) -> anyhow::Result<Topology> {
        let top = self.0.read_topology()?;
        Ok(Topology(top.into()))
    }

    fn read_state(&mut self) -> anyhow::Result<State> {
        if let Some(st) = self.0.read_state()? {
            Ok(State(st.into()))
        } else {
            Err(anyhow!("can't read state"))
        }
    }

    fn write(&mut self, data: Bound<'_, PyAny>) -> anyhow::Result<()> {
        if let Ok(s) = data.extract::<PyRef<'_, Source>>() {
            self.0.write(&s.0)?;
        } else if let Ok(s) = data.extract::<PyRef<'_, Sel>>() {
            self.0.write(&s.0)?;
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
            self.0.write(&(top.0.clone(), st.0.clone()))?;
        } else {
            return Err(anyhow!(
                "Invalid data type {} when writing to file",
                data.get_type()
            ));
        }
        Ok(())
    }

    fn write_topology(&mut self, data: Bound<'_, PyAny>) -> anyhow::Result<()> {
        if let Ok(s) = data.extract::<PyRef<'_, Source>>() {
            self.0.write_topology(&s.0)?;
        } else if let Ok(s) = data.extract::<PyRef<'_, Sel>>() {
            self.0.write_topology(&s.0)?;
        } else if let Ok(s) = data.extract::<PyRefMut<'_, Topology>>() {
            self.0.write_topology(&s.0)?;
        } else {
            return Err(anyhow!(
                "Invalid data type {} when writing to file",
                data.get_type()
            ));
        }
        Ok(())
    }

    fn write_state(&mut self, data: Bound<'_, PyAny>) -> anyhow::Result<()> {
        if let Ok(s) = data.extract::<PyRef<'_, Source>>() {
            self.0.write_state(&s.0)?;
        } else if let Ok(s) = data.extract::<PyRef<'_, Sel>>() {
            self.0.write_state(&s.0)?;
        } else if let Ok(s) = data.extract::<PyRefMut<'_, State>>() {
            self.0.write_state(&s.0)?;
        } else {
            return Err(anyhow!(
                "Invalid data type {} when writing to file",
                data.get_type()
            ));
        }
        Ok(())
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PyObject> {
        let st = slf.read_state();
        if st.is_ok() {
            Python::with_gil(|py| Some(st.unwrap().into_py_any(py)))
                .unwrap()
                .ok()
        } else {
            None
        }
    }

    fn skip_to_frame(&mut self, fr: usize) -> PyResult<()> {
        self.0.skip_to_frame(fr).map_err(|e| anyhow!(e))?;
        Ok(())
    }

    fn skip_to_time(&mut self, t: f32) -> anyhow::Result<()> {
        self.0.skip_to_time(t)?;
        Ok(())
    }

    fn tell_first(&self) -> anyhow::Result<(usize, f32)> {
        Ok(self.0.tell_first()?)
    }

    fn tell_current(&self) -> anyhow::Result<(usize, f32)> {
        Ok(self.0.tell_current()?)
    }

    fn tell_last(&self) -> anyhow::Result<(usize, f32)> {
        Ok(self.0.tell_last()?)
    }

    #[getter]
    fn stats(&self) -> FileStats {
        FileStats(self.0.stats.clone())
    }

    #[getter]
    fn file_name(&self) -> &Path {
        &self.0.file_name
    }
}

#[pyclass(unsendable)]
struct _ParTrajReader {
    iter: molar::io::IoStateIterator,
}

#[pymethods]
impl _ParTrajReader {
    #[new]
    fn new(fname: &str) -> anyhow::Result<Self> {
        let iter = molar::io::FileHandler::open(fname)?.into_iter();
        Ok(Self { iter })
    }

    fn next_state(&mut self) -> Option<State> {
        self.iter.next().map(|st| State(st.into()))
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
struct Source(molar::core::Source<molar::core::BuilderSerial>);

#[pymethods]
impl Source {
    #[new]
    #[pyo3(signature = (*py_args))]
    fn new<'py>(py_args: &Bound<'py, PyTuple>) -> PyResult<Self> {
        if py_args.len() == 1 {
            // From file
            Ok(Source(
                molar::core::Source::builder_from_file(&py_args.get_item(0)?.extract::<String>()?)
                    .map_err(|e| anyhow!(e))?,
            ))
        } else if py_args.len() == 2 {
            let top = py_args
                .get_item(0)?
                .downcast::<Topology>()?
                .try_borrow_mut()?;
            let st = py_args.get_item(1)?.downcast::<State>()?.try_borrow_mut()?;
            Ok(Source(
                molar::core::Source::new_builder(top.0.clone(), st.0.clone())
                    .map_err(|e| anyhow!(e))?,
            ))
        } else {
            Err(anyhow!("wrong number of arguments: 1 or 2 reqired")).map_err(|e| anyhow!(e))?
        }
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn select_all(&mut self) -> PyResult<Sel> {
        Ok(Sel(self.0.select_all().map_err(|e| anyhow!(e))?))
    }

    fn select(&mut self, sel_str: &str) -> PyResult<Sel> {
        Ok(Sel(self.0.select(sel_str).map_err(|e| anyhow!(e))?))
    }

    #[pyo3(signature = (arg=None))]
    fn __call__(&self, arg: Option<&Bound<'_, PyAny>>) -> PyResult<Sel> {
        if let Some(arg) = arg {
            if let Ok(val) = arg.extract::<String>() {
                if val.is_empty() {
                    Ok(Sel(self.0.select_all().map_err(|e| anyhow!(e))?))
                } else {
                    Ok(Sel(self.0.select(val).map_err(|e| anyhow!(e))?))
                }
            } else if let Ok(val) = arg.extract::<(usize, usize)>() {
                Ok(Sel(self
                    .0
                    .select(val.0..val.1)
                    .map_err(|e| anyhow!(e))?))
            } else if let Ok(val) = arg.extract::<Vec<usize>>() {
                Ok(Sel(self.0.select(val).map_err(|e| anyhow!(e))?))
            } else {
                Err(anyhow!(
                    "Invalid argument type {} when creating selection",
                    arg.get_type()
                )
                .into())
            }
        } else {
            Ok(Sel(self.0.select_all().map_err(|e| anyhow!(e))?))
        }
    }

    fn set_state<'py>(&mut self, st: &Bound<'py, State>) -> PyResult<Bound<'py, State>> {
        let mut st_ref = st.borrow_mut();
        // In Python we can pass by value, so we have to release State from the
        // Python object. To do this firt swap it with new dummy Holder
        // which is uniquilly owned and then release it from this holder
        let mut dum_holder = Holder::new(molar::core::State::default());
        unsafe { dum_holder.swap_unchecked(&mut st_ref.0) }; // st_ref is empty at this point
                                                             // dum_holder is uniquelly owned, so this never fails
        let dum_st = dum_holder.release().unwrap();
        // Now call set_state as usual
        let old_st = self.0.set_state(dum_st).map_err(|e| anyhow!(e))?;
        // We should not leave st empty, it should point to the same state as self
        unsafe { st_ref.0.replace_arc(self.0.get_state()) };
        // Pack old_st and return it
        Bound::new(st.py(), State(old_st.into()))
    }

    fn get_state(&self) -> State {
        State(self.0.get_state())
    }

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
        } else {
            unreachable!()
        }
    }

    fn append(&self, arg: &Bound<'_, PyAny>) -> anyhow::Result<()> {
        // In the future other types can be used as well
        if let Ok(sel) = arg.downcast::<Sel>() {
            self.0.append(&sel.borrow().0);
        } else if let Ok(sel) = arg.downcast::<Source>() {
            self.0.append(&sel.borrow().0);
        } else {
            anyhow::bail!("Unsupported type to append a Source")
        }
        Ok(())
    }
}

//====================================

#[pyclass(sequence, unsendable)]
struct Sel(molar::core::Sel<BuilderSerial>);

#[pymethods]
impl Sel {
    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn __call__(&self, arg: &Bound<'_, PyAny>) -> PyResult<Sel> {
        if let Ok(val) = arg.extract::<String>() {
            Ok(Sel(self.0.subsel(val).map_err(|e| anyhow!(e))?))
        } else if let Ok(val) = arg.extract::<(usize, usize)>() {
            Ok(Sel(self
                .0
                .subsel(val.0..=val.1)
                .map_err(|e| anyhow!(e))?))
        } else if let Ok(val) = arg.extract::<Vec<usize>>() {
            Ok(Sel(self
                .0
                .subsel(val)
                .map_err(|e| anyhow!(e))?))
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
        let p = s.0.nth_particle_mut(ind).unwrap();
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

    fn get_coord<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2::<f32>> {
        // Instead we allocate an uninitialized PyArray manually and fill it by data.
        // By doing this we save on unnecessary initiallization and extra allocation
        unsafe {
            let arr = numpy::PyArray2::<f32>::new(py,[3,self.0.len()],true);
            let arr_ptr = arr.data();
            for i in 0..self.0.len() {
                let pos_ptr = self.0.nth_pos_unchecked(i).coords.as_ptr();
                // This is faster than copying by element with uget_raw()
                std::ptr::copy_nonoverlapping(pos_ptr,arr_ptr.offset(i as isize * 3), 3);
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
                let pos_ptr = self.0.nth_pos_mut_unchecked(i).coords.as_mut_ptr();
                std::ptr::copy_nonoverlapping(ptr.offset(i as isize * 3), pos_ptr, 3);
            }
        }

        Ok(())
    }

    fn set_state(&mut self, st: Bound<'_, State>) -> PyResult<()> {
        let _ = self
            .0
            .set_state(st.try_borrow_mut()?.0.clone())
            .map_err(|e| anyhow!(e));
        Ok(())
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

    fn split_resindex(&self) -> anyhow::Result<Vec<Sel>> {
        Ok(self.0.split_resindex_iter()?.map(|sel| Sel(sel)).collect())
    }

    fn split_chain(&self) -> anyhow::Result<Vec<Sel>> {
        Ok(self.0.split_iter(|p| Some(p.atom.chain))?.map(|sel| Sel(sel)).collect())
    }

    fn split_molecule(&self) -> anyhow::Result<Vec<Sel>> {
        Ok(self.0.split_molecules_iter()?.map(|sel| Sel(sel)).collect())
    }

    fn to_gromacs_ndx(&self, name: &str) -> String {
        self.0.to_gromacs_ndx(name)
    }

    fn add(&mut self, arg: &Bound<'_, PyAny>) -> anyhow::Result<()> {
        // In the future other types can be used as well
        if let Ok(sel) = arg.downcast::<Sel>() {
            self.0.add(&sel.borrow().0)?;
        } else {
            anyhow::bail!("Unsupported type to append a Sel")
        }
        Ok(())
    }

    /// += opeator (append in place)
    fn __iadd__(&mut self, arg: &Bound<'_, PyAny>) -> anyhow::Result<()> {
        self.add(arg)
    }

    /// Invert in place
    fn invert(&mut self) {
        self.0.invert();
    }

    /// remove other from self
    fn remove_global(&mut self, other: &Sel) -> anyhow::Result<()> {
        self.0.remove_global(&other.0)?;
        Ok(())
    }

    /// -= (remove other from self in place)
    fn __isub__(&mut self, other: &Sel) -> anyhow::Result<()> {
        self.remove_global(other)
    }

    /// operator |
    fn __or__(&self, rhs: &Sel) -> Sel {
        Sel(self.0.union(&rhs.0))
    }

    /// operator &
    fn __and__(&self, rhs: &Sel) -> anyhow::Result<Sel> {
        Ok(Sel(self.0.intersection(&rhs.0)?))
    }

    /// -= (remove other from self)
    fn __sub__(&self, rhs: &Sel) -> anyhow::Result<Sel> {
        Ok(Sel(self.0.difference(&rhs.0)?))
    }

    /// ~ operator
    fn __invert__(&self) -> anyhow::Result<Sel> {
        Ok(Sel(self.0.complement()?))
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

#[pyfunction]
fn fit_transform(sel1: &Sel, sel2: &Sel) -> anyhow::Result<IsometryTransform> {
    let tr = molar::core::Sel::fit_transform(&sel1.0, &sel2.0)?;
    Ok(IsometryTransform(tr))
}

#[pyfunction]
fn rmsd(sel1: &Sel, sel2: &Sel) -> anyhow::Result<f32> {
    Ok(molar::core::Sel::rmsd(&sel1.0, &sel2.0)?)
}

#[pyfunction]
fn rmsd_mw(sel1: &Sel, sel2: &Sel) -> anyhow::Result<f32> {
    Ok(molar::core::Sel::rmsd_mw(&sel1.0, &sel2.0)?)
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
    let sel1 = &data1.borrow().0;    

    if let Ok(d) = cutoff.extract::<f32>() {
        // Distance cutoff
        if let Some(d2) = data2 {
            let sel2 = &d2.borrow().0;
            if pbc_dims.any() {
                res = molar::core::distance_search_double_pbc(
                    d,
                    sel1.iter_pos(),
                    sel2.iter_pos(),
                    sel1.iter_index(),
                    sel2.iter_index(),
                    sel1.get_box().ok_or_else(|| anyhow!("no periodic box"))?,
                    pbc_dims,
                );
            } else {
                res = molar::core::distance_search_double(
                    d,
                    sel1.iter_pos(),
                    sel2.iter_pos(),
                    sel1.iter_index(),
                    sel2.iter_index(),
                );
            }
        } else {
            if pbc_dims.any() {
                res = molar::core::distance_search_single_pbc(
                    d,
                    sel1.iter_pos(),
                    sel1.iter_index(),                    
                    sel1.get_box().ok_or_else(|| anyhow!("no periodic box"))?,
                    pbc_dims,
                );
            } else {
                res = molar::core::distance_search_single(
                    d,
                    sel1.iter_pos(),
                    sel1.iter_index(),                    
                );
            }
        }

        
    } else if let Ok(s) = cutoff.extract::<String>() {
        if s != "vdw" {
            bail!("Unknown cutoff type {s}");
        }

        // VdW cutof
        let vdw1: Vec<f32> = sel1.iter_atoms().map(|a| a.vdw()).collect();

        if sel1.len() != vdw1.len() {
            bail!("Size mismatch 1: {} {}",sel1.len(),vdw1.len());
        }

        if let Some(d2) = data2 {
            let sel2 = &d2.borrow().0;
            let vdw2: Vec<f32> = sel2.iter_atoms().map(|a| a.vdw()).collect();

            if sel2.len() != vdw2.len() {
                bail!("Size mismatch 2: {} {}",sel2.len(),vdw2.len());
            }

            if pbc_dims.any() {
                res = molar::core::distance_search_double_vdw(
                    sel1.iter_pos(),
                    sel2.iter_pos(),
                    &vdw1,
                    &vdw2,
                );
            } else {
                res = molar::core::distance_search_double_vdw_pbc(
                    sel1.iter_pos(),
                    sel2.iter_pos(),
                    &vdw1,
                    &vdw2,
                    sel1.get_box().ok_or_else(|| anyhow!("no periodic box"))?,
                    pbc_dims,
                );
            }
            
            // Convert local indices to global
            for el in &mut res {
                el.0 = sel1.nth_index(el.0).unwrap();
                el.1 = sel2.nth_index(el.1).unwrap();
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
        let pairs_arr = numpy::PyArray2::<usize>::new(py,[res.len(),2],true);
        for i in 0..res.len() {
            pairs_arr.uget_raw([i,0]).write(res[i].0);
            pairs_arr.uget_raw([i,1]).write(res[i].1);
        }
        
        // Distances array
        let dist_arr = numpy::PyArray1::<f32>::new(py,[res.len()],true);
        for i in 0..res.len() {
            dist_arr.uget_raw(i).write(res[i].2);
        }

        Ok((pairs_arr,dist_arr).into_bound_py_any(py)?)
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
    m.add_class::<Source>()?;
    m.add_class::<Sel>()?;
    m.add_class::<SasaResults>()?;
    m.add_class::<_ParTrajReader>()?;
    m.add_function(wrap_pyfunction!(greeting, m)?)?;
    m.add_function(wrap_pyfunction!(fit_transform, m)?)?;
    m.add_function(wrap_pyfunction!(rmsd, m)?)?;
    m.add_function(wrap_pyfunction!(rmsd_mw, m)?)?;
    m.add_function(wrap_pyfunction!(distance_search, m)?)?;
    Ok(())
}
