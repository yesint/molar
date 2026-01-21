use std::path::PathBuf;

use anyhow::{anyhow, bail};
use molar::prelude::*;
//use molar::prelude::*;
use numpy::{
    nalgebra::{self, Const, Dyn, VectorView},
    PyArray1, PyArrayLike1, PyArrayMethods, PyReadonlyArray2, PyUntypedArrayMethods, ToPyArray,
};
use pyo3::{
    exceptions::PyTypeError,
    prelude::*,
    types::{PySlice, PyTuple},
    IntoPyObjectExt,
};

mod utils;
use triomphe::Arc;
use utils::*;

mod atom;
use atom::AtomPy;

mod particle;
use particle::ParticlePy;

mod periodic_box;
use periodic_box::PeriodicBoxPy;

// mod membrane;
// use membrane::*;

mod topology_state;
use topology_state::*;
//-------------------------------------------

#[pyclass(unsendable, name = "FileHandler")]
struct FileHandlerPy(Option<FileHandler>, Option<IoStateIterator>);

const ALREADY_TRANDFORMED: &str = "file handler is already transformed to state iterator";

#[pymethods]
impl FileHandlerPy {
    #[new]
    fn new(fname: &str, mode: &str) -> anyhow::Result<Self> {
        match mode {
            "r" => Ok(FileHandlerPy(Some(FileHandler::open(fname)?), None)),
            "w" => Ok(FileHandlerPy(Some(FileHandler::create(fname)?), None)),
            _ => Err(anyhow!("Wrong file open mode")),
        }
    }

    fn read(&mut self) -> anyhow::Result<(TopologyPy, StatePy)> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        let (top, st) = h.read()?;

        Ok((top.into(), st.into()))
    }

    fn read_topology(&mut self) -> anyhow::Result<TopologyPy> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        let top = h.read_topology()?;
        Ok(top.into())
    }

    fn read_state(&mut self) -> anyhow::Result<StatePy> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        let st = h.read_state()?;
        Ok(st.into())
    }

    fn write(&mut self, data: Bound<'_, PyAny>) -> anyhow::Result<()> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        if let Ok(s) = data.extract::<PyRef<'_, SystemPy>>() {
            h.write(&*s)?;
        } else if let Ok(s) = data.extract::<PyRef<'_, SelPy>>() {
            h.write(&*s)?;
        } else if let Ok(s) = data.cast::<PyTuple>() {
            if s.len() != 2 {
                return Err(anyhow!("tuple must have two elements, not {}", s.len()));
            }
            let top = s
                .get_item(0)?
                .cast::<TopologyPy>()
                .map_err(|_| anyhow!("first tuple element should be a Topology"))?
                .as_ptr() as *const TopologyPy;

            let st = s
                .get_item(1)?
                .cast::<StatePy>()
                .map_err(|_| anyhow!("second tuple element should be a State"))?
                .as_ptr() as *const StatePy;
            h.write_topology(unsafe { &*top })?;
            h.write_state(unsafe { &*st })?;
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
        if let Ok(s) = data.cast::<SystemPy>() {
            h.write_topology(&s.borrow().top)?;
        } else if let Ok(s) = data.cast::<SelPy>() {
            h.write_topology(&s.borrow().sys.top)?;
        } else if let Ok(s) = data.cast::<TopologyPy>() {
            h.write_topology(&*s.borrow())?;
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
        if let Ok(s) = data.cast::<SystemPy>() {
            h.write_state(&s.borrow().st)?;
        } else if let Ok(s) = data.cast::<SelPy>() {
            h.write_state(&s.borrow().sys.st)?;
        } else if let Ok(s) = data.cast::<StatePy>() {
            h.write_state(&*s.borrow())?;
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

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<Py<PyAny>> {
        let st = slf
            .1
            .as_mut()
            .expect("Not transformed to state iterator yet")
            .next()
            .map(|st| StatePy::from(st));
        if st.is_some() {
            Some(st.unwrap().into_py_any(slf.py()).unwrap())
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
    fn stats(&self) -> anyhow::Result<FileStatsPy> {
        let h = self
            .0
            .as_ref()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        Ok(FileStatsPy(h.stats.clone()))
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

#[pyclass(name = "FileStats")]
struct FileStatsPy(FileStats);

#[pymethods]
impl FileStatsPy {
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

#[pyclass(name = "System")]
struct SystemPy {
    top: TopologyPy,
    st: StatePy,
}

impl SystemPy {
    pub(crate) fn clone_ref(&self) -> SystemPy {
        Self {
            top: self.top.clone_ref(),
            st: self.st.clone_ref(),
        }
    }
}

impl LenProvider for SystemPy {
    fn len(&self) -> usize {
        self.top.len()
    }
}

impl IndexProvider for SystemPy {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        i
    }
}

impl AtomPosAnalysis for SystemPy {
    fn atoms_ptr(&self) -> *const Atom {
        self.top.get().atoms.as_ptr()
    }

    fn coords_ptr(&self) -> *const Pos {
        self.st.get().coords.as_ptr()
    }
}

impl AtomPosAnalysisMut for SystemPy {}

impl BoxProvider for SystemPy {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.st.get().get_box()
    }
}

impl RandomBondProvider for SystemPy {
    fn num_bonds(&self) -> usize {
        0
    }

    unsafe fn get_bond_unchecked(&self, _i: usize) -> &[usize; 2] {
        unreachable!()
    }
}

impl TimeProvider for SystemPy {
    fn get_time(&self) -> f32 {
        self.st.get().get_time()
    }
}

impl SaveTopology for SystemPy {
    fn iter_atoms_dyn<'a>(&'a self) -> Box<dyn Iterator<Item = &'a Atom> + 'a> {
        Box::new(self.iter_atoms())
    }
}
impl SaveState for SystemPy {
    fn iter_pos_dyn<'a>(&'a self) -> Box<dyn Iterator<Item = &'a Pos> + 'a> {
        Box::new(self.iter_pos())
    }
}
impl SaveTopologyState for SystemPy {}

#[pymethods]
impl SystemPy {
    #[new]
    #[pyo3(signature = (*py_args))]
    fn new(py_args: &Bound<'_, PyTuple>) -> PyResult<Self> {
        if py_args.len() == 1 {
            // From file
            let fname = py_args.get_item(0)?.extract::<String>()?;
            let mut fh = molar::io::FileHandler::open(&fname).map_err(|e| anyhow!(e))?;
            let (top, st) = fh.read().map_err(|e| anyhow!(e))?;
            Ok(SystemPy {
                top: top.into(),
                st: st.into(),
            })
        } else if py_args.len() == 2 {
            // From existing Topology and State python objects
            let arg1 = py_args.get_item(0)?;
            let arg2 = py_args.get_item(1)?;
            let top = arg1.cast::<TopologyPy>()?;
            let st = arg2.cast::<StatePy>()?;
            let n1 = top.borrow().len();
            let n2 = st.borrow().len();
            if n1 != n2 {
                Err(PyTypeError::new_err(
                    "topology and state are of different size",
                ))
            } else {
                Ok(SystemPy {
                    top: top.borrow().clone_ref(),
                    st: st.borrow().clone_ref(),
                })
            }
        } else {
            // Empty System
            Ok(SystemPy {
                top: TopologyPy(Default::default()),
                st: StatePy(Default::default()),
            })
        }
    }

    fn __len__(&self) -> usize {
        self.top.get().len()
    }

    // fn select_all(slf: &Bound<Self>) -> anyhow::Result<Sel> {
    //     let bs = slf.borrow();
    //     let top = bs.top.borrow(slf.py());
    //     let st = bs.st.borrow(slf.py());
    //     let v = (0..top.0.len()).into_sel_index(&top.0, &st.0, None)?;
    //     Ok(Sel {
    //         top: Py::clone_ref(&bs.top, slf.py()),
    //         st: Py::clone_ref(&bs.st, slf.py()),
    //         sel: Sel::from_svec(v).unwrap(),
    //     })
    // }

    // fn select(slf: &Bound<Self>, sel_str: &str) -> anyhow::Result<Sel> {
    //     let bs = slf.borrow();
    //     let top = bs.top.borrow(slf.py());
    //     let st = bs.st.borrow(slf.py());
    //     let v = sel_str.into_sel_index(&top.0, &st.0, None)?;
    //     Ok(Sel {
    //         top: Py::clone_ref(&bs.top, slf.py()),
    //         st: Py::clone_ref(&bs.st, slf.py()),
    //         sel: Sel::from_svec(v).unwrap(),
    //     })
    // }

    #[pyo3(signature = (arg=None))]
    fn __call__(slf: &Bound<Self>, arg: Option<&Bound<'_, PyAny>>) -> anyhow::Result<SelPy> {
        let sys = slf.borrow();

        let index = if let Some(arg) = arg {
            // Argument present
            if let Ok(val) = arg.extract::<String>() {
                if val.is_empty() {
                    // Select all on empty string
                    (0..sys.len()).into_sel_index(&*sys, None)?
                } else {
                    // Otherwise do normal textual selection
                    val.into_sel_index(&*sys, None)?
                }
            } else if let Ok(val) = arg.extract::<(usize, usize)>() {
                // Range selection
                (val.0..val.1).into_sel_index(&*sys, None)?
            } else if let Ok(val) = arg.extract::<Vec<usize>>() {
                // Vector of indices
                val.into_sel_index(&*sys, None)?
            } else {
                bail!(
                    "Invalid argument type {} when creating selection",
                    arg.get_type()
                )
            }
        } else {
            // No argument, select all
            (0..sys.len()).into_sel_index(&*sys, None)?
        };

        Ok(SelPy {
            sys: slf.borrow().clone_ref(),
            index,
        })
    }

    fn set_state(&mut self, st: &Bound<'_, StatePy>) -> anyhow::Result<StatePy> {
        if self.st.get().interchangeable(st.borrow().get()) {
            let ret = self.st.clone_ref();
            self.st = st.borrow().clone_ref();
            Ok(ret)
        } else {
            bail!("incompatible state")
        }
    }

    fn set_topology(&mut self, top: &Bound<'_, TopologyPy>) -> anyhow::Result<TopologyPy> {
        if self.top.get().interchangeable(top.borrow().get()) {
            let ret = self.top.clone_ref();
            self.top = top.borrow().clone_ref();
            Ok(ret)
        } else {
            bail!("incompatible topology")
        }
    }

    fn save(&self, fname: &str) -> anyhow::Result<()> {
        Ok(SaveTopologyState::save(self, fname)?)
    }

    fn remove<'py>(slf: &Bound<'py, Self>, arg: &Bound<'py, PyAny>) -> anyhow::Result<()> {
        if let Ok(sel) = arg.cast::<SelPy>() {
            // Selection provided
            let sb = sel.borrow();
            sb.sys.top.get_mut().remove_atoms(sb.iter_index())?;
            sb.sys.st.get_mut().remove_coords(sb.iter_index())?;
            Ok(())
        } else {
            let sel = Self::__call__(slf, Some(arg))?;
            sel.sys.top.get_mut().remove_atoms(sel.iter_index())?;
            sel.sys.st.get_mut().remove_coords(sel.iter_index())?;
            Ok(())
        }
    }

    #[pyo3(signature = (*args))]
    fn append<'py>(slf: &Bound<'py, Self>, args: &Bound<'py, PyTuple>) -> anyhow::Result<()> {
        let slf_b = slf.borrow();

        if args.len() == 1 {
            let arg = args.get_item(0)?;
            let sel = if let Ok(sel) = arg.cast::<SelPy>() {
                sel
            } else {
                &Bound::new(slf.py(), Self::__call__(slf, Some(&arg))?)?
            };

            slf_b.top.get_mut().add_atoms(sel.borrow().iter_atoms().cloned());
            slf_b.st.get_mut().add_coords(sel.borrow().iter_pos().cloned());

            Ok(())
        } else if args.len() == 2 {
            let arg1 = args.get_item(0)?;
            let ab = arg1
                .cast::<AtomPy>()
                .map_err(|_| anyhow!("expected atom"))?
                .borrow();
            let pos = args
                .get_item(1)?
                .extract::<PyArrayLike1<f32>>()
                .map_err(|_| anyhow!("expected pos"))?;
            let v: VectorView<f32, Const<3>> = pos.try_as_matrix().unwrap();
            slf_b.top.get_mut().add_atoms(std::iter::once(&ab.0).cloned());
            slf_b.st.get_mut().add_coords(std::iter::once(Pos::new(v.x, v.y, v.z)));
            Ok(())
        } else {
            bail!("1 or 2 arguments expected");
        }
    }

    // fn append_atom_pos(slf: &Bound<'_, Self>, at: &AtomPy, pos: PyArrayLike1<'_, f32>) {
    //     let v: VectorView<f32, Const<3>> = pos.try_as_matrix().unwrap();
    //     slf.borrow_mut()
    //         .top
    //         .borrow_mut(slf.py())
    //         .0
    //         .add_atoms(std::iter::once(&at.0).cloned());
    //     slf.borrow_mut()
    //         .st
    //         .borrow_mut(slf.py())
    //         .0
    //         .add_coords(std::iter::once(Pos::new(v.x, v.y, v.z)));
    // }

    #[getter]
    fn get_time(&self) -> f32 {
        TimeProvider::get_time(self)
    }

    #[getter]
    fn get_box(&self, py: Python<'_>) -> PeriodicBoxPy {
        self.st
            .borrow(py)
            .pbox
            .as_ref()
            .map(|b| PeriodicBoxPy(b.clone()))
            .unwrap()
    }

    #[setter]
    fn set_box(&mut self, py: Python<'_>, b: &PeriodicBoxPy) {
        *self.st.borrow_mut(py).pbox.as_mut().unwrap() = b.0.clone();
    }

    #[setter]
    fn set_time(&mut self, py: Python<'_>, t: f32) {
        self.st.borrow_mut(py).time = t;
    }

    fn set_box_from(&mut self, sys: Bound<'_, SystemPy>) {
        let sys_ref = sys.borrow();
        let py = sys.py();
        self.st.bind(py).borrow_mut().pbox = sys_ref.st.bind(py).borrow().pbox.clone();
    }
}

//====================================

#[pyclass(sequence, name = "Sel")]
struct SelPy {
    sys: SystemPy,
    index: SVec,
}

impl LenProvider for SelPy {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SelPy {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        self.index.get_index_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> {
        self.index.iter_index()
    }
}

impl AtomPosAnalysis for SelPy {
    fn atoms_ptr(&self) -> *const Atom {
        Python::attach(|py| self.top.borrow(py).0.atoms.as_ptr())
    }

    fn coords_ptr(&self) -> *const Pos {
        Python::attach(|py| self.st.borrow(py).coords.bind(py).data() as *const Pos)
    }
}

impl AtomPosAnalysisMut for SelPy {}

impl BoxProvider for SelPy {
    fn get_box(&self) -> Option<&PeriodicBox> {
        Python::attach(|py| {
            let ptr = self
                .sys
                .borrow(py)
                .st
                .borrow(py)
                .pbox
                .as_ref()
                .map(|b| b as *const PeriodicBox)?;
            unsafe { Some(&*ptr) }
        })
    }
}

impl RandomBondProvider for SelPy {
    fn num_bonds(&self) -> usize {
        0
    }

    unsafe fn get_bond_unchecked(&self, _i: usize) -> &[usize; 2] {
        unreachable!()
    }
}

impl TimeProvider for SelPy {
    fn get_time(&self) -> f32 {
        Python::attach(|py| self.st.borrow(py).time)
    }
}

impl SaveTopology for SelPy {
    fn iter_atoms_dyn<'a>(&'a self) -> Box<dyn Iterator<Item = &'a Atom> + 'a> {
        Box::new(self.iter_atoms())
    }
}
impl SaveState for SelPy {
    fn iter_pos_dyn<'a>(&'a self) -> Box<dyn Iterator<Item = &'a Pos> + 'a> {
        Box::new(self.iter_pos())
    }
}
impl SaveTopologyState for SelPy {}

impl MeasurePeriodic for SelPy {}

impl RandomMoleculeProvider for SelPy {
    unsafe fn get_molecule_unchecked(&self, i: usize) -> &[usize; 2] {
        let p = Python::attach(|py| {
            self.top.borrow(py).0.get_molecule_unchecked(i) as *const [usize; 2]
        });
        unsafe { &*p }
    }

    fn num_molecules(&self) -> usize {
        Python::attach(|py| self.top.borrow(py).0.num_molecules())
    }
}

impl SelPy {
    fn from_svec(&self, index: SVec) -> Self {
        Python::attach(|py| Self {
            sys: self.sys.clone_ref(py),
            top: self.top.clone_ref(py),
            st: self.st.clone_ref(py),
            index,
        })
    }
}

#[pymethods]
impl SelPy {
    fn __len__(&self) -> usize {
        self.index.len()
    }

    fn __call__(&self, arg: &Bound<'_, PyAny>) -> PyResult<SelPy> {
        let sys = self.sys.borrow(arg.py());
        if let Ok(val) = arg.extract::<String>() {
            let v = val
                .into_sel_index(&*sys, Some(self.index.as_slice()))
                .map_err(|e| anyhow!(e))?;
            Ok(self.from_svec(v))
        } else if let Ok(val) = arg.extract::<(usize, usize)>() {
            let v = (val.0..=val.1)
                .into_sel_index(&*sys, Some(self.index.as_slice()))
                .map_err(|e| anyhow!(e))?;
            Ok(self.from_svec(v))
        } else if let Ok(val) = arg.extract::<Vec<usize>>() {
            let v = val
                .into_sel_index(&*sys, Some(self.index.as_slice()))
                .map_err(|e| anyhow!(e))?;
            Ok(self.from_svec(v))
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
        let n = s.__len__();

        let mut ind = if i < 0 {
            if i.abs() > n as isize {
                return Err(
                    anyhow!("Negative index {i} is out of bounds {}:-1", -(n as isize)).into(),
                );
            }
            n - i.unsigned_abs()
        } else if i >= n as isize {
            return Err(anyhow!("Index {} is out of bounds 0:{}", i, n).into());
        } else {
            i as usize
        };

        // Obtained ind is local to selection, so make it global
        ind += unsafe { s.index.get_unchecked(0) };

        let py = slf.py();

        let all = PySlice::new(py, 0, n as isize, 1);
        let idx = (all, ind);
        let sys = s.sys.bind(py).borrow();
        let st = sys.st.borrow(py);
        let arr = st.coords.bind(py);
        let el = arr.get_item(idx)?.cast_into::<PyArray1<f32>>()?;

        Ok(ParticlePy {
            top: slf.borrow().sys.borrow(py).top.clone_ref(py),
            pos: el.unbind(),
            id: ind,
        }
        .into_py_any(py)?)
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
        numpy::PyArray1::from_iter(py, self.index.iter_index())
    }

    fn get_coord<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f32>> {
        let coord_ptr = self.coords_ptr() as *const f32;
        // We allocate an uninitialized PyArray manually and fill it with data.
        // By doing this we save on unnecessary initiallization and extra allocation
        unsafe {
            let arr = numpy::PyArray2::<f32>::new(py, [3, self.len()], true);
            let arr_ptr = arr.data();
            for i in self.index.iter_index() {
                let pos_ptr = coord_ptr.add(i * 3);
                // This is faster than copying by element with uget_raw()
                std::ptr::copy_nonoverlapping(pos_ptr, arr_ptr.add(i * 3), 3);
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
        let arr_ptr = arr.data();
        let coord_ptr = self.coords_ptr_mut() as *mut f32;

        unsafe {
            for i in self.index.iter_index() {
                let pos_ptr = coord_ptr.add(i * 3);
                std::ptr::copy_nonoverlapping(arr_ptr.add(i * 3), pos_ptr, 3);
            }
        }
        Ok(())
    }

    fn set_state(&mut self, st: &Bound<'_, StatePy>) -> anyhow::Result<Py<StatePy>> {
        let py = st.py();
        let n = self.st.borrow(py).len();
        //let st_ref = cur_st.borrow(py);
        if n == st.borrow().len() {
            let ret = self.st.clone_ref(py);
            self.sys.borrow_mut(py).st = st.as_unbound().clone_ref(py);
            Ok(ret)
        } else {
            bail!("incompatible state")
        }
    }

    // fn set_topology(&mut self, top: &Bound<'_, TopologyPy>) -> anyhow::Result<Py<TopologyPy>> {
    //     let py = top.py();
    //     let cur_top = &self.sys.borrow(py).top;
    //     if cur_top.len() == top.borrow().len() {
    //         self.sys.borrow_mut(py).top = top.as_unbound().clone_ref(py);
    //         Ok(cur_top.clone_ref(py))
    //     } else {
    //         bail!("incompatible topology")
    //     }
    // }

    fn set_state_from(&mut self, arg: &Bound<'_, PyAny>) -> anyhow::Result<Py<StatePy>> {
        let py = arg.py();
        if let Ok(val) = arg.cast::<SystemPy>() {
            let b = val.borrow();
            let s = b.st.bind(py);
            self.set_state(s)
        } else if let Ok(val) = arg.cast::<SelPy>() {
            let b = val.borrow();
            let s = &b.sys.borrow(py).st;
            self.set_state(s.bind(py))
        } else {
            Err(anyhow!(
                "Invalid argument type {} in set_state_from()",
                arg.get_type()
            )
            .into())
        }
    }

    pub fn set_same_chain(&mut self, val: char) {
        AtomIterMutProvider::set_same_chain(self, val)
    }

    pub fn set_same_resname(&mut self, val: &str) {
        AtomIterMutProvider::set_same_resname(self, val)
    }

    pub fn set_same_resid(&mut self, val: i32) {
        AtomIterMutProvider::set_same_resid(self, val)
    }

    pub fn set_same_name(&mut self, val: &str) {
        AtomIterMutProvider::set_same_name(self, val)
    }

    pub fn set_same_mass(&mut self, val: f32) {
        AtomIterMutProvider::set_same_mass(self, val)
    }

    pub fn set_same_bfactor(&mut self, val: f32) {
        AtomIterMutProvider::set_same_bfactor(self, val)
    }

    #[getter]
    fn get_time(&self) -> f32 {
        TimeProvider::get_time(self)
    }

    #[getter]
    fn get_box(&self) -> PeriodicBoxPy {
        Python::attach(|py| {
            let st = self.st.borrow(py);
            let b = st.pbox.as_ref().unwrap();
            PeriodicBoxPy(b.clone())
        })
    }

    #[setter]
    fn set_box(&mut self, b: &PeriodicBoxPy) {
        Python::attach(|py| {
            let mut st = self.st.borrow_mut(py);
            st.pbox = Some(b.0.clone());
        })
    }

    #[setter]
    fn set_time(&mut self, t: f32) {
        Python::attach(|py| {
            let mut st = self.st.borrow_mut(py);
            st.time = t;
        })
    }

    fn set_box_from(&self, sys: Bound<'_, SystemPy>) {
        Python::attach(|py| {
            let sb2 = sys.borrow();
            let st2 = sb2.st.borrow(py);
            let b = st2.pbox.as_ref().unwrap();

            let sb1 = self.sys.borrow(py);
            let mut st1 = sb1.st.borrow_mut(py);
            st1.pbox = Some(b.clone());
        })
    }

    #[pyo3(signature = (dims=[false,false,false]))]
    fn com<'py>(
        &self,
        py: Python<'py>,
        dims: [bool; 3],
    ) -> PyResult<Bound<'py, numpy::PyArray1<f32>>> {
        let pbc_dims = PbcDims::new(dims[0], dims[1], dims[2]);
        Ok(clone_vec_to_pyarray1(
            &MeasurePeriodic::center_of_mass_pbc_dims(self, pbc_dims)
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
            &MeasurePeriodic::center_of_geometry_pbc_dims(self, pbc_dims)
                .map_err(|e| anyhow!(e))?
                .coords,
            py,
        ))
    }

    fn principal_transform(&self) -> anyhow::Result<IsometryTransform> {
        let tr = MeasureMasses::principal_transform(self)?;
        Ok(IsometryTransform(tr))
    }

    fn principal_transform_pbc(&self) -> anyhow::Result<IsometryTransform> {
        let tr = MeasurePeriodic::principal_transform_pbc(self)?;
        Ok(IsometryTransform(tr))
    }

    fn apply_transform(&mut self, tr: &IsometryTransform) {
        ModifyPos::apply_transform(self, &tr.0);
    }

    fn gyration(&self) -> anyhow::Result<f32> {
        Ok(MeasureMasses::gyration(self)?)
    }

    fn gyration_pbc(&self) -> anyhow::Result<f32> {
        Ok(MeasurePeriodic::gyration_pbc(self)?)
    }

    fn min_max<'py>(
        &self,
        py: Python<'py>,
    ) -> (Bound<'py, PyArray1<f32>>, Bound<'py, PyArray1<f32>>) {
        let (min, max) = MeasurePos::min_max(self);
        let minpy = clone_vec_to_pyarray1(&min.coords, py);
        let maxpy = clone_vec_to_pyarray1(&max.coords, py);
        (minpy, maxpy)
    }

    fn inertia<'py>(
        &self,
        py: Python<'py>,
    ) -> anyhow::Result<(
        Bound<'py, numpy::PyArray1<f32>>,
        Bound<'py, numpy::PyArray2<f32>>,
    )> {
        let (moments, axes) = MeasureMasses::inertia(self)?;
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
        let (moments, axes) = MeasurePeriodic::inertia_pbc(self)?;
        let mom = clone_vec_to_pyarray1(&moments, py);
        let ax = axes.to_pyarray(py);
        Ok((mom, ax))
    }

    fn save(&self, fname: &str) -> anyhow::Result<()> {
        Ok(SaveTopologyState::save(self, fname)?)
    }

    fn translate<'py>(&mut self, arg: PyArrayLike1<'py, f32>) -> anyhow::Result<()> {
        let vec: VectorView<f32, Const<3>, Dyn> = arg
            .try_as_matrix()
            .ok_or_else(|| anyhow!("conversion to Vector3 has failed"))?;
        ModifyPos::translate(self, &vec);
        Ok(())
    }

    fn split_resindex(&self, py: Python<'_>) -> Vec<SelPy> {
        AtomPosAnalysis::split_resindex(self)
            .map(|s| SelPy {
                sys: self.sys.clone_ref(py),
                top: self.top.clone_ref(py),
                st: self.st.clone_ref(py),
                index: s.into_svec(),
            })
            .collect()
    }

    fn split_chain(&self, py: Python<'_>) -> Vec<SelPy> {
        self.split(|p| Some(p.atom.chain))
            .map(|s| SelPy {
                sys: self.sys.clone_ref(py),
                top: self.top.clone_ref(py),
                st: self.st.clone_ref(py),
                index: s.into_svec(),
            })
            .collect()
    }

    fn split_molecule(&self, py: Python<'_>) -> Vec<SelPy> {
        self.split_mol_iter()
            .map(|s| SelPy {
                sys: self.sys.clone_ref(py),
                top: self.top.clone_ref(py),
                st: self.st.clone_ref(py),
                index: s.into_svec(),
            })
            .collect()
    }

    fn to_gromacs_ndx(&self, name: &str) -> String {
        self.index.as_gromacs_ndx_str(name)
    }

    // /// operator |
    // fn __or__(&self, rhs: &Sel) -> Sel {
    //     Sel::new_owned(&self.sel | &rhs.sel)
    // }

    // /// operator &
    // fn __and__(&self, rhs: &Sel) -> Sel {
    //     Sel::new_owned(&self.0 & &rhs.0)
    // }

    // /// -= (remove other from self)
    // fn __sub__(&self, rhs: &Sel) -> Sel {
    //     Sel::new_owned(&self.0 - &rhs.0)
    // }

    // /// ~ operator
    // fn __invert__(&self) -> Sel {
    //     Sel::new_owned(!&self.0)
    // }

    // fn sasa(&self) -> SasaResults {
    //     SasaResults(AtomPosAnalysis::sasa(self))
    // }
}

#[pyclass(unsendable, name = "SasaResults")]
struct SasaResultsPy(SasaResults);

#[pymethods]
impl SasaResultsPy {
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
fn fit_transform_py(sel1: &SelPy, sel2: &SelPy) -> anyhow::Result<IsometryTransform> {
    let tr = molar::prelude::fit_transform(sel1, sel2)?;
    Ok(IsometryTransform(tr))
}

#[pyfunction(name = "fit_transform_matching")]
fn fit_transform_matching_py(
    py: Python<'_>,
    sel1: &SelPy,
    sel2: &SelPy,
) -> anyhow::Result<IsometryTransform> {
    let (ind1, ind2) = get_matching_atoms_by_name(sel1, sel2);

    let sub1 = ind1.into_sel_index(sel1, Some(&sel1.index))?;
    let sub2 = ind2.into_sel_index(sel2, Some(&sel2.index))?;

    let sub_sel1 = SelPy {
        sys: sel1.sys.clone_ref(py),
        top: sel1.top.clone_ref(py),
        st: sel1.st.clone_ref(py),
        index: sub1,
    };

    let sub_sel2 = SelPy {
        sys: sel2.sys.clone_ref(py),
        top: sel2.top.clone_ref(py),
        st: sel2.st.clone_ref(py),
        index: sub2,
    };

    let tr = fit_transform(&sub_sel1, &sub_sel2)?;
    Ok(IsometryTransform(tr))
}

#[pyfunction]
fn rmsd_py(sel1: &SelPy, sel2: &SelPy) -> anyhow::Result<f32> {
    Ok(rmsd(sel1, sel2)?)
}

#[pyfunction(name = "rmsd_mw")]
fn rmsd_mw_py(sel1: &SelPy, sel2: &SelPy) -> anyhow::Result<f32> {
    Ok(rmsd_mw(sel1, sel2)?)
}

#[pyclass]
struct ParticleIterator {
    sel: Py<SelPy>,
    cur: isize,
}

#[pymethods]
impl ParticleIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<Py<PyAny>> {
        let ret = Python::attach(|py| {
            let s = slf.sel.bind(py);
            SelPy::__getitem__(s.clone(), slf.cur)
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
    data1: &Bound<'py, SelPy>,
    data2: Option<&Bound<'py, SelPy>>,
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
                res = distance_search_double_pbc(
                    d,
                    sel1.iter_pos(),
                    sel2.iter_pos(),
                    sel1.iter_index(),
                    sel2.iter_index(),
                    &sel1.get_box().0,
                    pbc_dims,
                );
            } else {
                res = distance_search_double(
                    d,
                    &sel1 as &SelPy,
                    &sel2 as &SelPy,
                    sel1.iter_index(),
                    sel2.iter_index(),
                );
            }
        } else {
            if pbc_dims.any() {
                res = distance_search_single_pbc(
                    d,
                    sel1.iter_pos(),
                    sel1.iter_index(),
                    &sel1.get_box().0,
                    pbc_dims,
                );
            } else {
                res = distance_search_single(d, &sel1 as &SelPy, sel1.iter_index());
            }
        }
    } else if let Ok(s) = cutoff.extract::<String>() {
        if s != "vdw" {
            bail!("Unknown cutoff type {s}");
        }

        // VdW cutof
        let vdw1: Vec<f32> = sel1.iter_atoms().map(|a| a.vdw()).collect();

        if sel1.len() != vdw1.len() {
            bail!("Size mismatch 1: {} {}", sel1.len(), vdw1.len());
        }

        if let Some(d2) = data2 {
            let sel2 = d2.borrow();
            let vdw2: Vec<f32> = sel2.iter_atoms().map(|a| a.vdw()).collect();

            if sel2.len() != vdw2.len() {
                bail!("Size mismatch 2: {} {}", sel2.len(), vdw2.len());
            }

            if pbc_dims.any() {
                res = distance_search_double_vdw(&sel1 as &SelPy, &sel2 as &SelPy, &vdw1, &vdw2);
            } else {
                res = distance_search_double_vdw_pbc(
                    sel1.iter_pos(),
                    sel2.iter_pos(),
                    &vdw1,
                    &vdw2,
                    &sel1.get_box().0,
                    pbc_dims,
                );
            }

            // Convert local indices to global
            unsafe {
                for el in &mut res {
                    el.0 = sel1.get_index_unchecked(el.0);
                    el.1 = sel2.get_index_unchecked(el.1);
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

#[pyclass(name = "NdxFile")]
struct NdxFilePy(NdxFile);

#[pymethods]
impl NdxFilePy {
    #[new]
    fn new(fname: &str) -> anyhow::Result<Self> {
        Ok(NdxFilePy(NdxFile::new(fname)?))
    }

    // fn get_group_as_sel(&self, py: Python<'_>, gr_name: &str, sys: &System) -> anyhow::Result<Sel> {
    //     Ok(Sel{
    //         sel: self.0.get_group_as_sel(gr_name, &sys.0)?,
    //         st: Py::clone_ref(&sys.st, py),
    //         top: Py::clone_ref(&sys.top, py),
    //     })
    // }
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
    m.add_class::<AtomPy>()?;
    m.add_class::<ParticlePy>()?;
    m.add_class::<TopologyPy>()?;
    m.add_class::<StatePy>()?;
    m.add_class::<PeriodicBoxPy>()?;
    m.add_class::<FileHandlerPy>()?;
    m.add_class::<SystemPy>()?;
    m.add_class::<SelPy>()?;
    m.add_class::<SasaResultsPy>()?;
    m.add_class::<NdxFilePy>()?;
    m.add_function(wrap_pyfunction!(greeting, m)?)?;
    m.add_function(wrap_pyfunction!(fit_transform_py, m)?)?;
    m.add_function(wrap_pyfunction!(fit_transform_matching_py, m)?)?;
    m.add_function(wrap_pyfunction!(rmsd_py, m)?)?;
    m.add_function(wrap_pyfunction!(rmsd_mw_py, m)?)?;
    m.add_function(wrap_pyfunction!(distance_search, m)?)?;
    //m.add_class::<LipidMolecule>()?;
    //m.add_class::<Membrane>()?;
    //m.add_class::<Histogram1D>()?;
    Ok(())
}
