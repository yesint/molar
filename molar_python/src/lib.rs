use std::{path::PathBuf, rc::Rc};

use anyhow::{anyhow, bail};
use molar::{
    core::{
        AtomIterProvider, AtomPosAnalysis, AtomPosAnalysisMut, BoxMutProvider, BoxProvider, IndexProvider, LenProvider, MeasureMasses, MeasurePeriodic, MeasurePos, ModifyPos, NonAtomPosAnalysis, NonAtomPosAnalysisMut, PosIterProvider, SelectionDef, TimeMutProvider, TimeProvider, rmsd_mw
    },
    io::{StateWrite, TopologyStateWrite, TopologyWrite},
};
//use molar::prelude::*;
use numpy::{
    nalgebra::{self, Const, Dyn, VectorView},
    PyArrayLike1, PyArrayMethods, PyReadonlyArray2, PyUntypedArrayMethods, ToPyArray,
};
use pyo3::{exceptions::PyTypeError, prelude::*, types::PyTuple, IntoPyObjectExt};

mod utils;
use utils::*;

mod atom;
use atom::Atom;

mod particle;
use particle::Particle;

mod periodic_box;
use periodic_box::PeriodicBox;

// mod membrane;
// use membrane::*;

mod topology_state;
use topology_state::*;
//-------------------------------------------

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
            let top = &s.top.borrow(data.py()).0;
            let st = &s.st.borrow(data.py()).0;
            h.write_topology(top)?;
            h.write_state(st)?;
        } else if let Ok(s) = data.extract::<PyRef<'_, Sel>>() {
            let top = &s.top.borrow(data.py()).0;
            let st = &s.st.borrow(data.py()).0;
            h.write_topology(top)?;
            h.write_state(st)?;
        } else if let Ok(s) = data.downcast::<PyTuple>() {
            if s.len() != 2 {
                return Err(anyhow!("Tuple must have two elements"));
            }
            let top = &s
                .iter()
                .next()
                .unwrap()
                .downcast::<Topology>()
                .unwrap()
                .borrow()
                .0;
            let st = &s
                .iter()
                .next()
                .unwrap()
                .downcast::<State>()
                .unwrap()
                .borrow()
                .0;
            h.write_topology(top)?;
            h.write_state(st)?;
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
            let top = &s.top.borrow(data.py()).0;
            h.write_topology(top)?;
        } else if let Ok(s) = data.extract::<PyRef<'_, Sel>>() {
            h.write_topology(&s.top.borrow(data.py()).0)?;
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
            let st = &s.st.borrow(data.py()).0;
            h.write_state(st)?;
        } else if let Ok(s) = data.extract::<PyRef<'_, Sel>>() {
            h.write_state(&s.st.borrow(data.py()).0)?;
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
struct System {
    top: Py<Topology>,
    st: Py<State>,
}

impl LenProvider for System {
    fn len(&self) -> usize {
        Python::with_gil(|py| self.top.borrow(py).0.len())
    }
}

impl IndexProvider for System {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        i
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        0..self.len()
    }
}

impl AtomPosAnalysis for System {
    fn atom_ptr(&self) -> *const molar::core::Atom {
        Python::with_gil(|py| self.top.borrow(py).0.atoms.as_ptr())
    }

    fn pos_ptr(&self) -> *const molar::core::Pos {
        Python::with_gil(|py| self.st.borrow(py).0.coords.as_ptr())
    }
}

impl AtomPosAnalysisMut for System {
    fn atom_mut_ptr(&mut self) -> *mut molar::core::Atom {
        Python::with_gil(|py| self.top.borrow_mut(py).0.atoms.as_mut_ptr())
    }

    fn pos_mut_ptr(&mut self) -> *mut molar::core::Pos {
        Python::with_gil(|py| self.st.borrow_mut(py).0.coords.as_mut_ptr())
    }
}

impl NonAtomPosAnalysis for System {
    fn top_ref(&self) -> &molar::core::Topology {
        Python::with_gil(|py| unsafe { &*(&self.top.borrow(py).0 as *const molar::core::Topology) })
    }

    fn st_ref(&self) -> &molar::core::State {
        Python::with_gil(|py| unsafe { &*(&self.st.borrow(py).0 as *const molar::core::State) })
    }
}

impl NonAtomPosAnalysisMut for System {
    fn st_ref_mut(&mut self) -> &mut molar::core::State {
        Python::with_gil(|py| unsafe {
            &mut *(&mut self.st.borrow_mut(py).0 as *mut molar::core::State)
        })
    }

    fn top_ref_mut(&mut self) -> &mut molar::core::Topology {
        Python::with_gil(|py| unsafe {
            &mut *(&mut self.top.borrow_mut(py).0 as *mut molar::core::Topology)
        })
    }
}

impl TopologyWrite for System {}
impl StateWrite for System {}
impl TopologyStateWrite for System {}

#[pymethods]
impl System {
    #[new]
    #[pyo3(signature = (*py_args))]
    fn new(py_args: &Bound<'_, PyTuple>) -> PyResult<Self> {
        let py = py_args.py();
        if py_args.len() == 1 {
            // From file
            let mut fh = molar::io::FileHandler::open(&py_args.get_item(0)?.extract::<String>()?)
                .map_err(|e| anyhow!(e))?;
            let (top, st) = fh.read().map_err(|e| anyhow!(e))?;
            Ok(System {
                top: Py::new(py, Topology(top))?,
                st: Py::new(py, State(st))?,
            })
        } else if py_args.len() == 2 {
            // From existing Topology and State python objects
            let top = py_args.get_item(0)?;
            let st = py_args.get_item(1)?;
            let n1 = top.downcast::<Topology>()?.borrow().0.len();
            let n2 = st.downcast::<State>()?.borrow().0.len();
            if n1 != n2 {
                Err(PyTypeError::new_err(
                    "topology and state are of different size",
                ))
            } else {
                Ok(System {
                    top: Py::clone_ref(&top.downcast::<Topology>()?.as_unbound(), py),
                    st: Py::clone_ref(&st.downcast::<State>()?.as_unbound(), py),
                })
            }
        } else {
            // Empty System
            Ok(System {
                top: Py::new(py_args.py(), Topology(Default::default()))?,
                st: Py::new(py_args.py(), State(Default::default()))?,
            })
        }
    }

    fn __len__(&self, py: Python<'_>) -> usize {
        self.top.borrow(py).0.len()
    }

    // fn select_all(slf: &Bound<Self>) -> anyhow::Result<Sel> {
    //     let bs = slf.borrow();
    //     let top = bs.top.borrow(slf.py());
    //     let st = bs.st.borrow(slf.py());
    //     let v = (0..top.0.len()).into_sel_index(&top.0, &st.0, None)?;
    //     Ok(Sel {
    //         top: Py::clone_ref(&bs.top, slf.py()),
    //         st: Py::clone_ref(&bs.st, slf.py()),
    //         sel: molar::core::Sel::from_svec(v).unwrap(),
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
    //         sel: molar::core::Sel::from_svec(v).unwrap(),
    //     })
    // }

    #[pyo3(signature = (arg=None))]
    fn __call__(slf: &Bound<Self>, arg: Option<&Bound<'_, PyAny>>) -> anyhow::Result<Sel> {
        let bs = slf.borrow();
        let top = bs.top.borrow(slf.py());
        let st = bs.st.borrow(slf.py());

        let v = if let Some(arg) = arg {
            if let Ok(val) = arg.extract::<String>() {
                if val.is_empty() {
                    (0..top.0.len()).into_sel_index(&top.0, &st.0, None)?
                } else {
                    val.into_sel_index(&top.0, &st.0, None)?
                }
            } else if let Ok(val) = arg.extract::<(usize, usize)>() {
                (val.0..val.1).into_sel_index(&top.0, &st.0, None)?
            } else if let Ok(val) = arg.extract::<Vec<usize>>() {
                val.into_sel_index(&top.0, &st.0, None)?
            } else {
                bail!(
                    "Invalid argument type {} when creating selection",
                    arg.get_type()
                )
            }
        } else {
            (0..top.0.len()).into_sel_index(&top.0, &st.0, None)?
        };

        Ok(Sel {
            top: Py::clone_ref(&bs.top, slf.py()),
            st: Py::clone_ref(&bs.st, slf.py()),
            sel: molar::core::Sel::from_svec(v).unwrap(),
        })
    }

    fn set_state(&mut self, st: &Bound<'_, State>) -> anyhow::Result<Py<State>> {
        if self.st.borrow(st.py()).0.interchangeable(&st.borrow().0) {
            let ret = Py::clone_ref(&self.st, st.py());
            self.st = Py::clone_ref(st.as_unbound(), st.py());
            Ok(ret)
        } else {
            bail!("incompatible state")
        }
    }

    fn set_topology(&mut self, top: &Bound<'_, Topology>) -> anyhow::Result<Py<Topology>> {
        if self.top.borrow(top.py()).0.interchangeable(&top.borrow().0) {
            let ret = Py::clone_ref(&self.top, top.py());
            self.top = Py::clone_ref(top.as_unbound(), top.py());
            Ok(ret)
        } else {
            bail!("incompatible topology")
        }
    }

    // fn get_state(&self) -> State {
    //     State(self.0.get_state())
    // }

    // fn get_topology(&self) -> Topology {
    //     Topology(self.0.get_topology())
    // }

    fn save(&self, fname: &str) -> anyhow::Result<()> {
        Ok(TopologyStateWrite::save(self, fname)?)
    }

    fn remove(slf: &Bound<'_, Self>, arg: &Bound<'_, PyAny>) -> anyhow::Result<()> {
        // In the future other types can be used as well
        if let Ok(sel) = arg.downcast::<Sel>() {
            let sb = sel.borrow();
            slf.borrow_mut()
                .top
                .borrow_mut(arg.py())
                .0
                .remove_atoms(sb.iter_index())?;
            slf.borrow_mut()
                .st
                .borrow_mut(arg.py())
                .0
                .remove_coords(sb.iter_index())?;
            Ok(())
        } else {
            let sel = Self::__call__(slf, Some(arg))?;
            slf.borrow_mut()
                .top
                .borrow_mut(arg.py())
                .0
                .remove_atoms(sel.iter_index())?;
            slf.borrow_mut()
                .st
                .borrow_mut(arg.py())
                .0
                .remove_coords(sel.iter_index())?;
            Ok(())
        }
    }

    fn append(slf: &Bound<'_, Self>, arg: &Bound<'_, PyAny>) -> anyhow::Result<()> {
        // In the future other types can be used as well
        if let Ok(sel) = arg.downcast::<Sel>() {
            let sb = sel.borrow();
            slf.borrow_mut()
                .top
                .borrow_mut(arg.py())
                .0
                .add_atoms(sb.iter_atoms().cloned());
            slf.borrow_mut()
                .st
                .borrow_mut(arg.py())
                .0
                .add_coords(sb.iter_pos().cloned());
            Ok(())
        } else {
            let sel = Self::__call__(slf, Some(arg))?;
            slf.borrow_mut()
                .top
                .borrow_mut(arg.py())
                .0
                .add_atoms(sel.iter_atoms().cloned());
            slf.borrow_mut()
                .st
                .borrow_mut(arg.py())
                .0
                .add_coords(sel.iter_pos().cloned());
            Ok(())
        }
    }

    #[getter]
    fn get_time(&self) -> f32 {
        TimeProvider::get_time(self)
    }

    #[setter]
    fn set_time(&mut self, t: f32) {
        TimeMutProvider::set_time(self,t);
    }

    fn set_box_from(&self, sys: Bound<'_, System>) {
        *self.st.borrow_mut(sys.py()).0.get_box_mut().unwrap() = sys
            .borrow()
            .st
            .borrow(sys.py())
            .0
            .get_box()
            .unwrap()
            .clone();
    }
}

//====================================

#[pyclass(sequence)]
struct Sel {
    top: Py<Topology>,
    st: Py<State>,
    sel: molar::core::Sel,
}

impl LenProvider for Sel {
    fn len(&self) -> usize {
        self.sel.len()
    }
}

impl IndexProvider for Sel {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        self.sel.get_index_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.sel.iter_index()
    }
}

impl AtomPosAnalysis for Sel {
    fn atom_ptr(&self) -> *const molar::core::Atom {
        Python::with_gil(|py| self.top.borrow(py).0.atoms.as_ptr())
    }

    fn pos_ptr(&self) -> *const molar::core::Pos {
        Python::with_gil(|py| self.st.borrow(py).0.coords.as_ptr())
    }
}

impl AtomPosAnalysisMut for Sel {
    fn atom_mut_ptr(&mut self) -> *mut molar::core::Atom {
        Python::with_gil(|py| self.top.borrow_mut(py).0.atoms.as_mut_ptr())
    }

    fn pos_mut_ptr(&mut self) -> *mut molar::core::Pos {
        Python::with_gil(|py| self.st.borrow_mut(py).0.coords.as_mut_ptr())
    }
}

impl NonAtomPosAnalysis for Sel {
    fn top_ref(&self) -> &molar::core::Topology {
        Python::with_gil(|py| unsafe { &*(&self.top.borrow(py).0 as *const molar::core::Topology) })
    }

    fn st_ref(&self) -> &molar::core::State {
        Python::with_gil(|py| unsafe { &*(&self.st.borrow(py).0 as *const molar::core::State) })
    }
}

impl NonAtomPosAnalysisMut for Sel {
    fn st_ref_mut(&mut self) -> &mut molar::core::State {
        Python::with_gil(|py| unsafe {
            &mut *(&mut self.st.borrow_mut(py).0 as *mut molar::core::State)
        })
    }

    fn top_ref_mut(&mut self) -> &mut molar::core::Topology {
        Python::with_gil(|py| unsafe {
            &mut *(&mut self.top.borrow_mut(py).0 as *mut molar::core::Topology)
        })
    }
}

impl TopologyWrite for Sel {}
impl StateWrite for Sel {}
impl TopologyStateWrite for Sel {}

impl Sel {
    fn from_svec(&self, py: Python<'_>, v: molar::core::SVec) -> Self {
        Self {
            top: Py::clone_ref(&self.top, py),
            st: Py::clone_ref(&self.st, py),
            sel: molar::core::Sel::from_svec(v).unwrap(),
        }
    }
}

#[pymethods]
impl Sel {
    fn __len__(&self) -> usize {
        self.sel.len()
    }

    fn __call__(&self, arg: &Bound<'_, PyAny>) -> PyResult<Sel> {
        let top = self.top.borrow(arg.py());
        let st = self.st.borrow(arg.py());
        if let Ok(val) = arg.extract::<String>() {
            let v = val
                .into_sel_index(&top.0, &st.0, None)
                .map_err(|e| anyhow!(e))?;
            Ok(self.from_svec(arg.py(), v))
        } else if let Ok(val) = arg.extract::<(usize, usize)>() {
            let v = (val.0..=val.1)
                .into_sel_index(&top.0, &st.0, None)
                .map_err(|e| anyhow!(e))?;
            Ok(self.from_svec(arg.py(), v))
        } else if let Ok(val) = arg.extract::<Vec<usize>>() {
            let v = val
                .into_sel_index(&top.0, &st.0, None)
                .map_err(|e| anyhow!(e))?;
            Ok(self.from_svec(arg.py(), v))
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

        Ok(Particle {
            top: Py::clone_ref(&s.top, slf.py()),
            st: Py::clone_ref(&s.st, slf.py()),
            id: ind,
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
        numpy::PyArray1::from_iter(py, self.sel.iter_index())
    }

    fn get_coord<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f32>> {
        let coord_ptr = self.st.borrow(py).0.coords.as_ptr() as *const f32;
        // We allocate an uninitialized PyArray manually and fill it with data.
        // By doing this we save on unnecessary initiallization and extra allocation
        unsafe {
            let arr = numpy::PyArray2::<f32>::new(py, [3, self.len()], true);
            let arr_ptr = arr.data();
            for i in self.sel.iter_index() {
                let pos_ptr = coord_ptr.add(i * 3);
                // This is faster than copying by element with uget_raw()
                std::ptr::copy_nonoverlapping(pos_ptr, arr_ptr.add(i * 3), 3);
            }
            arr
        }
    }

    fn set_coord(&mut self, py: Python<'_>, arr: PyReadonlyArray2<f32>) -> PyResult<()> {
        // Check if the shape is correct
        if arr.shape() != [3, self.__len__()] {
            return Err(anyhow!(
                "Array shape must be [3, {}], not {:?}",
                self.__len__(),
                arr.shape()
            ))?;
        }
        let arr_ptr = arr.data();
        let coord_ptr = self.st.borrow(py).0.coords.as_ptr() as *mut f32;

        unsafe {
            for i in self.sel.iter_index() {
                let pos_ptr = coord_ptr.add(i * 3);
                std::ptr::copy_nonoverlapping(arr_ptr.add(i * 3), pos_ptr, 3);
            }
        }

        Ok(())
    }

    fn set_state(&mut self, st: &Bound<'_, State>) -> anyhow::Result<Py<State>> {
        if self.st.borrow(st.py()).0.interchangeable(&st.borrow().0) {
            let ret = Py::clone_ref(&self.st, st.py());
            self.st = Py::clone_ref(st.as_unbound(), st.py());
            Ok(ret)
        } else {
            bail!("incompatible state")
        }
    }

    fn set_topology(&mut self, top: &Bound<'_, Topology>) -> anyhow::Result<Py<Topology>> {
        if self.top.borrow(top.py()).0.interchangeable(&top.borrow().0) {
            let ret = Py::clone_ref(&self.top, top.py());
            self.top = Py::clone_ref(top.as_unbound(), top.py());
            Ok(ret)
        } else {
            bail!("incompatible topology")
        }
    }

    fn set_state_from(&mut self, arg: &Bound<'_, PyAny>) -> anyhow::Result<Py<State>> {
        if let Ok(val) = arg.downcast::<System>() {
            let b = val.borrow();
            let s = b.st.bind(arg.py());
            self.set_state(s)
        } else if let Ok(val) = arg.downcast::<Sel>() {
            let b = val.borrow();
            let s = b.st.bind(arg.py());
            self.set_state(s)
        } else {
            Err(anyhow!(
                "Invalid argument type {} in set_state_from()",
                arg.get_type()
            )
            .into())
        }
    }

    pub fn set_same_chain(&mut self, val: char) {
        molar::core::AtomIterMutProvider::set_same_chain(self, val)
    }

    pub fn set_same_resname(&mut self, val: &str) {
        molar::core::AtomIterMutProvider::set_same_resname(self, val)
    }

    pub fn set_same_resid(&mut self, val: i32) {
        molar::core::AtomIterMutProvider::set_same_resid(self, val)
    }

    pub fn set_same_name(&mut self, val: &str) {
        molar::core::AtomIterMutProvider::set_same_name(self, val)
    }

    pub fn set_same_mass(&mut self, val: f32) {
        molar::core::AtomIterMutProvider::set_same_mass(self, val)
    }

    pub fn set_same_bfactor(&mut self, val: f32) {
        molar::core::AtomIterMutProvider::set_same_bfactor(self, val)
    }

    #[getter]
    fn get_time(&self) -> f32 {
        TimeProvider::get_time(self)
    }

    #[pyo3(signature = (dims=[false,false,false]))]
    fn com<'py>(
        &self,
        py: Python<'py>,
        dims: [bool; 3],
    ) -> PyResult<Bound<'py, numpy::PyArray1<f32>>> {
        let pbc_dims = molar::core::PbcDims::new(dims[0], dims[1], dims[2]);
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
        let pbc_dims = molar::core::PbcDims::new(dims[0], dims[1], dims[2]);
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
        Ok(TopologyStateWrite::save(self, fname)?)
    }

    fn translate<'py>(&mut self, arg: PyArrayLike1<'py, f32>) -> anyhow::Result<()> {
        let vec: VectorView<f32, Const<3>, Dyn> = arg
            .try_as_matrix()
            .ok_or_else(|| anyhow!("conversion to Vector3 has failed"))?;
        ModifyPos::translate(self, &vec);
        Ok(())
    }

    fn split_resindex(&self, py: Python<'_>) -> Vec<Sel> {
        self.split_resindex_iter()
            .map(|s| Sel {
                top: Py::clone_ref(&self.top, py),
                st: Py::clone_ref(&self.st, py),
                sel: s,
            })
            .collect()
    }

    fn split_chain(&self, py: Python<'_>) -> Vec<Sel> {
        self.split_iter(|p| Some(p.atom.chain))
            .map(|s| Sel {
                top: Py::clone_ref(&self.top, py),
                st: Py::clone_ref(&self.st, py),
                sel: s,
            })
            .collect()
    }

    fn split_molecule(&self, py: Python<'_>) -> Vec<Sel> {
        self.split_mol_iter()
            .map(|s| Sel {
                top: Py::clone_ref(&self.top, py),
                st: Py::clone_ref(&self.st, py),
                sel: s,
            })
            .collect()
    }

    fn to_gromacs_ndx(&self, name: &str) -> String {
        self.sel.as_gromacs_ndx_str(name)
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
    let tr = molar::prelude::fit_transform(sel1, sel2)?;
    Ok(IsometryTransform(tr))
}

#[pyfunction(name = "fit_transform_matching")]
fn fit_transform_matching_py(sel1: &Sel, sel2: &Sel) -> anyhow::Result<IsometryTransform> {
    let tr = molar::core::fit_transform_matching(sel1, sel2)?;
    Ok(IsometryTransform(tr))
}

#[pyfunction]
fn rmsd(sel1: &Sel, sel2: &Sel) -> anyhow::Result<f32> {
    Ok(MeasurePos::rmsd(sel1, sel2)?)
}

#[pyfunction(name = "rmsd_mw")]
fn rmsd_mw_py(sel1: &Sel, sel2: &Sel) -> anyhow::Result<f32> {
    Ok(rmsd_mw(sel1, sel2)?)
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
    let pbc_dims = molar::core::PbcDims::new(dims[0], dims[1], dims[2]);
    let sel1 = data1.borrow();

    if let Ok(d) = cutoff.extract::<f32>() {
        // Distance cutoff
        if let Some(d2) = data2 {
            let sel2 = d2.borrow();
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
                res = molar::core::distance_search_single(d, sel1.iter_pos(), sel1.iter_index());
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

#[pyclass]
struct NdxFile(molar::core::NdxFile);

#[pymethods]
impl NdxFile {
    #[new]
    fn new(fname: &str) -> anyhow::Result<Self> {
        Ok(NdxFile(molar::core::NdxFile::new(fname)?))
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
    m.add_function(wrap_pyfunction!(greeting, m)?)?;
    m.add_function(wrap_pyfunction!(fit_transform_py, m)?)?;
    m.add_function(wrap_pyfunction!(fit_transform_matching_py, m)?)?;
    m.add_function(wrap_pyfunction!(rmsd, m)?)?;
    m.add_function(wrap_pyfunction!(rmsd_mw_py, m)?)?;
    m.add_function(wrap_pyfunction!(distance_search, m)?)?;
    //m.add_class::<LipidMolecule>()?;
    //m.add_class::<Membrane>()?;
    //m.add_class::<Histogram1D>()?;
    Ok(())
}
