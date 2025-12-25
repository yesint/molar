use std::path::PathBuf;

use anyhow::{anyhow, bail};
use molar::prelude::*;
//use molar::prelude::*;
use numpy::{
    nalgebra::{self, Const, Dyn, VectorView},
    PyArray1, PyArrayLike1, PyArrayMethods, PyReadonlyArray2, PyUntypedArrayMethods, ToPyArray,
};
use pyo3::{exceptions::PyTypeError, prelude::*, types::PyTuple, IntoPyObjectExt};

mod utils;
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
        Ok((TopologyPy(top.into()), StatePy(st.into())))
    }

    fn read_topology(&mut self) -> anyhow::Result<TopologyPy> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        let top = h.read_topology()?;
        Ok(TopologyPy(top.into()))
    }

    fn read_state(&mut self) -> anyhow::Result<StatePy> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        let st = h.read_state()?;
        Ok(StatePy(st.into()))
    }

    fn write(&mut self, data: Bound<'_, PyAny>) -> anyhow::Result<()> {
        let h = self
            .0
            .as_mut()
            .ok_or_else(|| anyhow!(ALREADY_TRANDFORMED))?;
        if let Ok(s) = data.extract::<PyRef<'_, SystemPy>>() {
            let top = &s.top.borrow(data.py()).0;
            let st = &s.st.borrow(data.py()).0;
            h.write_topology(top)?;
            h.write_state(st)?;
        } else if let Ok(s) = data.extract::<PyRef<'_, SelPy>>() {
            let top = &s.top.borrow(data.py()).0;
            let st = &s.st.borrow(data.py()).0;
            h.write_topology(top)?;
            h.write_state(st)?;
        } else if let Ok(s) = data.cast::<PyTuple>() {
            if s.len() != 2 {
                return Err(anyhow!("tuple must have two elements, not {}", s.len()));
            }
            let top = &s
                .get_item(0)?
                .cast::<TopologyPy>()
                .map_err(|_| anyhow!("first tuple element should be a Topology"))?
                .borrow()
                .0;
            let st = &s
                .get_item(1)?
                .cast::<StatePy>()
                .map_err(|_| anyhow!("second tuple element should be a State"))?
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
        if let Ok(s) = data.extract::<PyRef<'_, SystemPy>>() {
            let top = &s.top.borrow(data.py()).0;
            h.write_topology(top)?;
        } else if let Ok(s) = data.extract::<PyRef<'_, SelPy>>() {
            h.write_topology(&s.top.borrow(data.py()).0)?;
        } else if let Ok(s) = data.extract::<PyRefMut<'_, TopologyPy>>() {
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
        if let Ok(s) = data.extract::<PyRef<'_, SystemPy>>() {
            let st = &s.st.borrow(data.py()).0;
            h.write_state(st)?;
        } else if let Ok(s) = data.extract::<PyRef<'_, SelPy>>() {
            h.write_state(&s.st.borrow(data.py()).0)?;
        } else if let Ok(s) = data.extract::<PyRefMut<'_, StatePy>>() {
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

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<Py<PyAny>> {
        let st = slf.1.as_mut().unwrap().next().map(|st| StatePy(st.into()));
        if st.is_some() {
            Python::attach(|py| Some(st.unwrap().into_py_any(py)))
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

#[pyclass(unsendable, sequence, name = "System")]
struct SystemPy {
    top: Py<TopologyPy>,
    st: Py<StatePy>,
}

impl LenProvider for SystemPy {
    fn len(&self) -> usize {
        Python::attach(|py| self.top.borrow(py).0.len())
    }
}

impl IndexProvider for SystemPy {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        i
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> {
        0..self.len()
    }
}

impl AtomPosAnalysis for SystemPy {
    fn atoms_ptr(&self) -> *const Atom {
        Python::attach(|py| self.top.borrow(py).0.atoms.as_ptr())
    }

    fn coords_ptr(&self) -> *const Pos {
        Python::attach(|py| self.st.borrow(py).0.coords.as_ptr())
    }
}

impl AtomPosAnalysisMut for SystemPy {
    fn atoms_ptr_mut(&mut self) -> *mut Atom {
        Python::attach(|py| self.top.borrow_mut(py).0.atoms.as_mut_ptr())
    }

    fn coords_ptr_mut(&mut self) -> *mut Pos {
        Python::attach(|py| self.st.borrow_mut(py).0.coords.as_mut_ptr())
    }
}

impl NonAtomPosAnalysis for SystemPy {
    fn top_ptr(&self) -> *const Topology {
        Python::attach(|py| &self.top.borrow(py).0 as *const Topology)
    }

    fn st_ptr(&self) -> *const State {
        Python::attach(|py| &self.st.borrow(py).0 as *const State)
    }
}

impl NonAtomPosAnalysisMut for SystemPy {
    fn st_ptr_mut(&mut self) -> *mut State {
        Python::attach(|py: Python<'_>| &mut self.st.borrow_mut(py).0 as *mut State)
    }

    fn top_ptr_mut(&mut self) -> *mut Topology {
        Python::attach(|py| &mut self.top.borrow_mut(py).0 as *mut Topology)
    }
}

impl SaveTopology for SystemPy {}
impl SaveState for SystemPy {}
impl SaveTopologyState for SystemPy {}

#[pymethods]
impl SystemPy {
    #[new]
    #[pyo3(signature = (*py_args))]
    fn new(py_args: &Bound<'_, PyTuple>) -> PyResult<Self> {
        let py = py_args.py();
        if py_args.len() == 1 {
            // From file
            let mut fh = molar::io::FileHandler::open(&py_args.get_item(0)?.extract::<String>()?)
                .map_err(|e| anyhow!(e))?;
            let (top, st) = fh.read().map_err(|e| anyhow!(e))?;
            Ok(SystemPy {
                top: Py::new(py, TopologyPy(top))?,
                st: Py::new(py, StatePy(st))?,
            })
        } else if py_args.len() == 2 {
            // From existing Topology and State python objects
            let top = py_args.get_item(0)?;
            let st = py_args.get_item(1)?;
            let n1 = top.cast::<TopologyPy>()?.borrow().0.len();
            let n2 = st.cast::<StatePy>()?.borrow().0.len();
            if n1 != n2 {
                Err(PyTypeError::new_err(
                    "topology and state are of different size",
                ))
            } else {
                Ok(SystemPy {
                    top: Py::clone_ref(&top.cast::<TopologyPy>()?.as_unbound(), py),
                    st: Py::clone_ref(&st.cast::<StatePy>()?.as_unbound(), py),
                })
            }
        } else {
            // Empty System
            Ok(SystemPy {
                top: Py::new(py_args.py(), TopologyPy(Default::default()))?,
                st: Py::new(py_args.py(), StatePy(Default::default()))?,
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
        let py = slf.py();
        let bs = slf.borrow();
        let top = bs.top.borrow(py);
        let st = bs.st.borrow(py);

        let index = if let Some(arg) = arg {
            // Argument present
            if let Ok(val) = arg.extract::<String>() {
                if val.is_empty() {
                    // Select all on empty string
                    (0..top.0.len()).into_sel_index(&top.0, &st.0, None)?
                } else {
                    // Otherwise do normal textual selection
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
            // No argument, select all
            (0..top.0.len()).into_sel_index(&top.0, &st.0, None)?
        };

        Ok(SelPy {
            top: Py::clone_ref(&bs.top, py),
            st: Py::clone_ref(&bs.st, py),
            index,
        })
    }

    fn set_state(&mut self, st: &Bound<'_, StatePy>) -> anyhow::Result<Py<StatePy>> {
        if self.st.borrow(st.py()).0.interchangeable(&st.borrow().0) {
            let ret = Py::clone_ref(&self.st, st.py());
            self.st = Py::clone_ref(st.as_unbound(), st.py());
            Ok(ret)
        } else {
            bail!("incompatible state")
        }
    }

    fn set_topology(&mut self, top: &Bound<'_, TopologyPy>) -> anyhow::Result<Py<TopologyPy>> {
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
        Ok(SaveTopologyState::save(self, fname)?)
    }

    fn remove<'py>(slf: &Bound<'py, Self>, arg: &Bound<'py, PyAny>) -> anyhow::Result<()> {
        let py = slf.py();
        let slf_b = slf.borrow();
        if let Ok(sel) = arg.cast::<SelPy>() {
            // Selection provided
            let sb = sel.borrow();
            slf_b.top.borrow_mut(py).0.remove_atoms(sb.iter_index())?;
            slf_b.st.borrow_mut(py).0.remove_coords(sb.iter_index())?;
            Ok(())
        } else {
            let sel = Self::__call__(slf, Some(arg))?;
            slf_b.top.borrow_mut(py).0.remove_atoms(sel.iter_index())?;
            slf_b.st.borrow_mut(py).0.remove_coords(sel.iter_index())?;
            Ok(())
        }
    }

    #[pyo3(signature = (*args))]
    fn append<'py>(slf: &Bound<'py, Self>, args: &Bound<'py, PyTuple>) -> anyhow::Result<()> {
        let py = slf.py();
        let slf_b = slf.borrow();
        let topb = &mut slf_b.top.borrow_mut(py).0;
        let stb = &mut slf_b.st.borrow_mut(py).0;

        if args.len() == 1 {
            let arg = args.get_item(0)?;
            // In the future other types can be used as well
            let sel = if let Ok(sel) = arg.cast::<SelPy>() {
                sel
            } else {
                &Bound::new(py, Self::__call__(slf, Some(&arg))?)?
            };

            let sb = sel.borrow();
            topb.add_atoms(sb.iter_atoms().cloned());
            stb.add_coords(sb.iter_pos().cloned());

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
            topb.add_atoms(std::iter::once(&ab.0).cloned());
            stb.add_coords(std::iter::once(Pos::new(v.x, v.y, v.z)));
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
        BoxProvider::get_box(&self.st.borrow(py).0)
            .map(|b| PeriodicBoxPy(b.clone()))
            .unwrap()
    }

    #[setter]
    fn set_box(&mut self, py: Python<'_>, b: &PeriodicBoxPy) {
        *self.st.borrow_mut(py).0.get_box_mut().unwrap() = b.0.clone();
    }

    #[setter]
    fn set_time(&mut self, t: f32) {
        TimeMutProvider::set_time(self, t);
    }

    fn set_box_from(&self, sys: Bound<'_, SystemPy>) {
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

#[pyclass(sequence, name = "Sel")]
struct SelPy {
    top: Py<TopologyPy>,
    st: Py<StatePy>,
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
        Python::attach(|py| self.st.borrow(py).0.coords.as_ptr())
    }
}

impl AtomPosAnalysisMut for SelPy {
    fn atoms_ptr_mut(&mut self) -> *mut Atom {
        Python::attach(|py| self.top.borrow_mut(py).0.atoms.as_mut_ptr())
    }

    fn coords_ptr_mut(&mut self) -> *mut Pos {
        Python::attach(|py| self.st.borrow_mut(py).0.coords.as_mut_ptr())
    }
}

impl NonAtomPosAnalysis for SelPy {
    fn top_ptr(&self) -> *const Topology {
        Python::attach(|py| &self.top.borrow(py).0 as *const Topology)
    }

    fn st_ptr(&self) -> *const State {
        Python::attach(|py| &self.st.borrow(py).0 as *const State)
    }
}

impl NonAtomPosAnalysisMut for SelPy {
    fn st_ptr_mut(&mut self) -> *mut State {
        Python::attach(|py| &mut self.st.borrow_mut(py).0 as *mut State)
    }

    fn top_ptr_mut(&mut self) -> *mut Topology {
        Python::attach(|py| &mut self.top.borrow_mut(py).0 as *mut Topology)
    }
}

impl SaveTopology for SelPy {}
impl SaveState for SelPy {}
impl SaveTopologyState for SelPy {}

impl SelPy {
    fn from_svec(&self, py: Python<'_>, index: SVec) -> Self {
        Self {
            top: Py::clone_ref(&self.top, py),
            st: Py::clone_ref(&self.st, py),
            index,
        }
    }
}

#[pymethods]
impl SelPy {
    fn __len__(&self) -> usize {
        self.index.len()
    }

    fn __call__(&self, arg: &Bound<'_, PyAny>) -> PyResult<SelPy> {
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
        let n = s.__len__();
        let ind = if i < 0 {
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

        let py = slf.py();

        Ok(ParticlePy {
            top: Py::clone_ref(&s.top, py),
            st: Py::clone_ref(&s.st, py),
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
        let coord_ptr = self.st.borrow(py).0.coords.as_ptr() as *const f32;
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
            for i in self.index.iter_index() {
                let pos_ptr = coord_ptr.add(i * 3);
                std::ptr::copy_nonoverlapping(arr_ptr.add(i * 3), pos_ptr, 3);
            }
        }

        Ok(())
    }

    fn set_state(&mut self, st: &Bound<'_, StatePy>) -> anyhow::Result<Py<StatePy>> {
        if self.st.borrow(st.py()).0.interchangeable(&st.borrow().0) {
            let ret = Py::clone_ref(&self.st, st.py());
            self.st = Py::clone_ref(st.as_unbound(), st.py());
            Ok(ret)
        } else {
            bail!("incompatible state")
        }
    }

    fn set_topology(&mut self, top: &Bound<'_, TopologyPy>) -> anyhow::Result<Py<TopologyPy>> {
        if self.top.borrow(top.py()).0.interchangeable(&top.borrow().0) {
            let ret = Py::clone_ref(&self.top, top.py());
            self.top = Py::clone_ref(top.as_unbound(), top.py());
            Ok(ret)
        } else {
            bail!("incompatible topology")
        }
    }

    fn set_state_from(&mut self, arg: &Bound<'_, PyAny>) -> anyhow::Result<Py<StatePy>> {
        if let Ok(val) = arg.cast::<SystemPy>() {
            let b = val.borrow();
            let s = b.st.bind(arg.py());
            self.set_state(s)
        } else if let Ok(val) = arg.cast::<SelPy>() {
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
    fn get_box(&self, py: Python<'_>) -> PeriodicBoxPy {
        BoxProvider::get_box(&self.st.borrow(py).0)
            .map(|b| PeriodicBoxPy(b.clone()))
            .unwrap()
    }

    #[setter]
    fn set_box(&mut self, py: Python<'_>, b: &PeriodicBoxPy) {
        *self.st.borrow_mut(py).0.get_box_mut().unwrap() = b.0.clone();
    }

    #[setter]
    fn set_time(&mut self, t: f32) {
        TimeMutProvider::set_time(self, t);
    }

    fn set_box_from(&self, sys: Bound<'_, SystemPy>) {
        *self.st.borrow_mut(sys.py()).0.get_box_mut().unwrap() = sys
            .borrow()
            .st
            .borrow(sys.py())
            .0
            .get_box()
            .unwrap()
            .clone();
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
                top: Py::clone_ref(&self.top, py),
                st: Py::clone_ref(&self.st, py),
                index: s.into_svec(),
            })
            .collect()
    }

    fn split_chain(&self, py: Python<'_>) -> Vec<SelPy> {
        self.split(|p| Some(p.atom.chain))
            .map(|s| SelPy {
                top: Py::clone_ref(&self.top, py),
                st: Py::clone_ref(&self.st, py),
                index: s.into_svec(),
            })
            .collect()
    }

    fn split_molecule(&self, py: Python<'_>) -> Vec<SelPy> {
        self.split_mol_iter()
            .map(|s| SelPy {
                top: Py::clone_ref(&self.top, py),
                st: Py::clone_ref(&self.st, py),
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
fn fit_transform_matching_py(sel1: &SelPy, sel2: &SelPy) -> anyhow::Result<IsometryTransform> {
    let tr = fit_transform_matching(sel1, sel2)?;
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
                    &sel1.get_box(py).0,
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
                    &sel1.get_box(py).0,
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
                    &sel1.get_box(py).0,
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
