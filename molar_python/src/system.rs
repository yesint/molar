use std::cell::UnsafeCell;

use molar::prelude::*;
use numpy::nalgebra::{Const, VectorView};
use numpy::{PyArray1, PyArrayLike1};
use pyo3::exceptions::PyValueError;
use pyo3::{exceptions::PyTypeError, prelude::*, types::PyTuple};

use crate::atom::AtomView;
use crate::periodic_box::PeriodicBoxPy;
use crate::topology_state::{StatePy, TopologyPy};
use crate::utils::*;
use crate::SelPy;

#[pyclass(name = "System", frozen)]
pub struct SystemPy {
    // Since we have to replace top and st we need to wrap them into UnsafeCell
    pub(crate) top: UnsafeCell<Py<TopologyPy>>,
    pub(crate) st: UnsafeCell<Py<StatePy>>,
}

unsafe impl Send for SystemPy {}
unsafe impl Sync for SystemPy {}

impl SystemPy {
    pub(crate) fn top(&self) -> &Topology {
        unsafe{&*self.top.get()}.get().inner()
    }

    pub(crate) fn top_mut(&self) -> &mut Topology {
        unsafe{&*self.top.get()}.get().inner_mut()
    }

    pub(crate) fn st(&self) -> &State {
        unsafe{&*self.st.get()}.get().inner()
    }

    pub(crate) fn st_mut(&self) -> &mut State {
        unsafe{&*self.st.get()}.get().inner_mut()
    }

    pub(crate) fn top_py(&self) -> &Py<TopologyPy> {
        unsafe{&*self.top.get()}
    }

    pub(crate) fn top_py_mut(&self) -> &mut Py<TopologyPy> {
        unsafe{&mut *self.top.get()}
    }

    pub(crate) fn st_py(&self) -> &Py<StatePy> {
        unsafe{&*self.st.get()}
    }

    pub(crate) fn st_py_mut(&self) -> &mut Py<StatePy> {
        unsafe{&mut *self.st.get()}
    }
}

impl LenProvider for SystemPy {
    fn len(&self) -> usize {
        self.top().len()
    }
}

impl IndexProvider for SystemPy {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        i
    }
}

impl AtomPosAnalysis for SystemPy {
    fn atoms_ptr(&self) -> *const Atom {
        self.top().atoms.as_ptr()
    }

    fn coords_ptr(&self) -> *const Pos {
        self.st().coords.as_ptr()
    }
}

impl AtomPosAnalysisMut for SystemPy {}

impl BoxProvider for SystemPy {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.st().get_box()
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
        self.st().get_time()
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
            let mut fh = molar::io::FileHandler::open(&fname).map_err(to_py_io_err)?;
            let (top, st) = fh.read().map_err(to_py_io_err)?;
            Ok(SystemPy {
                top: UnsafeCell::new(TopologyPy::from(top).into_py()),
                st: UnsafeCell::new(StatePy::from(st).into_py()),
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
                    top: UnsafeCell::new(top.clone().unbind()),
                    st: UnsafeCell::new(st.clone().unbind()),
                })
            }
        } else {
            // Empty System
            Ok(SystemPy {
                top: UnsafeCell::new(TopologyPy(Default::default()).into_py()),
                st: UnsafeCell::new(StatePy(Default::default()).into_py()),
            })
        }
    }

    fn __len__(&self) -> usize {
        self.top().len()
    }

    #[pyo3(signature = (arg=None))]
    fn __call__(slf: &Bound<Self>, arg: Option<&Bound<'_, PyAny>>) -> PyResult<SelPy> {
        let sys = slf.get();

        let index = if let Some(arg) = arg {
            // Argument present
            if let Ok(val) = arg.extract::<String>() {
                if val.is_empty() {
                    // Select all on empty string
                    (0..sys.len()).into_sel_index(&*sys, None)
                } else {
                    // Otherwise do normal textual selection
                    val.into_sel_index(&*sys, None)
                }
            } else if let Ok(val) = arg.extract::<(usize, usize)>() {
                // Range selection
                (val.0..val.1).into_sel_index(&*sys, None)
            } else if let Ok(val) = arg.extract::<Vec<usize>>() {
                // Vector of indices
                val.into_sel_index(&*sys, None)
            } else {
                let ty_name = arg.get_type().name()?.to_string();
                return Err(PyTypeError::new_err(format!(
                    "Invalid argument type {ty_name} when creating selection"
                )));
            }
        } else {
            // No argument, select all
            (0..sys.len()).into_sel_index(&*sys, None)
        };

        Ok(SelPy {
            top: UnsafeCell::new(sys.top_py().clone_ref(slf.py())),
            st: UnsafeCell::new(sys.st_py().clone_ref(slf.py())),
            index: index.map_err(|e| PyTypeError::new_err(e.to_string()))?,
        })
    }

    fn replace_state_deep(&self, st: &Bound<StatePy>) -> PyResult<()> {
        if self.st().interchangeable(st.get().inner()) {
            unsafe{std::ptr::swap(self.st_mut(), st.get().inner_mut())};
            Ok(())
        } else {
            return Err(PyValueError::new_err("incompatible state"));
        }
    }

    // fn replace_state_from(&mut self, arg: &Bound<'_, PyAny>) -> PyResult<StatePy> {
    //     if let Ok(sys) = arg.cast::<SystemPy>() {
    //         let st = sys.borrow().st.clone_ref();
    //         self.replace_state(&st)
    //     } else if let Ok(sel) = arg.cast::<SelPy>() {
    //         let st = sel.borrow().sys.st.clone_ref();
    //         self.replace_state(&st)
    //     } else {
    //         Err(PyTypeError::new_err(format!(
    //             "Invalid argument type {} in set_state_from()",
    //             arg.get_type()
    //         )))
    //     }
    // }

    // fn replace_topology(&mut self, top: &TopologyPy) -> PyResult<TopologyPy> {
    //     if self.top.inner().interchangeable(top.inner()) {
    //         let ret = self.top.clone_ref();
    //         self.top = top.clone_ref();
    //         Ok(ret)
    //     } else {
    //         return Err(PyValueError::new_err("incompatible topology"));
    //     }
    // }

    // fn replace_topology_deep(&mut self, top: &mut TopologyPy) -> PyResult<()> {
    //     if self.top.inner().interchangeable(top.inner()) {
    //         mem::swap(self.top.inner_mut(), top.inner_mut());
    //         Ok(())
    //     } else {
    //         return Err(PyValueError::new_err("incompatible topology"));
    //     }
    // }

    #[getter("state")]
    fn get_state(slf: Bound<Self>) -> Bound<StatePy> {
        slf.get().st_py().bind(slf.py()).clone()
    }

    #[setter("state")]
    fn set_state(&self, st: &Bound<StatePy>) -> PyResult<()> {
        if self.st().interchangeable(st.get().inner()) {
            *self.st_py_mut() = st.clone().unbind();
            Ok(())
        } else {
            return Err(PyValueError::new_err("incompatible state"));
        }
    }

    #[getter("topology")]
    fn get_topology(slf: Bound<Self>) -> Bound<TopologyPy> {
        slf.get().top_py().bind(slf.py()).clone()
    }

    #[setter("topology")]
    fn set_topology(&self, top: &Bound<TopologyPy>) -> PyResult<()> {
        if self.top().interchangeable(top.get().inner()) {
            *self.top_py_mut() = top.clone().unbind();
            Ok(())
        } else {
            return Err(PyValueError::new_err("incompatible topology"));
        }
    }

    fn save(&self, fname: &str) -> PyResult<()> {
        Ok(SaveTopologyState::save(self, fname).map_err(to_py_io_err)?)
    }

    fn remove<'py>(slf: &Bound<'py, Self>, arg: &Bound<'py, PyAny>) -> PyResult<()> {
        if let Ok(sel) = arg.cast::<SelPy>() {
            // Selection provided
            let sb = sel.get();
            sb
                .top_mut()
                .remove_atoms(sb.iter_index())
                .map_err(to_py_runtime_err)?;
            sb.st_mut()
                .remove_coords(sb.iter_index())
                .map_err(to_py_runtime_err)?;
            Ok(())
        } else {
            let sel = Self::__call__(slf, Some(arg))?;
            sel
                .top_mut()
                .remove_atoms(sel.iter_index())
                .map_err(to_py_runtime_err)?;
            sel
                .st_mut()
                .remove_coords(sel.iter_index())
                .map_err(to_py_runtime_err)?;
            Ok(())
        }
    }

    #[pyo3(signature = (*args))]
    fn append<'py>(slf: &Bound<'py, Self>, args: &Bound<'py, PyTuple>) -> PyResult<()> {
        let slf_b = slf.get();

        if args.len() == 1 {
            let arg = args.get_item(0)?;
            let sel = if let Ok(sel) = arg.cast::<SelPy>() {
                sel
            } else {
                &Bound::new(slf.py(), Self::__call__(slf, Some(&arg))?)?
            };

            slf_b
                .top_mut()
                .add_atoms(sel.borrow().iter_atoms().cloned());
            slf_b
                .st_mut()
                .add_coords(sel.borrow().iter_pos().cloned());

            Ok(())
        } else if args.len() == 2 {
            let arg1 = args.get_item(0)?;
            let ab = arg1.cast::<crate::atom::AtomPy>()?.borrow();
            let pos = args.get_item(1)?.extract::<PyArrayLike1<f32>>()?;
            let v: VectorView<f32, Const<3>> = pos.try_as_matrix().unwrap();
            slf_b
                .top_mut()
                .add_atoms(std::iter::once(&ab.0).cloned());
            slf_b
                .st_mut()
                .add_coords(std::iter::once(Pos::new(v.x, v.y, v.z)));
            Ok(())
        } else {
            Err(PyValueError::new_err("1 or 2 arguments expected"))
        }
    }

    #[getter]
    fn get_time(&self) -> f32 {
        TimeProvider::get_time(self)
    }

    #[getter]
    fn get_box(&self) -> PeriodicBoxPy {
        self.st()
            .pbox
            .as_ref()
            .map(|b| PeriodicBoxPy(b.clone()))
            .unwrap()
    }

    #[setter]
    fn set_box(&self, b: &PeriodicBoxPy) {
        *self.st_mut().pbox.as_mut().unwrap() = b.0.clone();
    }

    #[setter]
    fn set_time(&self, t: f32) {
        self.st_mut().time = t;
    }

    fn set_box_from(&self, src: Bound<'_, PyAny>) -> PyResult<()> {
        let st_ref = if let Ok(sys) = src.cast::<SystemPy>() {
            sys.get().st()
        } else if let Ok(sel) = src.cast::<SelPy>() {
            sel.get().st()
        } else {
            return Err(PyTypeError::new_err(format!(
                "Invalid argument type {} in set_box_from()",
                src.get_type().name()?.to_string()
            )));
        };
        self.st_mut().pbox = st_ref.pbox.clone();
        Ok(())
    }

    fn iter_pos(slf: PyRef<'_, Self>) -> PyRef<'_, SysPosIterator> {
        Bound::new(
            slf.py(),
            SysPosIterator {
                st: slf.st_py().clone_ref(slf.py()),
                cur: 0,
            },
        )
        .unwrap()
        .borrow()
    }

    fn iter_atoms(slf: PyRef<'_, Self>) -> PyRef<'_, SysAtomIterator> {
        Bound::new(
            slf.py(),
            SysAtomIterator {
                top: slf.top_py().clone_ref(slf.py()),
                cur: 0,
            },
        )
        .unwrap()
        .borrow()
    }
}

#[pyclass]
pub struct SysPosIterator {
    pub(crate) st: Py<StatePy>,
    pub(crate) cur: usize,
}

#[pymethods]
impl SysPosIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__<'py>(slf: &Bound<'py, Self>) -> Option<Bound<'py, PyArray1<f32>>> {
        let mut s = slf.borrow_mut();
        if s.cur >= s.st.get().len() {
            return None;
        }
        s.cur += 1;
        unsafe { Some(map_pyarray_to_pos(&s.st.bind(slf.py()), s.cur)) }
    }
}

#[pyclass]
pub struct SysAtomIterator {
    pub(crate) top: Py<TopologyPy>,
    pub(crate) cur: usize,
}

#[pymethods]
impl SysAtomIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<AtomView> {
        if self.cur >= self.top.get().len() {
            return None;
        }
        let atom_ptr = unsafe { self.top.get().inner().atoms.as_ptr().add(self.cur) as *mut Atom };
        self.cur += 1;
        Some(AtomView(atom_ptr))
    }
}
