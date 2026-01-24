use molar::prelude::*;
use numpy::PyArrayLike1;
use numpy::nalgebra::{Const, VectorView};
use pyo3::exceptions::PyValueError;
use pyo3::{exceptions::PyTypeError, prelude::*, types::PyTuple};

use crate::utils::*;
use crate::topology_state::{TopologyPy, StatePy};
use crate::periodic_box::PeriodicBoxPy;
use crate::SelPy;

#[pyclass(name = "System")]
pub struct SystemPy {
    pub top: TopologyPy,
    pub st: StatePy,
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
            let mut fh = molar::io::FileHandler::open(&fname).map_err(to_py_io_err)?;
            let (top, st) = fh.read().map_err(to_py_io_err)?;
            Ok(SystemPy {
                top: TopologyPy::from(top),
                st: StatePy::from(st),
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

    #[pyo3(signature = (arg=None))]
    fn __call__(slf: &Bound<Self>, arg: Option<&Bound<'_, PyAny>>) -> PyResult<SelPy> {
        let sys = slf.borrow();

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
            sys: slf.borrow().clone_ref(),
            index: index.map_err(|e| PyTypeError::new_err(e.to_string()))?,
        })
    }

    fn set_state(&mut self, st: &Bound<'_, StatePy>) -> PyResult<StatePy> {
        if self.st.get().interchangeable(st.borrow().get()) {
            let ret = self.st.clone_ref();
            self.st = st.borrow().clone_ref();
            Ok(ret)
        } else {
            Err(PyValueError::new_err("incompatible state"))
        }
    }

    fn set_topology(&mut self, top: &Bound<'_, TopologyPy>) -> PyResult<TopologyPy> {
        if self.top.get().interchangeable(top.borrow().get()) {
            let ret = self.top.clone_ref();
            self.top = top.borrow().clone_ref();
            Ok(ret)
        } else {
            Err(PyValueError::new_err("incompatible topology"))
        }
    }

    fn save(&self, fname: &str) -> PyResult<()> {
        Ok(SaveTopologyState::save(self, fname).map_err(to_py_io_err)?)
    }

    fn remove<'py>(slf: &Bound<'py, Self>, arg: &Bound<'py, PyAny>) -> PyResult<()> {
        if let Ok(sel) = arg.cast::<SelPy>() {
            // Selection provided
            let sb = sel.borrow();
            sb.sys
                .top
                .get_mut()
                .remove_atoms(sb.iter_index())
                .map_err(to_py_runtime_err)?;
            sb.sys
                .st
                .get_mut()
                .remove_coords(sb.iter_index())
                .map_err(to_py_runtime_err)?;
            Ok(())
        } else {
            let sel = Self::__call__(slf, Some(arg))?;
            sel.sys
                .top
                .get_mut()
                .remove_atoms(sel.iter_index())
                .map_err(to_py_runtime_err)?;
            sel.sys
                .st
                .get_mut()
                .remove_coords(sel.iter_index())
                .map_err(to_py_runtime_err)?;
            Ok(())
        }
    }

    #[pyo3(signature = (*args))]
    fn append<'py>(slf: &Bound<'py, Self>, args: &Bound<'py, PyTuple>) -> PyResult<()> {
        let slf_b = slf.borrow();

        if args.len() == 1 {
            let arg = args.get_item(0)?;
            let sel = if let Ok(sel) = arg.cast::<SelPy>() {
                sel
            } else {
                &Bound::new(slf.py(), Self::__call__(slf, Some(&arg))?)?
            };

            slf_b
                .top
                .get_mut()
                .add_atoms(sel.borrow().iter_atoms().cloned());
            slf_b
                .st
                .get_mut()
                .add_coords(sel.borrow().iter_pos().cloned());

            Ok(())
        } else if args.len() == 2 {
            let arg1 = args.get_item(0)?;
            let ab = arg1.cast::<crate::atom::AtomPy>()?.borrow();
            let pos = args.get_item(1)?.extract::<PyArrayLike1<f32>>()?;
            let v: VectorView<f32, Const<3>> = pos.try_as_matrix().unwrap();
            slf_b
                .top
                .get_mut()
                .add_atoms(std::iter::once(&ab.0).cloned());
            slf_b
                .st
                .get_mut()
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
        self.st
            .get()
            .pbox
            .as_ref()
            .map(|b| PeriodicBoxPy(b.clone()))
            .unwrap()
    }

    #[setter]
    fn set_box(&mut self, b: &PeriodicBoxPy) {
        *self.st.get_mut().pbox.as_mut().unwrap() = b.0.clone();
    }

    #[setter]
    fn set_time(&mut self, t: f32) {
        self.st.get_mut().time = t;
    }

    fn set_box_from(&mut self, src: Bound<'_, PyAny>) -> PyResult<()> {
        let st_ref = if let Ok(sys) = src.cast::<SystemPy>() {
            &sys.borrow().st
        } else if let Ok(sel) = src.cast::<SelPy>() {
            &sel.borrow().sys.st
        } else {
            return Err(PyTypeError::new_err(format!(
                "Invalid argument type {} in set_box_from()",
                src.get_type().name()?.to_string()
            )));
        };
        self.st.get_mut().pbox = st_ref.get().pbox.clone();
        Ok(())
    }
}
