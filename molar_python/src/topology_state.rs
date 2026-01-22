use std::cell::UnsafeCell;

use crate::{SelPy, SystemPy};

use super::periodic_box::PeriodicBoxPy;
use molar::prelude::*;
use pyo3::{exceptions::{PyAttributeError, PyTypeError}, prelude::*};
use triomphe::Arc;

#[pyclass(name = "State")]
pub(crate) struct StatePy (pub(crate) Arc<UnsafeCell<State>>);

unsafe impl Send for StatePy {}
unsafe impl Sync for StatePy {}

impl From<State> for StatePy {
    fn from(value: State) -> Self {
        Self(Arc::new(UnsafeCell::new(value)))
    }
}

impl StatePy {
    pub(crate) fn get(&self) -> &State {
        unsafe {&*self.0.get()}
    }

    pub(crate) fn get_mut(&self) -> &mut State {
        unsafe {&mut *self.0.get()}
    }

    pub(crate) fn clone_ref(&self) -> StatePy {
        Self(Arc::clone(&self.0))
    }
}

#[pymethods]
impl StatePy {
    fn __len__(&self) -> usize {
        self.get().len()
    }

    #[getter]
    fn get_box(&self) -> PyResult<PeriodicBoxPy> {
        Ok(PeriodicBoxPy(
            self.get().pbox
                .as_ref()
                .ok_or_else(|| PyAttributeError::new_err("No periodic box to get"))?
                .clone(),
        ))
    }

    #[setter]
    fn set_box(&mut self, val: Bound<'_, PeriodicBoxPy>) -> PyResult<()> {
        let b = self.get_mut()
            .pbox
            .as_mut()
            .ok_or_else(|| PyAttributeError::new_err("No periodic box to set"))?;
        *b = val.borrow().0.clone();
        Ok(())
    }

    #[getter]
    fn get_time(&mut self) -> f32 {
        self.get().time
    }

    #[setter]
    fn set_time(&mut self, t: f32) {
        self.get_mut().time = t;
    }

    fn set_box_from(&mut self, arg: Bound<'_, PyAny>) -> PyResult<()> {
        let st_ref = if let Ok(sys) = arg.cast::<SystemPy>() {
            &sys.borrow().st
        } else if let Ok(sel) = arg.cast::<SelPy>() {
            &sel.borrow().sys.st
        } else {
            let ty_name = arg.get_type().name()?.to_string();
            return Err(PyTypeError::new_err(format!(
                "Invalid argument type {ty_name} in set_box_from()"
            )));
        };
        self.get_mut().pbox = st_ref.get().pbox.clone();
        Ok(())
    }
}

impl SaveState for StatePy {
    fn iter_pos_dyn<'a>(&'a self) -> Box<dyn Iterator<Item = &'a Pos> + 'a> {
        self.get().iter_pos_dyn()
    }
}

impl LenProvider for StatePy {
    fn len(&self) -> usize {
        self.__len__()
    }
}

impl TimeProvider for StatePy {
    fn get_time(&self) -> f32 {
        self.get().time
    }
}

impl BoxProvider for StatePy {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.get().pbox.as_ref()
    }
}

impl RandomPosProvider for StatePy {
    unsafe fn get_pos_unchecked(&self, i: usize) -> &Pos {
        self.get().get_pos_unchecked(i)
    }
}

//----------------------------------------------------------------

#[pyclass(name = "Topology")]
pub(crate) struct TopologyPy(pub(crate) Arc<UnsafeCell<Topology>>);

unsafe impl Send for TopologyPy {}
unsafe impl Sync for TopologyPy {}

impl From<Topology> for TopologyPy {
    fn from(value: Topology) -> Self {
        Self(Arc::new(UnsafeCell::new(value)))
    }
}

impl TopologyPy {
    pub(crate) fn get(&self) -> &Topology {
        unsafe {&*self.0.get()}
    }

    pub(crate) fn get_mut(&self) -> &mut Topology {
        unsafe {&mut *self.0.get()}
    }

    pub(crate) fn clone_ref(&self) -> TopologyPy {
        Self(Arc::clone(&self.0))
    }
}

impl LenProvider for TopologyPy {
    fn len(&self)-> usize {
        self.get().len()
    }
}

#[pymethods]
impl TopologyPy {
    fn __len__(&self) -> usize {
        self.get().len()
    }
}

impl RandomBondProvider for TopologyPy {
    fn num_bonds(&self) -> usize {
        self.get().num_bonds()
    }

    unsafe fn get_bond_unchecked(&self, i: usize) -> &[usize; 2] {
        self.get().get_bond_unchecked(i)
    }
}

impl SaveTopology for TopologyPy {
    fn iter_atoms_dyn<'a>(&'a self) -> Box<dyn Iterator<Item = &'a Atom> + 'a> {
        self.get().iter_atoms_dyn()
    }
}