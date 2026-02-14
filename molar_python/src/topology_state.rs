use std::cell::UnsafeCell;

use crate::{SelPy, SystemPy};

use super::periodic_box::PeriodicBoxPy;
use molar::prelude::*;
use pyo3::{
    exceptions::{PyAttributeError, PyTypeError},
    prelude::*,
};
/// Coordinate frame with simulation time and periodic box.
///
/// Stores coordinates, optional periodic box, and time for one frame.

#[pyclass(name = "State", frozen)]
pub struct StatePy(pub(crate) UnsafeCell<State>);

unsafe impl Send for StatePy {}
unsafe impl Sync for StatePy {}

impl From<State> for StatePy {
    fn from(value: State) -> Self {
        Self(UnsafeCell::new(value))
    }
}

impl StatePy {
    pub(crate) fn inner(&self) -> &State {
        unsafe { &*self.0.get() }
    }

    pub(crate) fn inner_mut(&self) -> &mut State {
        unsafe { &mut *self.0.get() }
    }

    pub(crate) fn into_py(self) -> Py<StatePy> {
        Python::attach(|py| Py::new(py, self).unwrap())
    }
}

#[pymethods]
impl StatePy {
    /// Number of coordinates in this state.
    ///
    /// :returns: Number of atoms/coordinates.
    /// :rtype: int
    fn __len__(&self) -> usize {
        self.inner().len()
    }

    /// Periodic box of this state.
    ///
    /// :returns: Periodic box.
    /// :rtype: PeriodicBox
    #[getter]
    fn get_box(&self) -> PyResult<PeriodicBoxPy> {
        Ok(PeriodicBoxPy(
            self.inner()
                .pbox
                .as_ref()
                .ok_or_else(|| PyAttributeError::new_err("No periodic box to get"))?
                .clone(),
        ))
    }

    /// Set periodic box of this state.
    ///
    /// :param val: New periodic box.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter]
    fn set_box(&self, val: Bound<'_, PeriodicBoxPy>) -> PyResult<()> {
        let b = self
            .inner_mut()
            .pbox
            .as_mut()
            .ok_or_else(|| PyAttributeError::new_err("No periodic box to set"))?;
        *b = val.borrow().0.clone();
        Ok(())
    }

    /// Simulation time value.
    ///
    /// :returns: Frame time.
    /// :rtype: float
    #[getter]
    fn get_time(&self) -> f32 {
        self.inner().time
    }

    /// Set simulation time value.
    ///
    /// :param t: New frame time.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter]
    fn set_time(&self, t: f32) {
        self.inner_mut().time = t;
    }

    /// Copy periodic box from a ``System`` or ``Sel``.
    ///
    /// :param arg: Source object (``System`` or ``Sel``).
    /// :returns: ``None``.
    /// :rtype: None
    fn set_box_from(&self, arg: Bound<'_, PyAny>) -> PyResult<()> {
        let st_ref = if let Ok(sys) = arg.cast::<SystemPy>() {
            sys.get().r_st()
        } else if let Ok(sel) = arg.cast::<SelPy>() {
            sel.get().r_st()
        } else {
            let ty_name = arg.get_type().name()?.to_string();
            return Err(PyTypeError::new_err(format!(
                "Invalid argument type {ty_name} in set_box_from()"
            )));
        };
        self.inner_mut().pbox = st_ref.pbox.clone();
        Ok(())
    }
}

impl SaveState for StatePy {
    fn iter_pos_dyn<'a>(&'a self) -> Box<dyn Iterator<Item = &'a Pos> + 'a> {
        self.inner().iter_pos_dyn()
    }
}

impl LenProvider for StatePy {
    fn len(&self) -> usize {
        self.__len__()
    }
}

impl TimeProvider for StatePy {
    fn get_time(&self) -> f32 {
        self.inner().time
    }
}

impl BoxProvider for StatePy {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.inner().pbox.as_ref()
    }
}

impl RandomPosProvider for StatePy {
    unsafe fn get_pos_unchecked(&self, i: usize) -> &Pos {
        self.inner().get_pos_unchecked(i)
    }
}

//----------------------------------------------------------------
/// Molecular topology container (atoms and connectivity).

#[pyclass(name = "Topology", frozen)]
pub struct TopologyPy(pub(crate) UnsafeCell<Topology>);

unsafe impl Send for TopologyPy {}
unsafe impl Sync for TopologyPy {}

impl From<Topology> for TopologyPy {
    fn from(value: Topology) -> Self {
        Self(UnsafeCell::new(value))
    }
}

impl TopologyPy {
    pub(crate) fn inner(&self) -> &Topology {
        unsafe { &*self.0.get() }
    }

    pub(crate) fn inner_mut(&self) -> &mut Topology {
        unsafe { &mut *self.0.get() }
    }

    
    pub(crate) fn into_py(self) -> Py<TopologyPy> {
        Python::attach(|py| Py::new(py, self).unwrap())
    }
}

impl LenProvider for TopologyPy {
    fn len(&self) -> usize {
        self.inner().len()
    }
}

#[pymethods]
impl TopologyPy {
    /// Number of atoms in this topology.
    fn __len__(&self) -> usize {
        self.inner().len()
    }
}

impl RandomBondProvider for TopologyPy {
    fn num_bonds(&self) -> usize {
        self.inner().num_bonds()
    }

    unsafe fn get_bond_unchecked(&self, i: usize) -> &[usize; 2] {
        self.inner().get_bond_unchecked(i)
    }
}

impl SaveTopology for TopologyPy {
    fn iter_atoms_dyn<'a>(&'a self) -> Box<dyn Iterator<Item = &'a Atom> + 'a> {
        self.inner().iter_atoms_dyn()
    }
}
