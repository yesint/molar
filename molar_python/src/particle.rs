use molar::Float;
use super::topology_state::TopologyPy;
use crate::atom::AtomView;
use crate::utils::map_pyarray_to_pos;
use crate::{atom::AtomPy, topology_state::StatePy};
use molar::{AtomLike, AtomMutProvider, State, Topology};
use numpy::{PyArray1, PyArrayLike1, PyArrayMethods, PyUntypedArrayMethods};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
/// View over one atom and its coordinates inside a system/selection.
///
/// Exposes both coordinate fields and atom descriptors as mutable properties.
///
/// **Example**
///
/// .. code-block:: python
///
///    p = sys[0]
///    print(p.name, p.resid, p.pos)    # "CA", 1, [x, y, z]
///    p.x = 1.5    # move to 1.5 nm on x axis
///    p.mass = 12.0

#[pyclass(name = "Particle", frozen)]
pub(crate) struct ParticlePy {
    pub(crate) top: Py<TopologyPy>,
    pub(crate) st: Py<StatePy>,
    /// Global atom index within the parent system/selection (read-only).
    ///
    /// :returns: Global atom index.
    /// :rtype: int
    #[pyo3(get)]
    pub(crate) id: usize,
}

impl ParticlePy {
    pub(crate) fn top(&self) -> &Topology {
        self.top.get().inner()
    }

    pub(crate) fn top_mut(&self) -> &mut Topology {
        self.top.get().inner_mut()
    }

    pub(crate) fn st(&self) -> &State {
        self.st.get().inner()
    }

    pub(crate) fn st_mut(&self) -> &mut State {
        self.st.get().inner_mut()
    }
}

#[pymethods]
impl ParticlePy {
    /// Return atom position as a length-3 NumPy array view.
    ///
    /// :returns: Position vector ``[x, y, z]``.
    /// :rtype: numpy.ndarray
    #[getter(pos)]
    fn get_pos<'py>(slf: &'py Bound<'py, Self>) -> Bound<'py, PyArray1<Float>> {
        let s = slf.get();
        unsafe {
            map_pyarray_to_pos(s.st.bind(slf.py()), s.id)
        }
    }

    /// Set atom position from a length-3 vector.
    ///
    /// :param pos: New position vector.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter(pos)]
    fn set_pos(&self, pos: PyArrayLike1<Float>) -> PyResult<()> {
        if pos.len() != 3 {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "pos must have 3 elements",
            ));
        }
        let src = pos.data();
        let dst = self.st_mut().coords.as_mut_ptr() as *mut Float;
        if src != dst {
            unsafe { std::ptr::copy_nonoverlapping(src, dst, 3) };
        }
        Ok(())
    }

    /// X coordinate.
    ///
    /// :returns: X coordinate.
    /// :rtype: float
    #[getter(x)]
    fn get_x(&self) -> Float {
        unsafe { self.st().coords.get_unchecked(self.id).x }
    }

    /// Set X coordinate.
    ///
    /// :param value: New X coordinate.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter(x)]
    fn set_x(&self, value: Float) {
        unsafe { self.st_mut().coords.get_unchecked_mut(self.id).x = value }
    }

    /// Y coordinate.
    ///
    /// :returns: Y coordinate.
    /// :rtype: float
    #[getter(y)]
    fn get_y(&self) -> Float {
        unsafe { self.st().coords.get_unchecked(self.id).y }
    }

    /// Set Y coordinate.
    ///
    /// :param value: New Y coordinate.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter(y)]
    fn set_y(&self, value: Float) {
        unsafe { self.st_mut().coords.get_unchecked_mut(self.id).y = value }
    }

    /// Z coordinate.
    ///
    /// :returns: Z coordinate.
    /// :rtype: float
    #[getter(z)]
    fn get_z(&self) -> Float {
        unsafe { self.st().coords.get_unchecked(self.id).z }
    }

    /// Set Z coordinate.
    ///
    /// :param value: New Z coordinate.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter(z)]
    fn set_z(&self, value: Float) {
        unsafe { self.st_mut().coords.get_unchecked_mut(self.id).z = value }
    }

    fn __repr__(&self) -> String {
        let pos = unsafe { self.st().coords.get_unchecked(self.id) };
        let a = unsafe { self.top().atoms.get_unchecked(self.id) };
        format!(
            "Particle(id={}, name='{}', resname='{}', resid={}, pos=[{:.3}, {:.3}, {:.3}])",
            self.id, a.name, a.resname, a.resid, pos.x, pos.y, pos.z
        )
    }

    //atom
    /// Get mutable atom view.
    ///
    /// :returns: Mutable atom view.
    /// :rtype: Atom
    #[getter(atom)]
    fn get_atom(slf: &Bound<'_, Self>) -> AtomView {
        let s = slf.get();
        AtomView { top: s.top.clone_ref(slf.py()), index: s.id }
    }

    /// Replace atom from ``Atom`` or ``AtomView``.
    ///
    /// :param arg: Source atom object.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter(atom)]
    fn set_atom(&self, arg: &Bound<'_, PyAny>) -> PyResult<()> {
        let at = if let Ok(at) = arg.cast::<AtomPy>() {
            at.borrow().0.clone()
        } else if let Ok(v) = arg.cast::<AtomView>() {
            v.borrow().atom()?.clone()
        } else {
            let ty_name = arg.get_type().name()?.to_string();
            return Err(PyTypeError::new_err(format!(
                "Invalid argument type {ty_name} in set_atom()"
            )));
        };
        unsafe { *self.top_mut().get_atom_mut_unchecked(self.id) = at };
        Ok(())
    }

    /// Atom name.
    ///
    /// :returns: Atom name.
    /// :rtype: str
    #[getter(name)]
    fn get_name(&self) -> String {
        unsafe {
            self.top()
                .atoms
                .get_unchecked(self.id)
                .name
                .as_str()
                .to_owned()
        }
    }

    /// Set atom name.
    ///
    /// :param value: New atom name.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter(name)]
    fn set_name(&self, value: &str) {
        unsafe {
            self.top_mut().atoms.get_unchecked_mut(self.id).set_name(value)
        }
    }
    // resname
    /// Residue name.
    ///
    /// :returns: Residue name.
    /// :rtype: str
    #[getter(resname)]
    fn get_resname(&self) -> String {
        unsafe {
            self.top()
                .atoms
                .get_unchecked(self.id)
                .resname
                .as_str()
                .to_owned()
        }
    }

    /// Set residue name.
    ///
    /// :param value: New residue name.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter(resname)]
    fn set_resname(&self, value: &str) {
        unsafe {
            self.top_mut().atoms.get_unchecked_mut(self.id).set_resname(value)
        }
    }

    // resid
    /// Residue identifier.
    ///
    /// :returns: Residue identifier.
    /// :rtype: int
    #[getter(resid)]
    fn get_resid(&self) -> i32 {
        unsafe { self.top().atoms.get_unchecked(self.id).resid }
    }

    /// Set residue identifier.
    ///
    /// :param value: New residue identifier.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter(resid)]
    fn set_resid(&self, value: i32) {
        unsafe { self.top_mut().atoms.get_unchecked_mut(self.id).resid = value }
    }

    // resindex
    /// Residue index.
    ///
    /// :returns: Residue index.
    /// :rtype: int
    #[getter(resindex)]
    fn get_resindex(&self) -> usize {
        unsafe { self.top().atoms.get_unchecked(self.id).resindex }
    }

    /// Set residue index.
    ///
    /// :param value: New residue index.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter(resindex)]
    fn set_resindex(&self, value: usize) {
        unsafe { self.top_mut().atoms.get_unchecked_mut(self.id).resindex = value }
    }

    // atomic_number
    /// Atomic number.
    ///
    /// :returns: Atomic number.
    /// :rtype: int
    #[getter(atomic_number)]
    fn get_atomic_number(&self) -> u8 {
        unsafe { self.top().atoms.get_unchecked(self.id).atomic_number }
    }

    /// Set atomic number.
    ///
    /// :param value: New atomic number.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter(atomic_number)]
    fn set_atomic_number(&self, value: u8) {
        unsafe {
            self.top_mut()
                .atoms
                .get_unchecked_mut(self.id)
                .atomic_number = value
        }
    }

    // mass
    /// Atomic mass.
    ///
    /// :returns: Atomic mass.
    /// :rtype: float
    #[getter(mass)]
    fn get_mass(&self) -> Float {
        unsafe { self.top().atoms.get_unchecked(self.id).mass }
    }

    /// Set atomic mass.
    ///
    /// :param value: New atomic mass.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter(mass)]
    fn set_mass(&self, value: Float) {
        unsafe { self.top_mut().atoms.get_unchecked_mut(self.id).mass = value }
    }

    // charge
    /// Atom charge.
    ///
    /// :returns: Atom charge.
    /// :rtype: float
    #[getter(charge)]
    fn get_charge(&self) -> Float {
        unsafe { self.top().atoms.get_unchecked(self.id).charge }
    }

    /// Set atom charge.
    ///
    /// :param value: New atom charge.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter(charge)]
    fn set_charge(&self, value: Float) {
        unsafe { self.top_mut().atoms.get_unchecked_mut(self.id).charge = value }
    }

    // type_name
    /// Force-field atom type name.
    ///
    /// :returns: Force-field type name.
    /// :rtype: str
    #[getter(type_name)]
    fn get_type_name(&self) -> String {
        unsafe {
            self.top()
                .atoms
                .get_unchecked(self.id)
                .type_name
                .as_str()
                .to_owned()
        }
    }

    /// Set force-field atom type name.
    ///
    /// :param value: New force-field type name.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter(type_name)]
    fn set_type_name(&self, value: &str) {
        unsafe {
            self.top_mut()
                .atoms
                .get_unchecked_mut(self.id)
                .set_type_name(value)
        }
    }

    // type_id
    /// Force-field atom type identifier.
    ///
    /// :returns: Force-field type id.
    /// :rtype: int
    #[getter(type_id)]
    fn get_type_id(&self) -> u32 {
        unsafe { self.top().atoms.get_unchecked(self.id).type_id }
    }

    /// Set force-field atom type identifier.
    ///
    /// :param value: New force-field type id.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter(type_id)]
    fn set_type_id(&self, value: u32) {
        unsafe { self.top_mut().atoms.get_unchecked_mut(self.id).type_id = value }
    }

    // chain
    /// Chain identifier.
    ///
    /// :returns: Chain identifier.
    /// :rtype: str
    #[getter(chain)]
    fn get_chain(&self) -> char {
        unsafe { self.top().atoms.get_unchecked(self.id).chain }
    }

    /// Set chain identifier.
    ///
    /// :param value: New chain identifier.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter(chain)]
    fn set_chain(&self, value: char) {
        unsafe { self.top_mut().atoms.get_unchecked_mut(self.id).chain = value }
    }

    // bfactor
    /// Temperature factor (B-factor).
    ///
    /// :returns: B-factor value.
    /// :rtype: float
    #[getter(bfactor)]
    fn get_bfactor(&self) -> Float {
        unsafe { self.top().atoms.get_unchecked(self.id).bfactor }
    }

    /// Set temperature factor (B-factor).
    ///
    /// :param value: New B-factor.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter(bfactor)]
    fn set_bfactor(&self, value: Float) {
        unsafe { self.top_mut().atoms.get_unchecked_mut(self.id).bfactor = value }
    }

    // occupancy
    /// Occupancy value.
    ///
    /// :returns: Occupancy value.
    /// :rtype: float
    #[getter(occupancy)]
    fn get_occupancy(&self) -> Float {
        unsafe { self.top().atoms.get_unchecked(self.id).occupancy }
    }

    /// Set occupancy value.
    ///
    /// :param value: New occupancy value.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter(occupancy)]
    fn set_occupancy(&self, value: Float) {
        unsafe {
            self.top_mut()
                .atoms
                .get_unchecked_mut(self.id)
                .occupancy = value
        }
    }
}
