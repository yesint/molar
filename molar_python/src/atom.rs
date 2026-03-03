use crate::topology_state::TopologyPy;
use molar::prelude::{Atom, AtomLike};
use pyo3::{exceptions::PyIndexError, prelude::*};
/// Mutable atom container.
///
/// **Example**
///
/// .. code-block:: python
///
///    import pymolar
///    a = pymolar.Atom()
///    a.name    = "CA"
///    a.resname = "ALA"
///    a.resid   = 1
///    a.chain   = 'A'
///    a.mass    = 12.011

#[pyclass(name = "Atom")]
pub(crate) struct AtomPy(pub(crate) Atom);

#[pymethods]
impl AtomPy {
    #[new]
    /// Create an atom with default fields.
    fn new() -> Self {
        Self {
            0: Default::default(),
        }
    }

    /// Atom name.
    #[getter(name)]
    fn get_name(&self) -> &str {
        self.0.name.as_str()
    }

    #[setter(name)]
    fn set_name(&mut self, value: &str) {
        self.0.set_name(value);
    }

    // resname
    /// Residue name.
    #[getter(resname)]
    fn get_resname(&self) -> &str {
        self.0.resname.as_str()
    }

    #[setter(resname)]
    fn set_resname(&mut self, value: &str) {
        self.0.set_resname(value);
    }

    // resid
    /// Residue identifier.
    #[getter(resid)]
    fn get_resid(&self) -> i32 {
        self.0.resid
    }

    #[setter(resid)]
    fn set_resid(&mut self, value: i32) {
        self.0.resid = value;
    }

    // atomic_number
    /// Atomic number.
    #[getter(atomic_number)]
    fn get_atomic_number(&self) -> u8 {
        self.0.atomic_number
    }

    #[setter(atomic_number)]
    fn set_atomic_number(&mut self, value: u8) {
        self.0.atomic_number = value;
    }

    // mass
    /// Atomic mass.
    #[getter(mass)]
    fn get_mass(&self) -> f32 {
        self.0.mass
    }

    #[setter(mass)]
    fn set_mass(&mut self, value: f32) {
        self.0.mass = value;
    }

    // charge
    /// Atom charge.
    #[getter(charge)]
    fn get_charge(&self) -> f32 {
        self.0.charge
    }

    #[setter(charge)]
    fn set_charge(&mut self, value: f32) {
        self.0.charge = value;
    }

    // type_name
    /// Force-field atom type name.
    #[getter(type_name)]
    fn get_type_name(&self) -> &str {
        self.0.type_name.as_str()
    }

    #[setter(type_name)]
    fn set_type_name(&mut self, value: &str) {
        self.0.set_type_name(value);
    }

    // type_id
    /// Force-field atom type identifier.
    #[getter(type_id)]
    fn get_type_id(&self) -> u32 {
        self.0.type_id
    }

    #[setter(type_id)]
    fn set_type_id(&mut self, value: u32) {
        self.0.type_id = value;
    }

    // chain
    /// Chain identifier.
    #[getter(chain)]
    fn get_chain(&self) -> char {
        self.0.chain
    }

    #[setter(chain)]
    fn set_chain(&mut self, value: char) {
        self.0.chain = value;
    }

    // bfactor
    /// Temperature factor (B-factor).
    #[getter(bfactor)]
    fn get_bfactor(&self) -> f32 {
        self.0.bfactor
    }

    #[setter(bfactor)]
    fn set_bfactor(&mut self, value: f32) {
        self.0.bfactor = value;
    }

    // occupancy
    /// Occupancy value.
    #[getter(occupancy)]
    fn get_occupancy(&self) -> f32 {
        self.0.occupancy
    }

    #[setter(occupancy)]
    fn set_occupancy(&mut self, value: f32) {
        self.0.occupancy = value;
    }

    fn __repr__(&self) -> String {
        format!(
            "Atom(name='{}', resname='{}', resid={}, chain='{}')",
            self.0.name, self.0.resname, self.0.resid, self.0.chain
        )
    }
}

//----------------------------------------

/// View into an atom stored inside a System or Sel.
///
/// Holds a reference-counted handle to the topology and the atom's index,
/// so it remains safe even if the originating System is reallocated by
/// ``append()``. After ``remove()`` the index may address a different atom
/// (by design), but will raise ``IndexError`` if out of range.
///
/// **Example**
///
/// .. code-block:: python
///
///    for atom_view in sys.iter_atoms():
///        print(atom_view.name, atom_view.resid)
#[pyclass(name = "AtomView", unsendable)]
pub(crate) struct AtomView {
    pub(crate) top: Py<TopologyPy>,
    pub(crate) index: usize,
}

impl AtomView {
    pub(crate) fn atom(&self) -> PyResult<&Atom> {
        self.top
            .get()
            .inner()
            .atoms
            .get(self.index)
            .ok_or_else(|| PyIndexError::new_err(format!("atom index {} out of range", self.index)))
    }

    fn atom_mut(&self) -> PyResult<&mut Atom> {
        self.top
            .get()
            .inner_mut()
            .atoms
            .get_mut(self.index)
            .ok_or_else(|| PyIndexError::new_err(format!("atom index {} out of range", self.index)))
    }
}

#[pymethods]
impl AtomView {
    /// Atom name (e.g. ``'CA'``).
    ///
    /// :returns: Atom name.
    /// :rtype: str
    #[getter(name)]
    fn get_name(&self) -> PyResult<&str> {
        Ok(self.atom()?.name.as_str())
    }

    /// Set atom name.
    ///
    /// :param value: New atom name.
    #[setter(name)]
    fn set_name(&self, value: &str) -> PyResult<()> {
        self.atom_mut()?.set_name(value);
        Ok(())
    }

    /// Residue name (e.g. ``'ALA'``).
    ///
    /// :returns: Residue name.
    /// :rtype: str
    #[getter(resname)]
    fn get_resname(&self) -> PyResult<&str> {
        Ok(self.atom()?.resname.as_str())
    }

    /// Set residue name.
    ///
    /// :param value: New residue name.
    #[setter(resname)]
    fn set_resname(&self, value: &str) -> PyResult<()> {
        self.atom_mut()?.set_resname(value);
        Ok(())
    }

    /// Residue sequence number.
    ///
    /// :returns: Residue ID.
    /// :rtype: int
    #[getter(resid)]
    fn get_resid(&self) -> PyResult<i32> {
        Ok(self.atom()?.resid)
    }

    /// Set residue sequence number.
    ///
    /// :param value: New residue ID.
    #[setter(resid)]
    fn set_resid(&self, value: i32) -> PyResult<()> {
        self.atom_mut()?.resid = value;
        Ok(())
    }

    /// Atomic number (e.g. 6 for carbon).
    ///
    /// :returns: Atomic number.
    /// :rtype: int
    #[getter(atomic_number)]
    fn get_atomic_number(&self) -> PyResult<u8> {
        Ok(self.atom()?.atomic_number)
    }

    /// Set atomic number.
    ///
    /// :param value: New atomic number.
    #[setter(atomic_number)]
    fn set_atomic_number(&self, value: u8) -> PyResult<()> {
        self.atom_mut()?.atomic_number = value;
        Ok(())
    }

    /// Atomic mass in Da.
    ///
    /// :returns: Atomic mass.
    /// :rtype: float
    #[getter(mass)]
    fn get_mass(&self) -> PyResult<f32> {
        Ok(self.atom()?.mass)
    }

    /// Set atomic mass in Da.
    ///
    /// :param value: New atomic mass.
    #[setter(mass)]
    fn set_mass(&self, value: f32) -> PyResult<()> {
        self.atom_mut()?.mass = value;
        Ok(())
    }

    /// Partial charge in elementary charge units.
    ///
    /// :returns: Partial charge.
    /// :rtype: float
    #[getter(charge)]
    fn get_charge(&self) -> PyResult<f32> {
        Ok(self.atom()?.charge)
    }

    /// Set partial charge in elementary charge units.
    ///
    /// :param value: New partial charge.
    #[setter(charge)]
    fn set_charge(&self, value: f32) -> PyResult<()> {
        self.atom_mut()?.charge = value;
        Ok(())
    }

    /// Force-field atom type name.
    ///
    /// :returns: Type name.
    /// :rtype: str
    #[getter(type_name)]
    fn get_type_name(&self) -> PyResult<&str> {
        Ok(self.atom()?.type_name.as_str())
    }

    /// Set force-field atom type name.
    ///
    /// :param value: New type name.
    #[setter(type_name)]
    fn set_type_name(&self, value: &str) -> PyResult<()> {
        self.atom_mut()?.set_type_name(value);
        Ok(())
    }

    /// Force-field atom type integer ID.
    ///
    /// :returns: Type ID.
    /// :rtype: int
    #[getter(type_id)]
    fn get_type_id(&self) -> PyResult<u32> {
        Ok(self.atom()?.type_id)
    }

    /// Set force-field atom type integer ID.
    ///
    /// :param value: New type ID.
    #[setter(type_id)]
    fn set_type_id(&self, value: u32) -> PyResult<()> {
        self.atom_mut()?.type_id = value;
        Ok(())
    }

    /// Single-character chain identifier.
    ///
    /// :returns: Chain character.
    /// :rtype: str
    #[getter(chain)]
    fn get_chain(&self) -> PyResult<char> {
        Ok(self.atom()?.chain)
    }

    /// Set chain identifier.
    ///
    /// :param value: New chain character.
    #[setter(chain)]
    fn set_chain(&self, value: char) -> PyResult<()> {
        self.atom_mut()?.chain = value;
        Ok(())
    }

    /// B-factor (temperature factor).
    ///
    /// :returns: B-factor value.
    /// :rtype: float
    #[getter(bfactor)]
    fn get_bfactor(&self) -> PyResult<f32> {
        Ok(self.atom()?.bfactor)
    }

    /// Set B-factor (temperature factor).
    ///
    /// :param value: New B-factor.
    #[setter(bfactor)]
    fn set_bfactor(&self, value: f32) -> PyResult<()> {
        self.atom_mut()?.bfactor = value;
        Ok(())
    }

    /// Crystallographic occupancy [0, 1].
    ///
    /// :returns: Occupancy value.
    /// :rtype: float
    #[getter(occupancy)]
    fn get_occupancy(&self) -> PyResult<f32> {
        Ok(self.atom()?.occupancy)
    }

    /// Set crystallographic occupancy [0, 1].
    ///
    /// :param value: New occupancy.
    #[setter(occupancy)]
    fn set_occupancy(&self, value: f32) -> PyResult<()> {
        self.atom_mut()?.occupancy = value;
        Ok(())
    }

    fn __repr__(&self) -> String {
        match self.atom() {
            Ok(a) => format!(
                "AtomView(name='{}', resname='{}', resid={}, chain='{}', index={})",
                a.name, a.resname, a.resid, a.chain, self.index
            ),
            Err(_) => format!("AtomView(<index {} out of range>)", self.index),
        }
    }
}
