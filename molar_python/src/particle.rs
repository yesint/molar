use numpy::{npyffi, PyArrayMethods, PyReadonlyArray1};
use pyo3::prelude::*;

use crate::atom::Atom;

#[pyclass(unsendable)]
pub(crate) struct Particle {
    pub(crate) atom: &'static mut molar::core::Atom,
    // PyArray mapped to pos
    pub(crate) pos: *mut npyffi::PyArrayObject,
    // id is readonly
    #[pyo3(get)]
    pub(crate) id: usize,
}

#[pymethods]
impl Particle {
    //pos
    #[getter]
    fn get_pos<'py>(&self, py: Python<'py>) -> Bound<'py, PyAny> {
        // Return an owned ptr to avoid incorrect reference count.
        // I still don't quite understand why this is necessary.
        unsafe { Bound::from_owned_ptr(py, self.pos.cast()) }
    }

    #[setter]
    fn set_pos<'py>(&mut self, value: PyReadonlyArray1<'py, f32>) {
        unsafe {
            std::ptr::copy_nonoverlapping(value.data(), (*self.pos).data.cast(), 3);
        }
    }

    #[getter(x)]
    fn get_x(&self) -> f32 {
        unsafe { *(*self.pos).data.cast() }
    }

    #[setter(x)]
    fn set_x(&mut self, value: f32) {
        unsafe { *(*self.pos).data.cast() = value };
    }

    #[getter(y)]
    fn get_y(&self) -> f32 {
        unsafe { *(*self.pos).data.cast::<f32>().offset(1) }
    }

    #[setter(y)]
    fn set_y(&mut self, value: f32) {
        unsafe { *(*self.pos).data.cast::<f32>().offset(1) = value };
    }

    #[getter(z)]
    fn get_z(&self) -> f32 {
        unsafe { *(*self.pos).data.cast::<f32>().offset(2) }
    }

    #[setter(z)]
    fn set_z(&mut self, value: f32) {
        unsafe { *(*self.pos).data.cast::<f32>().offset(2) = value };
    }

    // atom
    #[getter(atom)]
    fn get_atom(&self, _py: Python) -> Atom {
        Atom(self.atom.clone())
    }

    #[setter(atom)]
    fn set_atom(&mut self, value: &Atom) {
        *self.atom = value.0.clone();
    }

    // name
    #[getter(name)]
    fn get_name(&self, _py: Python) -> &str {
        self.atom.name.as_str()
    }

    #[setter(name)]
    fn set_name(&mut self, value: &str) {
        self.atom.name = value.into();
    }

    // resname
    #[getter(resname)]
    fn get_resname(&self, _py: Python) -> &str {
        self.atom.resname.as_str()
    }

    #[setter(resname)]
    fn set_resname(&mut self, value: &str) {
        self.atom.resname = value.into();
    }

    // resid
    #[getter(resid)]
    fn get_resid(&self, _py: Python) -> i32 {
        self.atom.resid
    }

    #[setter(resid)]
    fn set_resid(&mut self, value: i32) {
        self.atom.resid = value;
    }

    // resindex
    #[getter(resindex)]
    fn get_resindex(&self, _py: Python) -> usize {
        self.id
    }

    #[setter(resindex)]
    fn set_resindex(&mut self, value: usize) {
        self.id = value;
    }

    // atomic_number
    #[getter(atomic_number)]
    fn get_atomic_number(&self, _py: Python) -> u8 {
        self.atom.atomic_number
    }

    #[setter(atomic_number)]
    fn set_atomic_number(&mut self, value: u8) {
        self.atom.atomic_number = value;
    }

    // mass
    #[getter(mass)]
    fn get_mass(&self, _py: Python) -> f32 {
        self.atom.mass
    }

    #[setter(mass)]
    fn set_mass(&mut self, value: f32) {
        self.atom.mass = value;
    }

    // charge
    #[getter(charge)]
    fn get_charge(&self, _py: Python) -> f32 {
        self.atom.charge
    }

    #[setter(charge)]
    fn set_charge(&mut self, value: f32) {
        self.atom.charge = value;
    }

    // type_name
    #[getter(type_name)]
    fn get_type_name(&self, _py: Python) -> &str {
        self.atom.type_name.as_str()
    }

    #[setter(type_name)]
    fn set_type_name(&mut self, value: &str) {
        self.atom.type_name = value.into();
    }

    // type_id
    #[getter(type_id)]
    fn get_type_id(&self, _py: Python) -> u32 {
        self.atom.type_id
    }

    #[setter(type_id)]
    fn set_type_id(&mut self, value: u32) {
        self.atom.type_id = value;
    }

    // chain
    #[getter(chain)]
    fn get_chain(&self, _py: Python) -> char {
        self.atom.chain
    }

    #[setter(chain)]
    fn set_chain(&mut self, value: char) {
        self.atom.chain = value;
    }

    // bfactor
    #[getter(bfactor)]
    fn get_bfactor(&self, _py: Python) -> f32 {
        self.atom.bfactor
    }

    #[setter(bfactor)]
    fn set_bfactor(&mut self, value: f32) {
        self.atom.bfactor = value;
    }

    // occupancy
    #[getter(occupancy)]
    fn get_occupancy(&self, _py: Python) -> f32 {
        self.atom.occupancy
    }

    #[setter(occupancy)]
    fn set_occupancy(&mut self, value: f32) {
        self.atom.occupancy = value;
    }
}
