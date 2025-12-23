use molar::prelude::Atom;
use pyo3::prelude::*;

#[pyclass(name = "Atom")]
pub(crate) struct AtomPy(pub(crate) Atom);

#[pymethods]
impl AtomPy {
    #[new]
    fn new() -> Self {
        Self {
            0: Default::default(),
        }
    }

    #[getter(name)]
    fn get_name(&self, _py: Python) -> &str {
        self.0.name.as_str()
    }

    #[setter(name)]
    fn set_name(&mut self, value: &str) {
        self.0.name = value.into();
    }

    // resname
    #[getter(resname)]
    fn get_resname(&self, _py: Python) -> &str {
        self.0.resname.as_str()
    }

    #[setter(resname)]
    fn set_resname(&mut self, value: &str) {
        self.0.resname = value.into();
    }

    // resid
    #[getter(resid)]
    fn get_resid(&self, _py: Python) -> i32 {
        self.0.resid
    }

    #[setter(resid)]
    fn set_resid(&mut self, value: i32) {
        self.0.resid = value;
    }

    // atomic_number
    #[getter(atomic_number)]
    fn get_atomic_number(&self, _py: Python) -> u8 {
        self.0.atomic_number
    }

    #[setter(atomic_number)]
    fn set_atomic_number(&mut self, value: u8) {
        self.0.atomic_number = value;
    }

    // mass
    #[getter(mass)]
    fn get_mass(&self, _py: Python) -> f32 {
        self.0.mass
    }

    #[setter(mass)]
    fn set_mass(&mut self, value: f32) {
        self.0.mass = value;
    }

    // charge
    #[getter(charge)]
    fn get_charge(&self, _py: Python) -> f32 {
        self.0.charge
    }

    #[setter(charge)]
    fn set_charge(&mut self, value: f32) {
        self.0.charge = value;
    }

    // type_name
    #[getter(type_name)]
    fn get_type_name(&self, _py: Python) -> &str {
        self.0.type_name.as_str()
    }

    #[setter(type_name)]
    fn set_type_name(&mut self, value: &str) {
        self.0.type_name = value.into();
    }

    // type_id
    #[getter(type_id)]
    fn get_type_id(&self, _py: Python) -> u32 {
        self.0.type_id
    }

    #[setter(type_id)]
    fn set_type_id(&mut self, value: u32) {
        self.0.type_id = value;
    }

    // chain
    #[getter(chain)]
    fn get_chain(&self, _py: Python) -> char {
        self.0.chain
    }

    #[setter(chain)]
    fn set_chain(&mut self, value: char) {
        self.0.chain = value;
    }

    // bfactor
    #[getter(bfactor)]
    fn get_bfactor(&self, _py: Python) -> f32 {
        self.0.bfactor
    }

    #[setter(bfactor)]
    fn set_bfactor(&mut self, value: f32) {
        self.0.bfactor = value;
    }

    // occupancy
    #[getter(occupancy)]
    fn get_occupancy(&self, _py: Python) -> f32 {
        self.0.occupancy
    }

    #[setter(occupancy)]
    fn set_occupancy(&mut self, value: f32) {
        self.0.occupancy = value;
    }
}

//----------------------------------------

#[pyclass(name = "AtomView", unsendable)]
pub(crate) struct AtomView(pub(crate) *mut Atom);

impl AtomView {
    pub(crate) fn new(atom: &mut Atom) -> Self {
        Self(atom)
    }
}

#[pymethods]
impl AtomView {
    #[getter(name)]
    fn get_name(&self, _py: Python) -> &str {
        unsafe { &*self.0 }.name.as_str()
    }

    #[setter(name)]
    fn set_name(&mut self, value: &str) {
        unsafe { &mut *self.0 }.name = value.into();
    }

    // resname
    #[getter(resname)]
    fn get_resname(&self, _py: Python) -> &str {
        unsafe { &*self.0 }.resname.as_str()
    }

    #[setter(resname)]
    fn set_resname(&mut self, value: &str) {
        unsafe { &mut *self.0 }.resname = value.into();
    }

    // resid
    #[getter(resid)]
    fn get_resid(&self, _py: Python) -> i32 {
        unsafe { &*self.0 }.resid
    }

    #[setter(resid)]
    fn set_resid(&mut self, value: i32) {
        unsafe { &mut *self.0 }.resid = value;
    }

    // atomic_number
    #[getter(atomic_number)]
    fn get_atomic_number(&self, _py: Python) -> u8 {
        unsafe { &*self.0 }.atomic_number
    }

    #[setter(atomic_number)]
    fn set_atomic_number(&mut self, value: u8) {
        unsafe { &mut *self.0 }.atomic_number = value;
    }

    // mass
    #[getter(mass)]
    fn get_mass(&self, _py: Python) -> f32 {
        unsafe { &*self.0 }.mass
    }

    #[setter(mass)]
    fn set_mass(&mut self, value: f32) {
        unsafe { &mut *self.0 }.mass = value;
    }

    // charge
    #[getter(charge)]
    fn get_charge(&self, _py: Python) -> f32 {
        unsafe { &*self.0 }.charge
    }

    #[setter(charge)]
    fn set_charge(&mut self, value: f32) {
        unsafe { &mut *self.0 }.charge = value;
    }

    // type_name
    #[getter(type_name)]
    fn get_type_name(&self, _py: Python) -> &str {
        unsafe { &*self.0 }.type_name.as_str()
    }

    #[setter(type_name)]
    fn set_type_name(&mut self, value: &str) {
        unsafe { &mut *self.0 }.type_name = value.into();
    }

    // type_id
    #[getter(type_id)]
    fn get_type_id(&self, _py: Python) -> u32 {
        unsafe { &*self.0 }.type_id
    }

    #[setter(type_id)]
    fn set_type_id(&mut self, value: u32) {
        unsafe { &mut *self.0 }.type_id = value;
    }

    // chain
    #[getter(chain)]
    fn get_chain(&self, _py: Python) -> char {
        unsafe { &*self.0 }.chain
    }

    #[setter(chain)]
    fn set_chain(&mut self, value: char) {
        unsafe { &mut *self.0 }.chain = value;
    }

    // bfactor
    #[getter(bfactor)]
    fn get_bfactor(&self, _py: Python) -> f32 {
        unsafe { &*self.0 }.bfactor
    }

    #[setter(bfactor)]
    fn set_bfactor(&mut self, value: f32) {
        unsafe { &mut *self.0 }.bfactor = value;
    }

    // occupancy
    #[getter(occupancy)]
    fn get_occupancy(&self, _py: Python) -> f32 {
        unsafe { &*self.0 }.occupancy
    }

    #[setter(occupancy)]
    fn set_occupancy(&mut self, value: f32) {
        unsafe { &mut *self.0 }.occupancy = value;
    }
}
