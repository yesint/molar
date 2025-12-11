use molar::prelude::*;
use numpy::{PyArray1, PyArrayMethods, PyReadonlyArray1, ndarray::ArrayView1};
use pyo3::prelude::*;
use super::{atom::AtomPy, topology_state::{StatePy, TopologyPy}};

#[pyclass(unsendable, name="Particle")]
pub(crate) struct ParticlePy {
    pub(crate) top: Py<TopologyPy>,
    pub(crate) st: Py<StatePy>,
    // id is readonly
    #[pyo3(get)]
    pub(crate) id: usize,
}

#[pymethods]
impl ParticlePy {
    //pos
    #[getter]
    fn get_pos<'py>(slf: &Bound<'py, Self>) -> Bound<'py, PyAny> {
        let s = slf.borrow();
        let mut st = s.st.borrow_mut(slf.py());
        let v = st.0.get_pos_mut(s.id).unwrap();
        let p = super::utils::map_pyarray_to_pos(v, slf);
        unsafe{Bound::from_borrowed_ptr(slf.py(), p as *mut pyo3::ffi::PyObject)}
    }

    #[setter]
    fn set_pos<'py>(&mut self, py: Python<'py>, value: PyReadonlyArray1<'py, f32>) {
        let mut p = self.st.borrow_mut(py).0.get_pos_mut(self.id).unwrap().coords;
        unsafe {
            std::ptr::copy_nonoverlapping(value.data(), p.as_mut_ptr(), 3);
        }
    }

    fn get_id(&self) -> usize {
        self.id
    }

    #[getter(x)]
    fn get_x(&self, py: Python<'_>) -> f32 {
        let st = self.st.borrow(py);
        st.0.get_pos(self.id).unwrap().x
    }

    #[setter(x)]
    fn set_x(&mut self, py: Python<'_>, value: f32) {
        let mut st = self.st.borrow_mut(py);
        st.0.get_pos_mut(self.id).unwrap().x = value;
    }

    #[getter(y)]
    fn get_y(&self, py: Python<'_>) -> f32 {
        let st = self.st.borrow(py);
        st.0.get_pos(self.id).unwrap().y
    }

    #[setter(y)]
    fn set_y(&mut self, py: Python<'_>, value: f32) {
        let mut st = self.st.borrow_mut(py);
        st.0.get_pos_mut(self.id).unwrap().y = value;
    }

    #[getter(z)]
    fn get_z(&self, py: Python<'_>) -> f32 {
        let st = self.st.borrow(py);
        st.0.get_pos(self.id).unwrap().z
    }

    #[setter(z)]
    fn set_z(&mut self, py: Python<'_>, value: f32) {
        let mut st = self.st.borrow_mut(py);
        st.0.get_pos_mut(self.id).unwrap().z = value;
    }

    // atom
    #[getter(atom)]
    fn get_atom(&self, py: Python<'_>) -> AtomPy {
        let top = self.top.borrow(py);
        AtomPy(top.0.get_atom(self.id).unwrap().clone())
    }

    #[setter(atom)]
    fn set_atom(&mut self, py: Python<'_>, value: &AtomPy) {
        let mut top = self.top.borrow_mut(py);
        *top.0.get_atom_mut(self.id).unwrap() = value.0.clone();
    }

    // name
    #[getter(name)]
    fn get_name(&self, py: Python<'_>) -> String {
        let top = self.top.borrow(py);
        top.0.get_atom(self.id).unwrap().name.clone()
    }

    #[setter(name)]
    fn set_name(&mut self, py: Python<'_>, value: &str) {
        let mut top = self.top.borrow_mut(py);
        top.0.get_atom_mut(self.id).unwrap().name = value.into();
    }

    // resname
    #[getter(resname)]
    fn get_resname(&self, py: Python<'_>) -> String {
        let top = self.top.borrow(py);
        top.0.get_atom(self.id).unwrap().resname.clone()
    }

    #[setter(name)]
    fn set_resname(&mut self, py: Python<'_>, value: &str) {
        let mut top = self.top.borrow_mut(py);
        top.0.get_atom_mut(self.id).unwrap().resname = value.into();
    }

    // resid
    #[getter(resid)]
    fn get_resid(&self, py: Python<'_>) -> i32 {
        let top = self.top.borrow(py);
        top.0.get_atom(self.id).unwrap().resid
    }

    #[setter(resid)]
    fn set_resid(&mut self, py: Python<'_>, value: i32) {
        let mut top = self.top.borrow_mut(py);
        top.0.get_atom_mut(self.id).unwrap().resid = value;
    }

    // resindex
    #[getter(resindex)]
    fn get_resindex(&self, py: Python<'_>) -> usize {
        let top = self.top.borrow(py);
        top.0.get_atom(self.id).unwrap().resindex
    }

    #[setter(resindex)]
    fn set_resindex(&mut self, py: Python<'_>, value: usize) {
        let mut top = self.top.borrow_mut(py);
        top.0.get_atom_mut(self.id).unwrap().resindex = value;
    }

    // atomic_number
    #[getter(atomic_number)]
    fn get_atomic_number(&self, py: Python<'_>) -> u8 {
        let top = self.top.borrow(py);
        top.0.get_atom(self.id).unwrap().atomic_number
    }

    #[setter(atomic_number)]
    fn set_atomic_number(&mut self, py: Python<'_>, value: u8) {
        let mut top = self.top.borrow_mut(py);
        top.0.get_atom_mut(self.id).unwrap().atomic_number = value;
    }

    // mass
    #[getter(mass)]
    fn get_mass(&self, py: Python<'_>) -> f32 {
        let top = self.top.borrow(py);
        top.0.get_atom(self.id).unwrap().mass
    }

    #[setter(mass)]
    fn set_mass(&mut self, py: Python<'_>, value: f32) {
        let mut top = self.top.borrow_mut(py);
        top.0.get_atom_mut(self.id).unwrap().mass = value;
    }

    // charge
    #[getter(charge)]
    fn get_charge(&self, py: Python<'_>) -> f32 {
        let top = self.top.borrow(py);
        top.0.get_atom(self.id).unwrap().charge
    }

    #[setter(charge)]
    fn set_charge(&mut self, py: Python<'_>, value: f32) {
        let mut top = self.top.borrow_mut(py);
        top.0.get_atom_mut(self.id).unwrap().charge = value;
    }

    // type_name
    #[getter(type_name)]
    fn get_type_name(&self, py: Python<'_>) -> String {
        let top = self.top.borrow(py);
        top.0.get_atom(self.id).unwrap().type_name.clone()
    }

    #[setter(type_name)]
    fn set_type_name(&mut self, py: Python<'_>, value: &str) {
        let mut top = self.top.borrow_mut(py);
        top.0.get_atom_mut(self.id).unwrap().type_name = value.into();
    }

    // type_id
    #[getter(type_id)]
    fn get_type_id(&self, py: Python<'_>) -> u32 {
        let top = self.top.borrow(py);
        top.0.get_atom(self.id).unwrap().type_id
    }

    #[setter(type_id)]
    fn set_type_id(&mut self, py: Python<'_>, value: u32) {
        let mut top = self.top.borrow_mut(py);
        top.0.get_atom_mut(self.id).unwrap().type_id = value;
    }

    // chain
    #[getter(chain)]
    fn get_chain(&self, py: Python<'_>) -> char {
        let top = self.top.borrow(py);
        top.0.get_atom(self.id).unwrap().chain
    }

    #[setter(chain)]
    fn set_chain(&mut self, py: Python<'_>, value: char) {
        let mut top = self.top.borrow_mut(py);
        top.0.get_atom_mut(self.id).unwrap().chain = value;
    }

    // bfactor
    #[getter(bfactor)]
    fn get_bfactor(&self, py: Python<'_>) -> f32 {
        let top = self.top.borrow(py);
        top.0.get_atom(self.id).unwrap().bfactor
    }

    #[setter(bfactor)]
    fn set_bfactor(&mut self, py: Python<'_>, value: f32) {
        let mut top = self.top.borrow_mut(py);
        top.0.get_atom_mut(self.id).unwrap().bfactor = value;
    }

    // occupancy
    #[getter(occupancy)]
    fn get_occupancy(&self, py: Python<'_>) -> f32 {
        let top = self.top.borrow(py);
        top.0.get_atom(self.id).unwrap().occupancy
    }

    #[setter(occupancy)]
    fn set_occupancy(&mut self, py: Python<'_>, value: f32) {
        let mut top = self.top.borrow_mut(py);
        top.0.get_atom_mut(self.id).unwrap().occupancy = value;
    }
}
