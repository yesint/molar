use super::topology_state::TopologyPy;
use crate::topology_state::StatePy;
use crate::utils::map_pyarray_to_pos;
use numpy::{PyArray1, PyArrayLike1, PyArrayMethods, PyUntypedArrayMethods};
use pyo3::prelude::*;

#[pyclass(name = "Particle")]
pub(crate) struct ParticlePy {
    pub(crate) top: TopologyPy,
    pub(crate) st: StatePy,
    // id is readonly
    #[pyo3(get)]
    pub(crate) id: usize,
}

#[pymethods]
impl ParticlePy {
    // //pos
    // #[getter]
    // fn get_pos<'py>(slf: &Bound<'py, Self>) -> Bound<'py, PyAny> {
    //     let s = slf.borrow();
    //     let mut st = s.st.borrow_mut(slf.py());
    //     let v = st.0.get_pos_mut(s.id).unwrap();
    //     let p = super::utils::map_pyarray_to_pos(v, slf);
    //     unsafe{Bound::from_borrowed_ptr(slf.py(), p as *mut pyo3::ffi::PyObject)}
    // }

    // #[setter]
    // fn set_pos<'py>(&mut self, py: Python<'py>, value: PyReadonlyArray1<'py, f32>) {
    //     let mut p = self.st.borrow_mut(py).0.get_pos_mut(self.id).unwrap().coords;
    //     unsafe {
    //         std::ptr::copy_nonoverlapping(value.data(), p.as_mut_ptr(), 3);
    //     }
    // }

    fn get_id(&self) -> usize {
        self.id
    }

    #[getter(pos)]
    fn get_pos<'py>(slf: &'py Bound<'py, Self>) -> Bound<'py, PyArray1<f32>> {
        let s = slf.borrow();
        unsafe {
            let data_ptr = s.st.get_mut().coords.as_mut_ptr().add(s.id);
            map_pyarray_to_pos(data_ptr, slf)
        }
    }

    #[setter(pos)]
    fn set_pos(&mut self, pos: PyArrayLike1<f32>) -> PyResult<()> {
        if pos.len() != 3 {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "pos must have 3 elements",
            ));
        }
        let src = pos.data();
        let dst = self.st.get_mut().coords.as_mut_ptr() as *mut f32;
        if src != dst {
            unsafe { std::ptr::copy_nonoverlapping(src, dst, 3) };
        }
        Ok(())
    }

    #[getter(x)]
    fn get_x(&self) -> f32 {
        unsafe { self.st.get().coords.get_unchecked(self.id).x }
    }

    #[setter(x)]
    fn set_x(&mut self, value: f32) {
        unsafe { self.st.get_mut().coords.get_unchecked_mut(self.id).x = value }
    }

    #[getter(y)]
    fn get_y(&self) -> f32 {
        unsafe { self.st.get().coords.get_unchecked(self.id).y }
    }

    #[setter(y)]
    fn set_y(&mut self, value: f32) {
        unsafe { self.st.get_mut().coords.get_unchecked_mut(self.id).y = value }
    }

    #[getter(z)]
    fn get_z(&self) -> f32 {
        unsafe { self.st.get().coords.get_unchecked(self.id).z }
    }

    #[setter(z)]
    fn set_z(&mut self, value: f32) {
        unsafe { self.st.get_mut().coords.get_unchecked_mut(self.id).z = value }
    }

    // atom
    // #[getter(atom)]
    // fn get_atom(&self) -> AtomPy {
    //     unsafe{ AtomPy(self.top.get_mut().atoms.as_mut_ptr().add(self.id)) }
    // }

    // #[setter(atom)]
    // fn set_atom(&mut self, py: Python<'_>, value: &AtomPy) {
    //     let mut top = self.top.borrow_mut(py);
    //     *top.0.get_atom_mut(self.id).unwrap() = value.0.clone();
    // }

    #[getter(name)]
    fn get_name(&self) -> String {
        unsafe {
            self.top
                .get()
                .atoms
                .get_unchecked(self.id)
                .name
                .as_str()
                .to_owned()
        }
    }

    #[setter(name)]
    fn set_name(&mut self, value: &str) {
        unsafe {
            self.top.get_mut().atoms.get_unchecked_mut(self.id).name = value.to_owned().into()
        }
    }
    // resname
    #[getter(resname)]
    fn get_resname(&self) -> String {
        unsafe {
            self.top
                .get()
                .atoms
                .get_unchecked(self.id)
                .resname
                .as_str()
                .to_owned()
        }
    }

    #[setter(resname)]
    fn set_resname(&mut self, value: &str) {
        unsafe {
            self.top.get_mut().atoms.get_unchecked_mut(self.id).resname = value.to_owned().into()
        }
    }

    // resid
    #[getter(resid)]
    fn get_resid(&self) -> i32 {
        unsafe { self.top.get().atoms.get_unchecked(self.id).resid }
    }

    #[setter(resid)]
    fn set_resid(&mut self, value: i32) {
        unsafe { self.top.get_mut().atoms.get_unchecked_mut(self.id).resid = value }
    }

    // resindex
    #[getter(resindex)]
    fn get_resindex(&self) -> usize {
        unsafe { self.top.get().atoms.get_unchecked(self.id).resindex }
    }

    #[setter(resindex)]
    fn set_resindex(&mut self, value: usize) {
        unsafe { self.top.get_mut().atoms.get_unchecked_mut(self.id).resindex = value }
    }

    // atomic_number
    #[getter(atomic_number)]
    fn get_atomic_number(&self) -> u8 {
        unsafe { self.top.get().atoms.get_unchecked(self.id).atomic_number }
    }

    #[setter(atomic_number)]
    fn set_atomic_number(&mut self, value: u8) {
        unsafe {
            self.top
                .get_mut()
                .atoms
                .get_unchecked_mut(self.id)
                .atomic_number = value
        }
    }

    // mass
    #[getter(mass)]
    fn get_mass(&self) -> f32 {
        unsafe { self.top.get().atoms.get_unchecked(self.id).mass }
    }

    #[setter(mass)]
    fn set_mass(&mut self, value: f32) {
        unsafe { self.top.get_mut().atoms.get_unchecked_mut(self.id).mass = value }
    }

    // charge
    #[getter(charge)]
    fn get_charge(&self) -> f32 {
        unsafe { self.top.get().atoms.get_unchecked(self.id).charge }
    }

    #[setter(charge)]
    fn set_charge(&mut self, value: f32) {
        unsafe { self.top.get_mut().atoms.get_unchecked_mut(self.id).charge = value }
    }

    // type_name
    #[getter(type_name)]
    fn get_type_name(&self) -> String {
        unsafe {
            self.top
                .get()
                .atoms
                .get_unchecked(self.id)
                .type_name
                .as_str()
                .to_owned()
        }
    }

    #[setter(type_name)]
    fn set_type_name(&mut self, value: &str) {
        unsafe {
            self.top
                .get_mut()
                .atoms
                .get_unchecked_mut(self.id)
                .type_name = value.to_owned().into()
        }
    }

    // type_id
    #[getter(type_id)]
    fn get_type_id(&self) -> u32 {
        unsafe { self.top.get().atoms.get_unchecked(self.id).type_id }
    }

    #[setter(type_id)]
    fn set_type_id(&mut self, value: u32) {
        unsafe { self.top.get_mut().atoms.get_unchecked_mut(self.id).type_id = value }
    }

    // chain
    #[getter(chain)]
    fn get_chain(&self) -> char {
        unsafe { self.top.get().atoms.get_unchecked(self.id).chain }
    }

    #[setter(chain)]
    fn set_chain(&mut self, value: char) {
        unsafe { self.top.get_mut().atoms.get_unchecked_mut(self.id).chain = value }
    }

    // bfactor
    #[getter(bfactor)]
    fn get_bfactor(&self) -> f32 {
        unsafe { self.top.get().atoms.get_unchecked(self.id).bfactor }
    }

    #[setter(bfactor)]
    fn set_bfactor(&mut self, value: f32) {
        unsafe { self.top.get_mut().atoms.get_unchecked_mut(self.id).bfactor = value }
    }

    // occupancy
    #[getter(occupancy)]
    fn get_occupancy(&self) -> f32 {
        unsafe { self.top.get().atoms.get_unchecked(self.id).occupancy }
    }

    #[setter(occupancy)]
    fn set_occupancy(&mut self, value: f32) {
        unsafe {
            self.top
                .get_mut()
                .atoms
                .get_unchecked_mut(self.id)
                .occupancy = value
        }
    }
}
