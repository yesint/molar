use super::topology_state::TopologyPy;
use crate::atom::AtomView;
use crate::utils::map_pyarray_to_pos;
use crate::{atom::AtomPy, topology_state::StatePy};
use molar::RandomAtomMutProvider;
use numpy::{PyArray1, PyArrayLike1, PyArrayMethods, PyUntypedArrayMethods};
use pyo3::exceptions::PyTypeError;
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
    #[getter(pos)]
    fn get_pos<'py>(slf: &'py Bound<'py, Self>) -> Bound<'py, PyArray1<f32>> {
        let s = slf.borrow();
        unsafe {
            map_pyarray_to_pos(&s.st, s.id, slf.py())
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
        let dst = self.st.inner_mut().coords.as_mut_ptr() as *mut f32;
        if src != dst {
            unsafe { std::ptr::copy_nonoverlapping(src, dst, 3) };
        }
        Ok(())
    }

    #[getter(x)]
    fn get_x(&self) -> f32 {
        unsafe { self.st.inner().coords.get_unchecked(self.id).x }
    }

    #[setter(x)]
    fn set_x(&mut self, value: f32) {
        unsafe { self.st.inner_mut().coords.get_unchecked_mut(self.id).x = value }
    }

    #[getter(y)]
    fn get_y(&self) -> f32 {
        unsafe { self.st.inner().coords.get_unchecked(self.id).y }
    }

    #[setter(y)]
    fn set_y(&mut self, value: f32) {
        unsafe { self.st.inner_mut().coords.get_unchecked_mut(self.id).y = value }
    }

    #[getter(z)]
    fn get_z(&self) -> f32 {
        unsafe { self.st.inner().coords.get_unchecked(self.id).z }
    }

    #[setter(z)]
    fn set_z(&mut self, value: f32) {
        unsafe { self.st.inner_mut().coords.get_unchecked_mut(self.id).z = value }
    }

    //atom
    #[getter(atom)]
    fn get_atom(&self) -> AtomView {
        unsafe { AtomView(self.top.inner_mut().get_atom_mut_unchecked(self.id)) }
    }

    #[setter(atom)]
    fn set_atom(&mut self, arg: &Bound<'_, PyAny>) -> PyResult<()> {
        let at = if let Ok(at) = arg.cast::<AtomPy>() {
            at.borrow().0.clone()
        } else if let Ok(v) = arg.cast::<AtomView>() {
            unsafe { (*v.borrow().0).clone() }
        } else {
            let ty_name = arg.get_type().name()?.to_string();
            return Err(PyTypeError::new_err(format!(
                "Invalid argument type {ty_name} in set_atom()"
            )));
        };
        unsafe { *self.top.inner_mut().get_atom_mut_unchecked(self.id) = at };
        Ok(())
    }

    #[getter(name)]
    fn get_name(&self) -> String {
        unsafe {
            self.top
                .inner()
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
            self.top.inner_mut().atoms.get_unchecked_mut(self.id).name = value.to_owned().into()
        }
    }
    // resname
    #[getter(resname)]
    fn get_resname(&self) -> String {
        unsafe {
            self.top
                .inner()
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
            self.top.inner_mut().atoms.get_unchecked_mut(self.id).resname = value.to_owned().into()
        }
    }

    // resid
    #[getter(resid)]
    fn get_resid(&self) -> i32 {
        unsafe { self.top.inner().atoms.get_unchecked(self.id).resid }
    }

    #[setter(resid)]
    fn set_resid(&mut self, value: i32) {
        unsafe { self.top.inner_mut().atoms.get_unchecked_mut(self.id).resid = value }
    }

    // resindex
    #[getter(resindex)]
    fn get_resindex(&self) -> usize {
        unsafe { self.top.inner().atoms.get_unchecked(self.id).resindex }
    }

    #[setter(resindex)]
    fn set_resindex(&mut self, value: usize) {
        unsafe { self.top.inner_mut().atoms.get_unchecked_mut(self.id).resindex = value }
    }

    // atomic_number
    #[getter(atomic_number)]
    fn get_atomic_number(&self) -> u8 {
        unsafe { self.top.inner().atoms.get_unchecked(self.id).atomic_number }
    }

    #[setter(atomic_number)]
    fn set_atomic_number(&mut self, value: u8) {
        unsafe {
            self.top
                .inner_mut()
                .atoms
                .get_unchecked_mut(self.id)
                .atomic_number = value
        }
    }

    // mass
    #[getter(mass)]
    fn get_mass(&self) -> f32 {
        unsafe { self.top.inner().atoms.get_unchecked(self.id).mass }
    }

    #[setter(mass)]
    fn set_mass(&mut self, value: f32) {
        unsafe { self.top.inner_mut().atoms.get_unchecked_mut(self.id).mass = value }
    }

    // charge
    #[getter(charge)]
    fn get_charge(&self) -> f32 {
        unsafe { self.top.inner().atoms.get_unchecked(self.id).charge }
    }

    #[setter(charge)]
    fn set_charge(&mut self, value: f32) {
        unsafe { self.top.inner_mut().atoms.get_unchecked_mut(self.id).charge = value }
    }

    // type_name
    #[getter(type_name)]
    fn get_type_name(&self) -> String {
        unsafe {
            self.top
                .inner()
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
                .inner_mut()
                .atoms
                .get_unchecked_mut(self.id)
                .type_name = value.to_owned().into()
        }
    }

    // type_id
    #[getter(type_id)]
    fn get_type_id(&self) -> u32 {
        unsafe { self.top.inner().atoms.get_unchecked(self.id).type_id }
    }

    #[setter(type_id)]
    fn set_type_id(&mut self, value: u32) {
        unsafe { self.top.inner_mut().atoms.get_unchecked_mut(self.id).type_id = value }
    }

    // chain
    #[getter(chain)]
    fn get_chain(&self) -> char {
        unsafe { self.top.inner().atoms.get_unchecked(self.id).chain }
    }

    #[setter(chain)]
    fn set_chain(&mut self, value: char) {
        unsafe { self.top.inner_mut().atoms.get_unchecked_mut(self.id).chain = value }
    }

    // bfactor
    #[getter(bfactor)]
    fn get_bfactor(&self) -> f32 {
        unsafe { self.top.inner().atoms.get_unchecked(self.id).bfactor }
    }

    #[setter(bfactor)]
    fn set_bfactor(&mut self, value: f32) {
        unsafe { self.top.inner_mut().atoms.get_unchecked_mut(self.id).bfactor = value }
    }

    // occupancy
    #[getter(occupancy)]
    fn get_occupancy(&self) -> f32 {
        unsafe { self.top.inner().atoms.get_unchecked(self.id).occupancy }
    }

    #[setter(occupancy)]
    fn set_occupancy(&mut self, value: f32) {
        unsafe {
            self.top
                .inner_mut()
                .atoms
                .get_unchecked_mut(self.id)
                .occupancy = value
        }
    }
}
