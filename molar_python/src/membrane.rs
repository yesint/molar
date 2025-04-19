use crate::{map_const_pyarray_to_pos, Sel, Source, State};
use anyhow::anyhow;
use molar::core::{BoxProvider, BuilderSerial, Holder};
use pyo3::prelude::*;

#[pyclass(unsendable)]
pub(crate) struct Membrane(molar_membrane::Membrane<BuilderSerial>);

#[pymethods]
impl Membrane {
    #[new]
    fn new(src: &Source, opt_str: &str) -> anyhow::Result<Self> {
        Ok(Self(molar_membrane::Membrane::new(&src.0, opt_str)?))
    }

    fn add_lipids_to_group(&mut self, gr_name: &str, ids: Vec<usize>) -> anyhow::Result<()> {
        Ok(self.0.add_lipids_to_group(gr_name, &ids)?)
    }

    fn set_state<'py>(&mut self, st: &Bound<'py, crate::State>) -> PyResult<()> {
        let mut st_ref = &st.borrow_mut().0;
        todo!("we don't need to release state here, just use existing holder");
        Ok(())
    }

    fn compute(&mut self) -> PyResult<()> {
        Ok(self.0.compute()?)
    }

    fn finalize(&mut self) -> PyResult<()> {
        Ok(self.0.finalize()?)
    }

    #[getter]
    fn get_lipids(&self) -> Vec<LipidMolecule> {
        type M = molar_membrane::LipidMolecule<BuilderSerial>;
        self.0
            .iter_lipids()
            .map(|(id, l)| LipidMolecule(l as *const M as *mut M, id))
            .collect()
    }
}

#[pyclass(unsendable)]
pub(crate) struct LipidMolecule(
    *mut molar_membrane::LipidMolecule<BuilderSerial>,
    usize,
);

impl LipidMolecule {
    // fn get_mut(&mut self) -> &mut molar_membrane::LipidMolecule<BuilderSerial> {
    //     unsafe {&mut *self.0}
    // }

    fn get(&self) -> &molar_membrane::LipidMolecule<BuilderSerial> {
        unsafe {&*self.0}
    }
}


#[pymethods]
impl LipidMolecule {
    #[getter]
    fn get_id(&self) -> usize {
        self.1
    }

    #[getter]
    fn get_sel(&self) -> Sel {
        Sel::new_ref(&self.get().sel)
    }

    #[getter]
    fn get_head_sel(&self) -> Sel {
        Sel::new_ref(&self.get().head_sel)
    }

    #[getter]
    fn get_mid_sel(&self) -> Sel {
        Sel::new_ref(&self.get().mid_sel)
    }

    #[getter]
    fn get_tail_end_sel(&self) -> Sel {
        Sel::new_ref(&self.get().tail_end_sel)
    }

    #[getter]
    fn get_head_marker(slf: Bound<Self>) -> Bound<'_, PyAny> {
        let ptr = map_const_pyarray_to_pos(slf.py(), &slf.borrow().get().head_marker, &slf);
        unsafe { Bound::from_owned_ptr(slf.py(), ptr.cast()) }
    }

    #[getter]
    fn get_mid_marker(slf: Bound<Self>) -> Bound<'_, PyAny> {
        let ptr = map_const_pyarray_to_pos(slf.py(), &slf.borrow().get().mid_marker, &slf);
        unsafe { Bound::from_owned_ptr(slf.py(), ptr.cast()) }
    }

    #[getter]
    fn get_tail_marker(slf: Bound<Self>) -> Bound<'_, PyAny> {
        let ptr = map_const_pyarray_to_pos(slf.py(), &slf.borrow().get().tail_marker, &slf);
        unsafe { Bound::from_owned_ptr(slf.py(), ptr.cast()) }
    }
}
