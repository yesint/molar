use crate::{map_const_pyarray_to_vec3, Sel, System};
use pyo3::prelude::*;
use triomphe::Arc;

#[pyclass(unsendable)]
pub(crate) struct Membrane(molar_membrane::Membrane);

#[pymethods]
impl Membrane {
    #[new]
    fn new(src: &System, opt_str: &str) -> anyhow::Result<Self> {
        Ok(Self(molar_membrane::Membrane::new(&src.0, opt_str)?))
    }

    fn add_ids_to_group(&mut self, gr_name: &str, ids: Vec<usize>) -> anyhow::Result<()> {
        Ok(self.0.add_ids_to_group(gr_name, &ids)?)
    }

    fn add_resindeces_to_group(&mut self, gr_name: &str, res: Vec<usize>) -> anyhow::Result<()> {
        Ok(self.0.add_resindeces_to_group(gr_name, &res)?)
    }

    fn set_state(&mut self, st: &crate::State) -> anyhow::Result<()> {
        self.0.set_state(Arc::clone(&st.0))?;
        Ok(())
    }

    fn compute(&mut self) -> PyResult<()> {
        Ok(self.0.compute()?)
    }

    fn finalize(&mut self) -> PyResult<()> {
        Ok(self.0.finalize()?)
    }

    fn smooth_curvature(&mut self, n_neib: usize) {
        self.0.smooth_curvature(n_neib)
    }

    #[getter]
    fn get_all_lipids(&self) -> Vec<LipidMolecule> {
        type M = molar_membrane::LipidMolecule;
        self.0
            .iter_all_lipids()
            .map(|l| LipidMolecule(l as *const M as *mut M))
            .collect()
    }

    #[getter]
    fn get_valid_lipids(&self) -> Vec<LipidMolecule> {
        type M = molar_membrane::LipidMolecule;
        self.0
            .iter_valid_lipids()
            .map(|l| LipidMolecule(l as *const M as *mut M))
            .collect()
    }

    fn group_lipids(&self, gr_name: &str) -> anyhow::Result<Vec<LipidMolecule>> {
        type M = molar_membrane::LipidMolecule;
        Ok(self
            .0
            .iter_group(gr_name)?
            .map(|l| LipidMolecule(l as *const M as *mut M))
            .collect())
    }

    fn group_ids(&self, gr_name: &str) -> anyhow::Result<Vec<usize>> {
        Ok(self.0.iter_group_ids(gr_name)?.collect())
    }

    fn group_lipids_valid(&self, gr_name: &str) -> anyhow::Result<Vec<LipidMolecule>> {
        type M = molar_membrane::LipidMolecule;
        Ok(self
            .0
            .iter_group_valid(gr_name)?
            .map(|l| LipidMolecule(l as *const M as *mut M))
            .collect())
    }

    fn group_ids_valid(&self, gr_name: &str) -> anyhow::Result<Vec<usize>> {
        Ok(self.0.iter_group_ids_valid(gr_name)?.collect())
    }

    fn reset_groups(&mut self) {
        self.0.reset_groups();
    }

    fn reset_valid_lipids(&mut self) {
        self.0.reset_valid_lipids();
    }

    fn write_vmd_visualization(&self, fname: &str) -> anyhow::Result<()> {
        self.0.write_vmd_visualization(fname)
    }
}

#[pyclass(unsendable)]
pub(crate) struct LipidMolecule(*mut molar_membrane::LipidMolecule);

impl LipidMolecule {
    #[inline(always)]
    fn get(&self) -> &molar_membrane::LipidMolecule {
        unsafe { &*self.0 }
    }

    #[inline(always)]
    fn get_mut(&mut self) -> &mut molar_membrane::LipidMolecule {
        unsafe { &mut *self.0 }
    }
}

#[pymethods]
impl LipidMolecule {
    #[getter]
    fn get_id(&self) -> usize {
        self.get().id
    }

    #[getter]
    fn get_valid(&self) -> bool {
        self.get().valid
    }

    #[setter]
    fn set_valid(&mut self, val: bool) {
        self.get_mut().valid = val;
    }

    #[getter]
    fn get_sel(&self) -> Sel {
        unsafe { Sel::new_ref(std::mem::transmute(&self.get().sel)) }
    }

    #[getter]
    fn get_head_sel(&self) -> Sel {
        unsafe { Sel::new_ref(std::mem::transmute(&self.get().head_sel)) }
    }

    #[getter]
    fn get_mid_sel(&self) -> Sel {
        unsafe { Sel::new_ref(std::mem::transmute(&self.get().mid_sel)) }
    }

    #[getter]
    fn get_tail_end_sel(&self) -> Sel {
        unsafe { Sel::new_ref(std::mem::transmute(&self.get().tail_end_sel)) }
    }

    #[getter]
    fn get_head_marker(slf: Bound<'_,Self>) -> Bound<'_, PyAny> {
        let ptr = map_const_pyarray_to_vec3(slf.py(), &slf.borrow().get().head_marker.coords, &slf);
        unsafe { Bound::from_owned_ptr(slf.py(), ptr.cast()) }
    }

    #[getter]
    fn get_mid_marker(slf: Bound<'_,Self>) -> Bound<'_, PyAny> {
        let ptr = map_const_pyarray_to_vec3(slf.py(), &slf.borrow().get().mid_marker.coords, &slf);
        unsafe { Bound::from_owned_ptr(slf.py(), ptr.cast()) }
    }

    #[getter]
    fn get_tail_marker(slf: Bound<'_,Self>) -> Bound<'_, PyAny> {
        let ptr = map_const_pyarray_to_vec3(slf.py(), &slf.borrow().get().tail_marker.coords, &slf);
        unsafe { Bound::from_owned_ptr(slf.py(), ptr.cast()) }
    }

    #[getter]
    fn get_area(&self) -> f32 {
        self.get().area
    }

    #[getter]
    fn get_mean_curv(&self) -> f32 {
        self.get().mean_curv
    }

    #[getter]
    fn get_gauss_curv(&self) -> f32 {
        self.get().gaussian_curv
    }

    #[getter]
    fn get_normal(slf: Bound<'_,Self>) -> Bound<'_, PyAny> {
        let ptr = map_const_pyarray_to_vec3(slf.py(), &slf.borrow().get().normal, &slf);
        unsafe { Bound::from_owned_ptr(slf.py(), ptr.cast()) }
    }
}

#[pyclass]
pub(crate) struct Histogram1D(molar_membrane::Histogram1D);

#[pymethods]
impl Histogram1D {
    #[new]
    fn new(min: f32, max: f32, n_bins: usize) -> Self {
        Self(molar_membrane::Histogram1D::new(min, max, n_bins))
    }

    pub fn add_one(&mut self, val: f32) {
        self.0.add_one(val);
    }

    pub fn normalize_density(&mut self) {
        self.0.normalize_density();
    }

    pub fn save_to_file(&self, fname: &str) -> anyhow::Result<()> {
        self.0.save_to_file(fname)
    }
}
