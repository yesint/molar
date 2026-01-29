use molar::prelude::*;
use numpy::nalgebra::{Const, Dyn, VectorView};
use numpy::{
    PyArray1, PyArrayLike1, PyArrayMethods, PyReadonlyArray2, PyUntypedArrayMethods, ToPyArray,
};
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::{exceptions::PyIndexError, prelude::*, types::PyAny};

use crate::periodic_box::PeriodicBoxPy;
use crate::system::SystemPy;
use crate::topology_state::{StatePy, TopologyPy};
use crate::utils::*;
use crate::atom::AtomView;
use crate::{ParticleIterator, ParticlePy};

#[pyclass(sequence, name = "Sel")]
pub struct SelPy {
    pub sys: SystemPy,
    pub index: SVec,
}

impl LenProvider for SelPy {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SelPy {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        self.index.get_index_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> {
        self.index.iter_index()
    }
}

impl AtomPosAnalysis for SelPy {
    fn atoms_ptr(&self) -> *const Atom {
        self.sys.top.get().atoms.as_ptr()
    }

    fn coords_ptr(&self) -> *const Pos {
        self.sys.st.get().coords.as_ptr()
    }
}

impl AtomPosAnalysisMut for SelPy {}

impl BoxProvider for SelPy {
    fn get_box(&self) -> Option<&PeriodicBox> {
        let ptr = self
            .sys
            .st
            .get()
            .pbox
            .as_ref()
            .map(|b| b as *const PeriodicBox)?;
        unsafe { Some(&*ptr) }
    }
}

impl RandomBondProvider for SelPy {
    fn num_bonds(&self) -> usize {
        0
    }

    unsafe fn get_bond_unchecked(&self, _i: usize) -> &[usize; 2] {
        unreachable!()
    }
}

impl TimeProvider for SelPy {
    fn get_time(&self) -> f32 {
        self.sys.st.get().time
    }
}

impl SaveTopology for SelPy {
    fn iter_atoms_dyn<'a>(&'a self) -> Box<dyn Iterator<Item = &'a Atom> + 'a> {
        Box::new(self.iter_atoms())
    }
}

impl SaveState for SelPy {
    fn iter_pos_dyn<'a>(&'a self) -> Box<dyn Iterator<Item = &'a Pos> + 'a> {
        Box::new(self.iter_pos())
    }
}

impl SaveTopologyState for SelPy {}

impl MeasurePeriodic for SelPy {}

impl RandomMoleculeProvider for SelPy {
    unsafe fn get_molecule_unchecked(&self, i: usize) -> &[usize; 2] {
        self.sys.top.get().get_molecule_unchecked(i)
    }

    fn num_molecules(&self) -> usize {
        self.sys.top.get().num_molecules()
    }
}

impl SelPy {
    pub fn from_svec(&self, index: SVec) -> Self {
        Self {
            sys: self.sys.clone_ref(),
            index,
        }
    }
}

impl Clone for SelPy {
    fn clone(&self) -> Self {
        SelPy {
            sys: self.sys.clone_ref(),
            index: self.index.clone(),
        }
    }
}

#[pyclass]
pub struct SelPosIterator {
    pub(crate) sel: SelPy,
    pub(crate) cur: usize,
}

#[pymethods]
impl SelPosIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__<'py>(slf: &Bound<'py, Self>) -> Option<Bound<'py, PyArray1<f32>>> {
        let mut s = slf.borrow_mut();
        if s.cur >= s.sel.len() {
            return None;
        }
        let idx = unsafe { s.sel.index.get_index_unchecked(s.cur) };
        let pos_ptr = unsafe { s.sel.coords_ptr().add(idx) as *mut Pos };
        s.cur += 1;

        unsafe { Some(map_pyarray_to_pos(pos_ptr, slf)) }
    }
}

#[pyclass]
pub struct SelAtomIterator {
    pub(crate) sel: SelPy,
    pub(crate) cur: usize,
}

#[pymethods]
impl SelAtomIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<AtomView> {
        if self.cur >= self.sel.len() {
            return None;
        }
        let idx = unsafe { self.sel.index.get_index_unchecked(self.cur) };
        let atom_ptr = unsafe { self.sel.atoms_ptr().add(idx) as *mut Atom };
        self.cur += 1;
        Some(AtomView(atom_ptr))
    }
}

#[pymethods]
impl SelPy {
    fn __len__(&self) -> usize {
        self.index.len()
    }

    fn __call__(&self, arg: &Bound<'_, PyAny>) -> PyResult<SelPy> {
        let v = if let Ok(val) = arg.extract::<String>() {
            val.into_sel_index(&self.sys, Some(self.index.as_slice()))
                .map_err(to_py_runtime_err)?
        } else if let Ok(val) = arg.extract::<(usize, usize)>() {
            (val.0..=val.1)
                .into_sel_index(&self.sys, Some(self.index.as_slice()))
                .map_err(to_py_runtime_err)?
        } else if let Ok(val) = arg.extract::<Vec<usize>>() {
            val.into_sel_index(&self.sys, Some(self.index.as_slice()))
                .map_err(to_py_runtime_err)?
        } else {
            return Err(PyTypeError::new_err(format!(
                "Invalid argument type {} when creating selection",
                arg.get_type()
            )));
        };
        Ok(self.from_svec(v))
    }

    pub(crate) fn __getitem__(&self, i: isize) -> PyResult<ParticlePy> {
        let n = self.len();

        let mut ind = if i < 0 {
            if i.abs() > n as isize {
                return Err(PyIndexError::new_err(format!(
                    "Negative index {i} is out of bounds {}:-1",
                    -(n as isize)
                )));
            }
            n - i.unsigned_abs()
        } else if i >= n as isize {
            return Err(PyIndexError::new_err(format!(
                "Index {} is out of bounds 0:{}",
                i, n
            )));
        } else {
            i as usize
        };

        ind += unsafe { self.index.get_unchecked(0) };

        Ok(ParticlePy {
            top: self.sys.top.clone_ref(),
            st: self.sys.st.clone_ref(),
            id: ind,
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, ParticleIterator> {
        Bound::new(
            slf.py(),
            ParticleIterator {
                sel: slf.clone(),
                cur: 0,
            },
        )
        .unwrap()
        .borrow()
    }

    #[getter("index")]
    fn get_index<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray1<usize>> {
        numpy::PyArray1::from_iter(py, self.index.iter_index())
    }

    #[getter("coords")]
    fn get_coords<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f32>> {
        let coord_ptr = self.coords_ptr() as *const f32;
        unsafe {
            let arr = numpy::PyArray2::<f32>::new(py, [3, self.len()], true);
            let arr_ptr = arr.data();
            for i in self.index.iter_index() {
                let pos_ptr = coord_ptr.add(i * 3);
                std::ptr::copy_nonoverlapping(pos_ptr, arr_ptr.add(i * 3), 3);
            }
            arr
        }
    }

    #[getter("system")]
    fn get_system(&self) -> SystemPy {
        self.sys.clone_ref()
    }

    #[setter("system")]
    fn set_system(&mut self, sys: &SystemPy) -> PyResult<()> {
        if self.sys.st.get().interchangeable(sys.st.get()) {
            self.sys = sys.clone_ref();
            Ok(())
        } else {
            return Err(PyValueError::new_err("incompatible system"));
        }
    }

    #[setter("coords")]
    fn set_coords(&mut self, arr: PyReadonlyArray2<f32>) -> PyResult<()> {
        if arr.shape() != [3, self.__len__()] {
            return Err(PyValueError::new_err(format!(
                "Array shape must be [3, {}], not {:?}",
                self.__len__(),
                arr.shape()
            )));
        }
        let arr_ptr = arr.data();
        let coord_ptr = self.coords_ptr_mut() as *mut f32;

        unsafe {
            for i in self.index.iter_index() {
                let pos_ptr = coord_ptr.add(i * 3);
                std::ptr::copy_nonoverlapping(arr_ptr.add(i * 3), pos_ptr, 3);
            }
        }
        Ok(())
    }

    #[getter("state")]
    fn get_state(&mut self) -> StatePy {
        self.sys.st.clone_ref()
    }

    #[setter("state")]
    fn set_state(&mut self, st: &StatePy) -> PyResult<()> {
        if self.sys.st.get().interchangeable(st.get()) {
            self.sys.st = st.clone_ref();
            Ok(())
        } else {
            return Err(PyValueError::new_err("incompatible state"));
        }
    }

    fn replace_state(&mut self, st: &StatePy) -> PyResult<StatePy> {
        if self.sys.st.get().interchangeable(st.get()) {
            let ret = self.sys.st.clone_ref();
            self.sys.st = st.clone_ref();
            Ok(ret)
        } else {
            return Err(PyValueError::new_err("incompatible state"));
        }
    }

    fn replace_state_from(&mut self, arg: &Bound<'_, PyAny>) -> PyResult<StatePy> {
        if let Ok(sys) = arg.cast::<SystemPy>() {
            let st = sys.borrow().st.clone_ref();
            self.replace_state(&st)
        } else if let Ok(sel) = arg.cast::<SelPy>() {
            let st = sel.borrow().sys.st.clone_ref();
            self.replace_state(&st)
        } else {
            Err(PyTypeError::new_err(format!(
                "Invalid argument type {} in set_state_from()",
                arg.get_type()
            )))
        }
    }

    fn replace_system(&mut self, sys: &SystemPy) -> PyResult<SystemPy> {
        let ret = self.sys.clone_ref();
        self.sys = sys.clone_ref();
        Ok(ret)
    }

    #[getter("topology")]
    fn get_topology(&self) -> TopologyPy {
        self.sys.top.clone_ref()
    }

    #[setter("topology")]
    fn set_topology(&mut self, top: &TopologyPy) -> PyResult<()> {
        if self.sys.top.get().interchangeable(top.get()) {
            self.sys.top = top.clone_ref();
            Ok(())
        } else {
            return Err(PyValueError::new_err("incompatible topology"));
        }
    }

    fn replace_topology(&mut self, top: &TopologyPy) -> PyResult<TopologyPy> {
        if self.sys.top.get().interchangeable(top.get()) {
            let ret = self.sys.top.clone_ref();
            self.sys.top = top.clone_ref();
            Ok(ret)
        } else {
            return Err(PyValueError::new_err("incompatible topology"));
        }
    }

    pub fn set_same_chain(&mut self, val: char) {
        AtomIterMutProvider::set_same_chain(self, val)
    }

    pub fn set_same_resname(&mut self, val: &str) {
        AtomIterMutProvider::set_same_resname(self, val)
    }

    pub fn set_same_resid(&mut self, val: i32) {
        AtomIterMutProvider::set_same_resid(self, val)
    }

    pub fn set_same_name(&mut self, val: &str) {
        AtomIterMutProvider::set_same_name(self, val)
    }

    pub fn set_same_mass(&mut self, val: f32) {
        AtomIterMutProvider::set_same_mass(self, val)
    }

    pub fn set_same_bfactor(&mut self, val: f32) {
        AtomIterMutProvider::set_same_bfactor(self, val)
    }

    #[getter]
    fn get_time(&self) -> f32 {
        TimeProvider::get_time(self)
    }

    #[getter]
    fn get_box(&self) -> PeriodicBoxPy {
        PeriodicBoxPy(self.sys.st.require_box().unwrap().clone())
    }

    #[setter]
    fn set_box(&mut self, b: &PeriodicBoxPy) {
        self.sys.st.get_mut().pbox = Some(b.0.clone());
    }

    #[setter]
    fn set_time(&mut self, t: f32) {
        self.sys.st.get_mut().time = t;
    }

    fn set_box_from(&self, sys: Bound<'_, SystemPy>) {
        self.sys.st.get_mut().pbox = Some(sys.borrow().st.get().require_box().unwrap().clone());
    }

    #[pyo3(signature = (dims=[false,false,false]))]
    fn com<'py>(
        &self,
        py: Python<'py>,
        dims: [bool; 3],
    ) -> PyResult<Bound<'py, numpy::PyArray1<f32>>> {
        let pbc_dims = PbcDims::new(dims[0], dims[1], dims[2]);
        Ok(clone_vec_to_pyarray1(
            &MeasurePeriodic::center_of_mass_pbc_dims(self, pbc_dims)
                .map_err(to_py_runtime_err)?
                .coords,
            py,
        ))
    }

    #[pyo3(signature = (dims=[false,false,false]))]
    fn cog<'py>(
        &self,
        py: Python<'py>,
        dims: [bool; 3],
    ) -> PyResult<Bound<'py, numpy::PyArray1<f32>>> {
        let pbc_dims = PbcDims::new(dims[0], dims[1], dims[2]);
        Ok(clone_vec_to_pyarray1(
            &MeasurePeriodic::center_of_geometry_pbc_dims(self, pbc_dims)
                .map_err(to_py_runtime_err)?
                .coords,
            py,
        ))
    }

    fn principal_transform(&self) -> PyResult<crate::IsometryTransform> {
        let tr = MeasureMasses::principal_transform(self).map_err(to_py_runtime_err)?;
        Ok(crate::IsometryTransform(tr))
    }

    fn principal_transform_pbc(&self) -> PyResult<crate::IsometryTransform> {
        let tr = MeasurePeriodic::principal_transform_pbc(self).map_err(to_py_runtime_err)?;
        Ok(crate::IsometryTransform(tr))
    }

    fn apply_transform(&mut self, tr: &crate::IsometryTransform) {
        ModifyPos::apply_transform(self, &tr.0);
    }

    fn gyration(&self) -> PyResult<f32> {
        Ok(MeasureMasses::gyration(self).map_err(to_py_runtime_err)?)
    }

    fn gyration_pbc(&self) -> PyResult<f32> {
        Ok(MeasurePeriodic::gyration_pbc(self).map_err(to_py_runtime_err)?)
    }

    fn min_max<'py>(
        &self,
        py: Python<'py>,
    ) -> (Bound<'py, PyArray1<f32>>, Bound<'py, PyArray1<f32>>) {
        let (min, max) = MeasurePos::min_max(self);
        let minpy = clone_vec_to_pyarray1(&min.coords, py);
        let maxpy = clone_vec_to_pyarray1(&max.coords, py);
        (minpy, maxpy)
    }

    fn inertia<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<(
        Bound<'py, numpy::PyArray1<f32>>,
        Bound<'py, numpy::PyArray2<f32>>,
    )> {
        let (moments, axes) = MeasureMasses::inertia(self).map_err(to_py_runtime_err)?;
        let mom = clone_vec_to_pyarray1(&moments, py);
        let ax = axes.to_pyarray(py);
        Ok((mom, ax))
    }

    fn inertia_pbc<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<(
        Bound<'py, numpy::PyArray1<f32>>,
        Bound<'py, numpy::PyArray2<f32>>,
    )> {
        let (moments, axes) = MeasurePeriodic::inertia_pbc(self).map_err(to_py_runtime_err)?;
        let mom = clone_vec_to_pyarray1(&moments, py);
        let ax = axes.to_pyarray(py);
        Ok((mom, ax))
    }

    fn save(&self, fname: &str) -> PyResult<()> {
        Ok(SaveTopologyState::save(self, fname).map_err(to_py_runtime_err)?)
    }

    fn translate<'py>(&mut self, arg: PyArrayLike1<'py, f32>) -> PyResult<()> {
        let vec: VectorView<f32, Const<3>, Dyn> = arg
            .try_as_matrix()
            .ok_or_else(|| PyValueError::new_err("conversion to Vector3 has failed"))?;
        ModifyPos::translate(self, &vec);
        Ok(())
    }

    fn split_resindex(&self) -> Vec<SelPy> {
        AtomPosAnalysis::split_resindex(self)
            .map(|s| SelPy {
                sys: self.sys.clone_ref(),
                index: s.into_svec(),
            })
            .collect()
    }

    fn split_chain(&self) -> Vec<SelPy> {
        self.split(|p| Some(p.atom.chain))
            .map(|s| SelPy {
                sys: self.sys.clone_ref(),
                index: s.into_svec(),
            })
            .collect()
    }

    fn split_molecule(&self) -> Vec<SelPy> {
        self.split_mol_iter()
            .map(|s| SelPy {
                sys: self.sys.clone_ref(),
                index: s.into_svec(),
            })
            .collect()
    }

    fn to_gromacs_ndx(&self, name: &str) -> String {
        self.index.as_gromacs_ndx_str(name)
    }

    fn iter_pos(slf: PyRef<'_, Self>) -> PyRef<'_, SelPosIterator> {
        Bound::new(
            slf.py(),
            SelPosIterator {
                sel: slf.clone(),
                cur: 0,
            },
        )
        .unwrap()
        .borrow()
    }

    fn iter_atoms(slf: PyRef<'_, Self>) -> PyRef<'_, SelAtomIterator> {
        Bound::new(
            slf.py(),
            SelAtomIterator {
                sel: slf.clone(),
                cur: 0,
            },
        )
        .unwrap()
        .borrow()
    }
}

