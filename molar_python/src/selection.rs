use std::cell::UnsafeCell;
use std::sync::atomic::{AtomicUsize, Ordering};

use molar::prelude::*;
use numpy::nalgebra::{Const, Dyn, VectorView};
use numpy::{
    PyArray1, PyArrayLike1, PyArrayMethods, PyReadonlyArray2, PyUntypedArrayMethods, ToPyArray,
};
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::{exceptions::PyIndexError, prelude::*, types::PyAny};

use crate::atom::AtomView;
use crate::periodic_box::PeriodicBoxPy;
use crate::system::SystemPy;
use crate::topology_state::{StatePy, TopologyPy};
use crate::utils::*;
use crate::{ParticleIterator, ParticlePy};

#[pyclass(name = "Sel", frozen)]
pub struct SelPy {
    top: UnsafeCell<Py<TopologyPy>>,
    st: UnsafeCell<Py<StatePy>>,
    index: SVec,
}

unsafe impl Send for SelPy {}
unsafe impl Sync for SelPy {}

impl SelPy {
    pub(crate) fn new(py_top: Py<TopologyPy>, py_st: Py<StatePy>, index: SVec) -> Self {
        Self {
            top: UnsafeCell::new(py_top),
            st: UnsafeCell::new(py_st),
            index,
        }
    }

    pub(crate) fn index(&self) -> &[usize] {
        &self.index
    }

    pub(crate) fn r_top(&self) -> &Topology {
        unsafe { &*self.top.get() }.get().inner()
    }

    pub(crate) fn r_top_mut(&self) -> &mut Topology {
        unsafe { &*self.top.get() }.get().inner_mut()
    }

    pub(crate) fn r_st(&self) -> &State {
        unsafe { &*self.st.get() }.get().inner()
    }

    pub(crate) fn r_st_mut(&self) -> &mut State {
        unsafe { &*self.st.get() }.get().inner_mut()
    }

    pub(crate) fn py_top(&self) -> &Py<TopologyPy> {
        unsafe { &*self.top.get() }
    }

    pub(crate) fn py_top_mut(&self) -> &mut Py<TopologyPy> {
        unsafe { &mut *self.top.get() }
    }

    pub(crate) fn py_st(&self) -> &Py<StatePy> {
        unsafe { &*self.st.get() }
    }

    pub(crate) fn py_st_mut(&self) -> &mut Py<StatePy> {
        unsafe { &mut *self.st.get() }
    }
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
        self.r_top().atoms.as_ptr()
    }

    fn coords_ptr(&self) -> *const Pos {
        self.r_st().coords.as_ptr()
    }
}

impl AtomPosAnalysisMut for SelPy {}

impl BoxProvider for SelPy {
    fn get_box(&self) -> Option<&PeriodicBox> {
        let ptr = self.r_st().pbox.as_ref().map(|b| b as *const PeriodicBox)?;
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
        self.r_st().time
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
        self.r_top().get_molecule_unchecked(i)
    }

    fn num_molecules(&self) -> usize {
        self.r_top().num_molecules()
    }
}

impl SelPy {
    pub fn from_svec(&self, index: SVec) -> Self {
        Python::attach(|py| Self {
            top: UnsafeCell::new(self.py_top().clone_ref(py)),
            st: UnsafeCell::new(self.py_st().clone_ref(py)),
            index,
        })
    }
}

impl Clone for SelPy {
    fn clone(&self) -> Self {
        Python::attach(|py| SelPy {
            top: UnsafeCell::new(self.py_top().clone_ref(py)),
            st: UnsafeCell::new(self.py_st().clone_ref(py)),
            index: self.index.clone(),
        })
    }
}

//-------------------------------------------
pub struct TmpSel<'a> {
    pub(crate) top: &'a Topology,
    pub(crate) st: &'a State,
    pub(crate) index: &'a [usize],
}

impl IndexSliceProvider for TmpSel<'_> {
    fn get_index_slice(&self) -> &[usize] {
        self.index
    }
}

impl AtomPosAnalysis for TmpSel<'_> {
    fn atoms_ptr(&self) -> *const Atom {
        self.top.atoms.as_ptr()
    }

    fn coords_ptr(&self) -> *const Pos {
        self.st.coords.as_ptr()
    }
}

impl AtomPosAnalysisMut for TmpSel<'_> {}

//-----------------------------------------

#[pyclass(frozen)]
pub struct SelPosIterator {
    pub(crate) sel: Py<SelPy>,
    pub(crate) cur: AtomicUsize,
}

#[pymethods]
impl SelPosIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__<'py>(slf: &Bound<'py, Self>) -> Option<Bound<'py, PyArray1<f32>>> {
        let s = slf.get();
        let sel = s.sel.get();
        if s.cur.load(Ordering::Relaxed) >= sel.len() {
            return None;
        }
        let idx = unsafe { sel.index.get_index_unchecked(s.cur.load(Ordering::Relaxed)) };
        s.cur.fetch_add(1, Ordering::Relaxed);
        unsafe { Some(map_pyarray_to_pos(sel.py_st().bind(slf.py()), idx)) }
    }
}

#[pyclass(frozen)]
pub struct SelAtomIterator {
    pub(crate) sel: Py<SelPy>,
    pub(crate) cur: AtomicUsize,
}

#[pymethods]
impl SelAtomIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__<'py>(slf: &Bound<'py, Self>) -> Option<AtomView> {
        let s = slf.get();
        let sel = s.sel.borrow(slf.py());
        if s.cur.load(Ordering::Relaxed) >= sel.len() {
            return None;
        }
        let idx = unsafe { sel.index.get_index_unchecked(s.cur.load(Ordering::Relaxed)) };
        s.cur.fetch_add(1, Ordering::Relaxed);
        let atom_ptr = unsafe { sel.atoms_ptr().add(idx) as *mut Atom };
        Some(AtomView(atom_ptr))
    }
}

#[pymethods]
impl SelPy {
    fn __len__(&self) -> usize {
        self.index.len()
    }

    fn __call__(&self, arg: &Bound<'_, PyAny>) -> PyResult<SelPy> {
        // Construct temp system
        let sys = SystemPy::new(
            self.py_top().clone_ref(arg.py()),
            self.py_st().clone_ref(arg.py()),
        );
        let v = if let Ok(val) = arg.extract::<String>() {
            val.into_sel_index(&sys, Some(self.index.as_slice()))
                .map_err(to_py_runtime_err)?
        } else if let Ok(val) = arg.extract::<(usize, usize)>() {
            (val.0..=val.1)
                .into_sel_index(&sys, Some(self.index.as_slice()))
                .map_err(to_py_runtime_err)?
        } else if let Ok(val) = arg.extract::<Vec<usize>>() {
            val.into_sel_index(&sys, Some(self.index.as_slice()))
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

        Python::attach(|py| {
            Ok(ParticlePy {
                top: self.py_top().clone_ref(py),
                st: self.py_st().clone_ref(py),
                id: ind,
            })
        })
    }

    fn __iter__(slf: Bound<'_, Self>) -> Bound<'_, ParticleIterator> {
        Bound::new(
            slf.py(),
            ParticleIterator {
                sel: slf.clone().unbind(),
                cur: AtomicUsize::new(0),
            },
        )
        .unwrap()
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
    fn get_system<'py>(slf: Bound<'py, Self>) -> Bound<'py, SystemPy> {
        Bound::new(
            slf.py(),
            SystemPy::new(
                slf.get().py_top().clone_ref(slf.py()),
                slf.get().py_st().clone_ref(slf.py()),
            ),
        )
        .unwrap()
    }

    #[setter("system")]
    fn set_system(&self, sys: &Bound<SystemPy>) -> PyResult<()> {
        let py = sys.py();
        self.set_topology(sys.get().py_top().bind(py))?;
        self.set_state(sys.get().py_st().bind(py))?;
        Ok(())
    }

    #[setter("coords")]
    fn set_coords(&self, arr: PyReadonlyArray2<f32>) -> PyResult<()> {
        if arr.shape() != [3, self.__len__()] {
            return Err(PyValueError::new_err(format!(
                "Array shape must be [3, {}], not {:?}",
                self.__len__(),
                arr.shape()
            )));
        }
        let arr_ptr = arr.data();
        let coord_ptr = self.coords_ptr() as *const f32 as *mut f32;

        unsafe {
            for i in self.index.iter_index() {
                let pos_ptr = coord_ptr.add(i * 3);
                std::ptr::copy_nonoverlapping(arr_ptr.add(i * 3), pos_ptr, 3);
            }
        }
        Ok(())
    }

    #[getter("state")]
    fn get_state(slf: Bound<Self>) -> Bound<StatePy> {
        slf.get().py_st().bind(slf.py()).clone()
    }

    #[setter("state")]
    fn set_state(&self, st: &Bound<StatePy>) -> PyResult<()> {
        if self.r_st().interchangeable(st.get().inner()) {
            *self.py_st_mut() = st.clone().unbind();
            Ok(())
        } else {
            return Err(PyValueError::new_err("incompatible state"));
        }
    }

    #[getter("topology")]
    fn get_topology(slf: Bound<Self>) -> Bound<TopologyPy> {
        slf.get().py_top().bind(slf.py()).clone()
    }

    #[setter("topology")]
    fn set_topology(&self, top: &Bound<TopologyPy>) -> PyResult<()> {
        if self.r_top().interchangeable(top.get().inner()) {
            *self.py_top_mut() = top.clone().unbind();
            Ok(())
        } else {
            return Err(PyValueError::new_err("incompatible topology"));
        }
    }

    fn replace_state_deep(&self, st: &Bound<StatePy>) -> PyResult<()> {
        if self.r_st().interchangeable(st.get().inner()) {
            unsafe { std::ptr::swap(self.r_st_mut(), st.get().inner_mut()) };
            Ok(())
        } else {
            return Err(PyValueError::new_err("incompatible state"));
        }
    }
    // fn replace_state_from(&self, arg: &Bound<'_, PyAny>) -> PyResult<StatePy> {
    //     if let Ok(sys) = arg.cast::<SystemPy>() {
    //         let st = sys.borrow().st.clone_ref();
    //         self.replace_state(&st)
    //     } else if let Ok(sel) = arg.cast::<SelPy>() {
    //         let st = sel.borrow().sys.st.clone_ref();
    //         self.replace_state(&st)
    //     } else {
    //         Err(PyTypeError::new_err(format!(
    //             "Invalid argument type {} in set_state_from()",
    //             arg.get_type()
    //         )))
    //     }
    // }

    // fn replace_system(&self, sys: &SystemPy) -> PyResult<SystemPy> {
    //     let ret = self.sys.clone_ref();
    //     self.sys = sys.clone_ref();
    //     Ok(ret)
    // }

    // fn replace_topology(&self, top: &TopologyPy) -> PyResult<TopologyPy> {
    //     if self.sys.top.inner().interchangeable(top.inner()) {
    //         let ret = self.sys.top.clone_ref();
    //         self.sys.top = top.clone_ref();
    //         Ok(ret)
    //     } else {
    //         return Err(PyValueError::new_err("incompatible topology"));
    //     }
    // }

    // fn replace_topology_deep(&self, top: &mut TopologyPy) -> PyResult<()> {
    //     if self.sys.top.inner().interchangeable(top.inner()) {
    //         mem::swap(self.sys.top.inner_mut(), top.inner_mut());
    //         Ok(())
    //     } else {
    //         return Err(PyValueError::new_err("incompatible topology"));
    //     }
    // }

    pub fn set_same_chain(&self, val: char) {
        AtomIterMutProvider::set_same_chain(self.r_top_mut(), val)
    }

    pub fn set_same_resname(&self, val: &str) {
        AtomIterMutProvider::set_same_resname(self.r_top_mut(), val)
    }

    pub fn set_same_resid(&self, val: i32) {
        AtomIterMutProvider::set_same_resid(self.r_top_mut(), val)
    }

    pub fn set_same_name(&self, val: &str) {
        AtomIterMutProvider::set_same_name(self.r_top_mut(), val)
    }

    pub fn set_same_mass(&self, val: f32) {
        AtomIterMutProvider::set_same_mass(self.r_top_mut(), val)
    }

    pub fn set_same_bfactor(&self, val: f32) {
        AtomIterMutProvider::set_same_bfactor(self.r_top_mut(), val)
    }

    #[getter]
    fn get_time(&self) -> f32 {
        TimeProvider::get_time(self)
    }

    #[setter]
    fn set_time(&self, t: f32) {
        self.r_st_mut().time = t;
    }

    #[getter]
    fn get_box(&self) -> PeriodicBoxPy {
        PeriodicBoxPy(self.r_st().require_box().unwrap().clone())
    }

    #[setter]
    fn set_box(&self, b: &PeriodicBoxPy) {
        self.r_st_mut().pbox = Some(b.0.clone());
    }

    fn set_box_from(&self, sys: Bound<'_, SystemPy>) {
        self.r_st_mut().pbox = Some(sys.get().r_st().require_box().unwrap().clone());
    }

    // Analysis functions

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

    fn apply_transform(&self, tr: &crate::IsometryTransform) {
        TmpSel {
            top: self.r_top(),
            st: self.r_st(),
            index: &self.index,
        }
        .apply_transform(&tr.0);
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

    fn translate<'py>(&self, arg: PyArrayLike1<'py, f32>) -> PyResult<()> {
        let vec: VectorView<f32, Const<3>, Dyn> = arg
            .try_as_matrix()
            .ok_or_else(|| PyValueError::new_err("conversion to Vector3 has failed"))?;
        TmpSel {
            top: self.r_top(),
            st: self.r_st(),
            index: &self.index,
        }
        .translate(&vec);
        Ok(())
    }

    fn split_resindex(&self) -> Vec<SelPy> {
        Python::attach(|py| {
            AtomPosAnalysis::split_resindex(self)
                .map(|s| SelPy {
                    top: UnsafeCell::new(self.py_top().clone_ref(py)),
                    st: UnsafeCell::new(self.py_st().clone_ref(py)),
                    index: s.into_svec(),
                })
                .collect()
        })
    }

    fn split_chain(&self) -> Vec<SelPy> {
        Python::attach(|py| {
            self.split(|p| Some(p.atom.chain))
                .map(|s| SelPy {
                    top: UnsafeCell::new(self.py_top().clone_ref(py)),
                    st: UnsafeCell::new(self.py_st().clone_ref(py)),
                    index: s.into_svec(),
                })
                .collect()
        })
    }

    fn split_molecule(&self) -> Vec<SelPy> {
        Python::attach(|py| {
            self.split_mol_iter()
                .map(|s| SelPy {
                    top: UnsafeCell::new(self.py_top().clone_ref(py)),
                    st: UnsafeCell::new(self.py_st().clone_ref(py)),
                    index: s.into_svec(),
                })
                .collect()
        })
    }

    fn to_gromacs_ndx(&self, name: &str) -> String {
        self.index.as_gromacs_ndx_str(name)
    }

    fn iter_pos(slf: Bound<'_, Self>) -> Bound<'_, SelPosIterator> {
        Bound::new(
            slf.py(),
            SelPosIterator {
                sel: slf.clone().unbind(),
                cur: AtomicUsize::new(0),
            },
        )
        .unwrap()
    }

    fn iter_atoms(slf: Bound<'_, Self>) -> Bound<'_, SelAtomIterator> {
        Bound::new(
            slf.py(),
            SelAtomIterator {
                sel: slf.clone().unbind(),
                cur: AtomicUsize::new(0),
            },
        )
        .unwrap()
    }
}
