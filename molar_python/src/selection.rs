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
/// Atom selection view with analysis and editing utilities.
///
/// Provides selection algebra, coordinate editing, and common analysis metrics.
///
/// **Example**
///
/// .. code-block:: python
///
///    sel = system("name CA")
///    com = sel.com()
///    rg = sel.gyration()

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
/// Iterator over selected atom positions.

#[pyclass(frozen)]
pub struct SelPosIterator {
    pub(crate) sel: Py<SelPy>,
    pub(crate) cur: AtomicUsize,
}

#[pymethods]
impl SelPosIterator {
    /// Return iterator object.
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    /// Return next selected position as NumPy array view.
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
/// Iterator over selected atoms.

#[pyclass(frozen)]
pub struct SelAtomIterator {
    pub(crate) sel: Py<SelPy>,
    pub(crate) cur: AtomicUsize,
}

#[pymethods]
impl SelAtomIterator {
    /// Return iterator object.
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    /// Return next selected atom view.
    fn __next__<'py>(slf: &Bound<'py, Self>) -> Option<AtomView> {
        let s = slf.get();
        let sel = s.sel.borrow(slf.py());
        if s.cur.load(Ordering::Relaxed) >= sel.len() {
            return None;
        }
        let idx = unsafe { sel.index.get_index_unchecked(s.cur.load(Ordering::Relaxed)) };
        s.cur.fetch_add(1, Ordering::Relaxed);
        Some(AtomView { top: sel.py_top().clone_ref(slf.py()), index: idx })
    }
}

/// Returns `true` if both selections share the same underlying topology Python object.
fn same_top(a: &SelPy, b: &SelPy) -> bool {
    a.py_top().as_ptr() == b.py_top().as_ptr()
}

#[pymethods]
impl SelPy {
    /// Number of atoms in this selection.
    ///
    /// :returns: Selection size.
    /// :rtype: int
    fn __len__(&self) -> usize {
        self.index.len()
    }

    fn __repr__(&self) -> String {
        format!("Sel(n={} atoms)", self.index.len())
    }

    /// Build sub-selection from query string, range, or explicit indices.
    ///
    /// :param arg: Selection expression, range tuple, or index list.
    /// :returns: Derived selection.
    /// :rtype: Sel
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

    /// Return one particle from the selection by index (supports negative indexing).
    ///
    /// :param i: Index within selection.
    /// :returns: Particle view.
    /// :rtype: Particle
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

    /// Iterate over particles in this selection.
    ///
    /// :returns: Particle iterator.
    /// :rtype: Iterator[Particle]
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

    /// Global atom indices of this selection.
    ///
    /// :returns: Index array.
    /// :rtype: numpy.ndarray
    #[getter("index")]
    fn get_index<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray1<usize>> {
        numpy::PyArray1::from_iter(py, self.index.iter_index())
    }

    /// Coordinates of selected atoms as an array of shape ``[3, n_atoms]``.
    ///
    /// :returns: Coordinate array.
    /// :rtype: numpy.ndarray
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

    /// A ``System`` sharing this selection's topology and state.
    ///
    /// :returns: System view.
    /// :rtype: System
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

    /// Replace backing topology and state from another system.
    ///
    /// :param sys: Source system.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter("system")]
    fn set_system(&self, sys: &Bound<SystemPy>) -> PyResult<()> {
        let py = sys.py();
        self.set_topology(sys.get().py_top().bind(py))?;
        self.set_state(sys.get().py_st().bind(py))?;
        Ok(())
    }

    /// Set coordinates from an array of shape ``[3, n_atoms]``.
    ///
    /// :param arr: Coordinate array.
    /// :returns: ``None``.
    /// :rtype: None
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

    /// Backing state object.
    ///
    /// :returns: Backing state.
    /// :rtype: State
    #[getter("state")]
    fn get_state(slf: Bound<Self>) -> Bound<StatePy> {
        slf.get().py_st().bind(slf.py()).clone()
    }

    /// Replace backing state with a compatible state.
    ///
    /// :param st: New state.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter("state")]
    fn set_state(&self, st: &Bound<StatePy>) -> PyResult<()> {
        if self.r_st().interchangeable(st.get().inner()) {
            *self.py_st_mut() = st.clone().unbind();
            Ok(())
        } else {
            return Err(PyValueError::new_err("incompatible state"));
        }
    }

    /// Backing topology object.
    ///
    /// :returns: Backing topology.
    /// :rtype: Topology
    #[getter("topology")]
    fn get_topology(slf: Bound<Self>) -> Bound<TopologyPy> {
        slf.get().py_top().bind(slf.py()).clone()
    }

    /// Replace backing topology with a compatible topology.
    ///
    /// :param top: New topology.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter("topology")]
    fn set_topology(&self, top: &Bound<TopologyPy>) -> PyResult<()> {
        if self.r_top().interchangeable(top.get().inner()) {
            *self.py_top_mut() = top.clone().unbind();
            Ok(())
        } else {
            return Err(PyValueError::new_err("incompatible topology"));
        }
    }

    /// Replace state data in-place by swapping with a compatible state.
    ///
    /// :param st: Compatible state object.
    /// :returns: ``None``.
    /// :rtype: None
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

    /// Set chain ID for all selected atoms in-place.
    ///
    /// :param val: New chain character.
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    sel.set_same_chain('A')
    pub fn set_same_chain(&self, val: char) {
        AtomIterMutProvider::set_same_chain(self.r_top_mut(), val)
    }

    /// Set residue name for all selected atoms in-place.
    ///
    /// :param val: New residue name.
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    sel.set_same_resname('ALA')
    pub fn set_same_resname(&self, val: &str) {
        AtomIterMutProvider::set_same_resname(self.r_top_mut(), val)
    }

    /// Set residue ID for all selected atoms in-place.
    ///
    /// :param val: New residue ID.
    pub fn set_same_resid(&self, val: i32) {
        AtomIterMutProvider::set_same_resid(self.r_top_mut(), val)
    }

    /// Set atom name for all selected atoms in-place.
    ///
    /// :param val: New atom name.
    pub fn set_same_name(&self, val: &str) {
        AtomIterMutProvider::set_same_name(self.r_top_mut(), val)
    }

    /// Set atomic mass for all selected atoms in-place.
    ///
    /// :param val: New mass in Da.
    pub fn set_same_mass(&self, val: f32) {
        AtomIterMutProvider::set_same_mass(self.r_top_mut(), val)
    }

    /// Set B-factor for all selected atoms in-place.
    ///
    /// :param val: New B-factor value.
    pub fn set_same_bfactor(&self, val: f32) {
        AtomIterMutProvider::set_same_bfactor(self.r_top_mut(), val)
    }

    /// Selection time value (proxied to backing state).
    ///
    /// :returns: Time value.
    /// :rtype: float
    #[getter]
    fn get_time(&self) -> f32 {
        TimeProvider::get_time(self)
    }

    /// Set selection time value (proxied to backing state).
    ///
    /// :param t: New time value.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter]
    fn set_time(&self, t: f32) {
        self.r_st_mut().time = t;
    }

    /// Periodic box (proxied to backing state).
    ///
    /// :returns: Periodic box.
    /// :rtype: PeriodicBox
    #[getter]
    fn get_box(&self) -> PeriodicBoxPy {
        PeriodicBoxPy(self.r_st().require_box().unwrap().clone())
    }

    /// Set periodic box (proxied to backing state).
    ///
    /// :param b: New periodic box.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter]
    fn set_box(&self, b: &PeriodicBoxPy) {
        self.r_st_mut().pbox = Some(b.0.clone());
    }

    /// Copy periodic box from a system.
    ///
    /// :param sys: Source system.
    /// :returns: ``None``.
    /// :rtype: None
    fn set_box_from(&self, sys: Bound<'_, SystemPy>) {
        self.r_st_mut().pbox = Some(sys.get().r_st().require_box().unwrap().clone());
    }

    // Analysis functions

    #[pyo3(signature = (dims=None))]
    #[pyo3(text_signature = "($self, dims=None)")]
    /// Center of mass, optionally using periodic dimensions.
    ///
    /// :param dims: Periodic dimensions ``[x, y, z]`` booleans.
    /// :returns: Center-of-mass vector ``[x, y, z]`` in nm.
    /// :rtype: numpy.ndarray
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    center = sel.com()                           # [x, y, z] in nm
    ///    center = sel.com(dims=[True, True, True])    # with PBC
    fn com<'py>(
        &self,
        py: Python<'py>,
        dims: Option<[bool; 3]>,
    ) -> PyResult<Bound<'py, numpy::PyArray1<f32>>> {
        let dims = dims.unwrap_or([false, false, false]);
        let pbc_dims = PbcDims::new(dims[0], dims[1], dims[2]);
        Ok(clone_vec_to_pyarray1(
            &MeasurePeriodic::center_of_mass_pbc_dims(self, pbc_dims)
                .map_err(to_py_runtime_err)?
                .coords,
            py,
        ))
    }

    #[pyo3(signature = (dims=None))]
    #[pyo3(text_signature = "($self, dims=None)")]
    /// Center of geometry, optionally using periodic dimensions.
    ///
    /// :param dims: Periodic dimensions ``[x, y, z]`` booleans.
    /// :returns: Center-of-geometry vector ``[x, y, z]`` in nm.
    /// :rtype: numpy.ndarray
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    center = sel.cog()                           # [x, y, z] in nm
    ///    center = sel.cog(dims=[True, True, True])    # with PBC
    fn cog<'py>(
        &self,
        py: Python<'py>,
        dims: Option<[bool; 3]>,
    ) -> PyResult<Bound<'py, numpy::PyArray1<f32>>> {
        let dims = dims.unwrap_or([false, false, false]);
        let pbc_dims = PbcDims::new(dims[0], dims[1], dims[2]);
        Ok(clone_vec_to_pyarray1(
            &MeasurePeriodic::center_of_geometry_pbc_dims(self, pbc_dims)
                .map_err(to_py_runtime_err)?
                .coords,
            py,
        ))
    }

    /// Principal-axes alignment transform (non-periodic).
    ///
    /// :returns: Rigid transform.
    /// :rtype: IsometryTransform
    fn principal_transform(&self) -> PyResult<crate::IsometryTransform> {
        let tr = MeasureMasses::principal_transform(self).map_err(to_py_runtime_err)?;
        Ok(crate::IsometryTransform(tr))
    }

    /// Principal-axes alignment transform with periodic handling.
    ///
    /// :returns: Rigid transform.
    /// :rtype: IsometryTransform
    fn principal_transform_pbc(&self) -> PyResult<crate::IsometryTransform> {
        let tr = MeasurePeriodic::principal_transform_pbc(self).map_err(to_py_runtime_err)?;
        Ok(crate::IsometryTransform(tr))
    }

    /// Apply rigid transform to selected coordinates.
    ///
    /// :param tr: Transform to apply.
    /// :returns: ``None``.
    /// :rtype: None
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    import pymolar
    ///    tr = pymolar.fit_transform(mobile, ref)
    ///    mobile.apply_transform(tr)
    fn apply_transform(&self, tr: &crate::IsometryTransform) {
        TmpSel {
            top: self.r_top(),
            st: self.r_st(),
            index: &self.index,
        }
        .apply_transform(&tr.0);
    }

    /// Radius of gyration (non-periodic).
    ///
    /// :returns: Radius of gyration in nm.
    /// :rtype: float
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    rg = sel.gyration()    # radius of gyration in nm
    fn gyration(&self) -> PyResult<f32> {
        Ok(MeasureMasses::gyration(self).map_err(to_py_runtime_err)?)
    }

    /// Radius of gyration with periodic handling.
    ///
    /// :returns: Radius of gyration.
    /// :rtype: float
    fn gyration_pbc(&self) -> PyResult<f32> {
        Ok(MeasurePeriodic::gyration_pbc(self).map_err(to_py_runtime_err)?)
    }

    /// Compute DSSP secondary structure assignment for each residue.
    ///
    /// Implements the Kabsch & Sander (1983) algorithm.  Pass a **protein-only**
    /// selection; non-protein residues lack backbone atoms and will appear as
    /// ``'='`` (break) in the output.
    ///
    /// :returns: List of single-character DSSP codes, one per residue.
    ///
    ///     ========= =================================
    ///     Character Meaning
    ///     ========= =================================
    ///     ``H``     Alpha helix
    ///     ``G``     3\ :sub:`10` helix
    ///     ``I``     π helix
    ///     ``P``     Poly-proline II helix
    ///     ``E``     Extended beta strand (sheet)
    ///     ``B``     Isolated beta bridge
    ///     ``T``     Hydrogen-bonded turn
    ///     ``S``     Bend (Cα angle ≥ 70°)
    ///     ``~``     Loop / coil
    ///     ``=``     Break (missing backbone atoms)
    ///     ========= =================================
    ///
    /// :rtype: list[str]
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    prot = sys("protein")
    ///    codes = prot.dssp()            # ['H', 'H', 'T', 'E', '~', ...]
    ///    helix_count = codes.count('H')
    fn dssp(&self) -> Vec<String> {
        MeasureAtomPos::dssp(self)
            .ss()
            .iter()
            .map(|ss| ss.to_char().to_string())
            .collect()
    }

    /// Compact DSSP secondary structure string, one character per residue.
    ///
    /// Equivalent to ``''.join(sel.dssp())``.  See :meth:`dssp` for the code
    /// table and caveats about non-protein residues.
    ///
    /// :returns: DSSP string, e.g. ``"HHHHHTTEEEE~~~~~"``.
    /// :rtype: str
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    prot = sys("protein")
    ///    print(prot.dssp_string())    # e.g. "HHHHHTTEEEE~~~~~"
    fn dssp_string(&self) -> String {
        MeasureAtomPos::dssp(self).ss_string()
    }

    /// Axis-aligned bounding-box min and max coordinates.
    ///
    /// :returns: Tuple ``(min_xyz, max_xyz)`` as NumPy arrays.
    /// :rtype: tuple[numpy.ndarray, numpy.ndarray]
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    lo, hi = sel.min_max()    # two [x, y, z] arrays in nm
    fn min_max<'py>(
        &self,
        py: Python<'py>,
    ) -> (Bound<'py, PyArray1<f32>>, Bound<'py, PyArray1<f32>>) {
        let (min, max) = MeasurePos::min_max(self);
        let minpy = clone_vec_to_pyarray1(&min.coords, py);
        let maxpy = clone_vec_to_pyarray1(&max.coords, py);
        (minpy, maxpy)
    }

    /// Principal moments and axes of inertia tensor.
    ///
    /// :returns: Tuple ``(moments, axes)`` where moments is length-3 and axes is 3×3.
    /// :rtype: tuple[numpy.ndarray, numpy.ndarray]
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    moments, axes = sel.inertia()
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

    /// Principal moments and axes of inertia tensor with periodic handling.
    ///
    /// :returns: Tuple ``(moments, axes)``.
    /// :rtype: tuple[numpy.ndarray, numpy.ndarray]
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

    /// Save topology and selected state coordinates to file.
    ///
    /// :param fname: Output file path.
    /// :returns: ``None``.
    /// :rtype: None
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    sel.save("subset.pdb")
    fn save(&self, fname: &str) -> PyResult<()> {
        Ok(SaveTopologyState::save(self, fname).map_err(to_py_runtime_err)?)
    }

    /// Translate selected coordinates by vector.
    ///
    /// :param arg: Translation vector ``[x, y, z]`` in nm.
    /// :returns: ``None``.
    /// :rtype: None
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    import numpy as np
    ///    sel.translate(np.array([0.5, 0.0, 0.0], dtype=np.float32))  # shift 0.5 nm along x
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

    /// Split selection into sub-selections by residue index.
    ///
    /// :returns: List of per-residue selections.
    /// :rtype: list[Sel]
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    residues = sel.split_resindex()    # list of Sel, one per residue
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

    /// Split selection into sub-selections by chain ID.
    ///
    /// :returns: List of per-chain selections.
    /// :rtype: list[Sel]
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    chains = sel.split_chain()    # list of Sel, one per chain
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

    /// Split selection into molecular connected components.
    ///
    /// :returns: List of per-molecule selections.
    /// :rtype: list[Sel]
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    mols = sel.split_molecule()    # list of Sel, one per molecule
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

    /// Export selection indices in GROMACS NDX group format.
    ///
    /// :param name: NDX group name.
    /// :returns: NDX formatted group block.
    /// :rtype: str
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    ndx_str = sel.to_gromacs_ndx("Protein")
    fn to_gromacs_ndx(&self, name: &str) -> String {
        self.index.as_gromacs_ndx_str(name)
    }

    /// Iterate over selected positions.
    ///
    /// :returns: Position iterator.
    /// :rtype: SelPosIterator
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

    /// Iterate over selected atoms.
    ///
    /// :returns: Atom iterator.
    /// :rtype: SelAtomIterator
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

    /// Test whether a global atom index is in this selection.
    ///
    /// :param idx: Global atom index.
    /// :returns: ``True`` if the index belongs to this selection.
    /// :rtype: bool
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    prot = sys("protein")
    ///    print(42 in prot)    # True or False
    fn __contains__(&self, idx: usize) -> bool {
        self.index.contains(&idx)
    }

    /// Union of two selections (``sel1 | sel2``); both must belong to the same system.
    ///
    /// :param other: Second selection.
    /// :returns: Union selection.
    /// :rtype: Sel
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    prot  = sys("protein")
    ///    bb    = sys("backbone")
    ///    union = prot | bb    # all protein + backbone atoms
    fn __or__(slf: &Bound<'_, Self>, other: &SelPy) -> PyResult<SelPy> {
        let s = slf.get();
        if !same_top(s, other) {
            return Err(PyValueError::new_err("selections must belong to the same system"));
        }
        let a = s.index();
        let b = other.index();
        let mut result: Vec<usize> = Vec::with_capacity(a.len() + b.len());
        let (mut i, mut j) = (0, 0);
        while i < a.len() && j < b.len() {
            if a[i] < b[j]      { result.push(a[i]); i += 1; }
            else if a[i] > b[j] { result.push(b[j]); j += 1; }
            else                { result.push(a[i]); i += 1; j += 1; }
        }
        result.extend_from_slice(&a[i..]);
        result.extend_from_slice(&b[j..]);
        let index: SVec = result.into_iter().collect();
        Ok(SelPy::new(s.py_top().clone_ref(slf.py()), s.py_st().clone_ref(slf.py()), index))
    }

    /// Intersection of two selections (``sel1 & sel2``); raises ``ValueError`` if empty.
    ///
    /// :param other: Second selection.
    /// :returns: Intersection selection.
    /// :rtype: Sel
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    prot  = sys("protein")
    ///    bb    = sys("backbone")
    ///    inter = prot & bb    # backbone atoms only
    fn __and__(slf: &Bound<'_, Self>, other: &SelPy) -> PyResult<SelPy> {
        let s = slf.get();
        if !same_top(s, other) {
            return Err(PyValueError::new_err("selections must belong to the same system"));
        }
        let a = s.index();
        let b = other.index();
        let mut result: Vec<usize> = Vec::new();
        let (mut i, mut j) = (0, 0);
        while i < a.len() && j < b.len() {
            if a[i] < b[j]      { i += 1; }
            else if a[i] > b[j] { j += 1; }
            else                { result.push(a[i]); i += 1; j += 1; }
        }
        if result.is_empty() {
            return Err(PyValueError::new_err("empty intersection"));
        }
        let index: SVec = result.into_iter().collect();
        Ok(SelPy::new(s.py_top().clone_ref(slf.py()), s.py_st().clone_ref(slf.py()), index))
    }

    /// Set difference (``sel1 - sel2``, atoms in sel1 not in sel2); raises ``ValueError`` if empty.
    ///
    /// :param other: Second selection.
    /// :returns: Difference selection.
    /// :rtype: Sel
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    prot  = sys("protein")
    ///    bb    = sys("backbone")
    ///    diff  = prot - bb    # sidechain atoms
    fn __sub__(slf: &Bound<'_, Self>, other: &SelPy) -> PyResult<SelPy> {
        let s = slf.get();
        if !same_top(s, other) {
            return Err(PyValueError::new_err("selections must belong to the same system"));
        }
        let a = s.index();
        let b = other.index();
        let mut result: Vec<usize> = Vec::new();
        let (mut i, mut j) = (0, 0);
        while i < a.len() {
            if j >= b.len() || a[i] < b[j] { result.push(a[i]); i += 1; }
            else if a[i] > b[j]            { j += 1; }
            else                           { i += 1; j += 1; }
        }
        if result.is_empty() {
            return Err(PyValueError::new_err("empty difference"));
        }
        let index: SVec = result.into_iter().collect();
        Ok(SelPy::new(s.py_top().clone_ref(slf.py()), s.py_st().clone_ref(slf.py()), index))
    }

    /// Complement: all system atoms NOT in this selection (``~sel``); raises ``ValueError`` if empty.
    ///
    /// :returns: Complement selection.
    /// :rtype: Sel
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    bb     = sys("backbone")
    ///    compl  = ~bb    # non-backbone atoms
    fn __invert__(slf: &Bound<'_, Self>) -> PyResult<SelPy> {
        let s = slf.get();
        let n = s.r_top().len();
        let sel_idx = s.index();
        let mut result: Vec<usize> = Vec::with_capacity(n.saturating_sub(sel_idx.len()));
        let mut j = 0usize;
        for i in 0..n {
            if j < sel_idx.len() && sel_idx[j] == i { j += 1; }
            else { result.push(i); }
        }
        if result.is_empty() {
            return Err(PyValueError::new_err(
                "empty complement (selection covers all atoms)",
            ));
        }
        let index: SVec = result.into_iter().collect();
        Ok(SelPy::new(s.py_top().clone_ref(slf.py()), s.py_st().clone_ref(slf.py()), index))
    }

    /// Expand selection to include all atoms in the same residues as any selected atom.
    ///
    /// :returns: Expanded selection.
    /// :rtype: Sel
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    ca  = sys("name CA")
    ///    res = ca.whole_residues()    # all atoms in CA-containing residues
    fn whole_residues(slf: &Bound<'_, Self>) -> SelPy {
        let s = slf.get();
        let sel = AtomPosAnalysis::whole_residues(s);
        let index: SVec = sel.iter_index().collect();
        SelPy::new(s.py_top().clone_ref(slf.py()), s.py_st().clone_ref(slf.py()), index)
    }

    /// Expand selection to include all atoms in the same chains as any selected atom.
    ///
    /// :returns: Expanded selection.
    /// :rtype: Sel
    ///
    /// **Example**
    ///
    /// .. code-block:: python
    ///
    ///    ca    = sys("name CA")
    ///    chain = ca.whole_chains()    # all atoms in those chains
    fn whole_chains(slf: &Bound<'_, Self>) -> SelPy {
        let s = slf.get();
        let sel = AtomPosAnalysis::whole_chains(s);
        let index: SVec = sel.iter_index().collect();
        SelPy::new(s.py_top().clone_ref(slf.py()), s.py_st().clone_ref(slf.py()), index)
    }
}
