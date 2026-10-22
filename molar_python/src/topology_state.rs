use std::cell::UnsafeCell;
use molar_ff::{ApplyFF, ApplyCharges};
use crate::force_field::{ChargeArg, FFArg, to_charge_error, to_ff_error};

use crate::{SelPy, SystemPy};

use super::periodic_box::PeriodicBoxPy;
use molar::prelude::*;
use numpy::{AllowTypeChange, PyArray2, PyArrayLike2, PyArrayMethods, PyUntypedArrayMethods};
use pyo3::{
    exceptions::{PyAttributeError, PyTypeError, PyValueError},
    prelude::*,
};

/// Python names of [`BondOrder`] variants.
fn order_name(o: BondOrder) -> &'static str {
    match o {
        BondOrder::Unspecified => "unspecified",
        BondOrder::Single => "single",
        BondOrder::Double => "double",
        BondOrder::Triple => "triple",
        BondOrder::Aromatic => "aromatic",
    }
}

fn order_from_name(name: &str) -> PyResult<BondOrder> {
    Ok(match name {
        "unspecified" => BondOrder::Unspecified,
        "single" => BondOrder::Single,
        "double" => BondOrder::Double,
        "triple" => BondOrder::Triple,
        "aromatic" => BondOrder::Aromatic,
        _ => {
            return Err(PyValueError::new_err(format!(
                "unknown bond order {name:?}; expected unspecified, single, double, triple or aromatic"
            )));
        }
    })
}

/// Bond atom-index pairs as a new ``(n_bonds, 2)`` integer array.
pub(crate) fn bonds_array<'py>(py: Python<'py>, bonds: &BondStorage) -> Bound<'py, PyArray2<usize>> {
    let arr = PyArray2::<usize>::zeros(py, [bonds.len(), 2], false);
    {
        let mut view = unsafe { arr.as_array_mut() };
        for (k, [i, j]) in bonds.iter_pairs().enumerate() {
            view[[k, 0]] = i;
            view[[k, 1]] = j;
        }
    }
    arr
}

/// Bond orders by name, or ``None`` when the source recorded no orders.
pub(crate) fn bond_order_names(bonds: &BondStorage) -> Option<Vec<&'static str>> {
    bonds
        .has_orders()
        .then(|| bonds.iter().map(|b| order_name(b.order())).collect())
}

/// Build a bond table from Python pairs and optional order names. Index
/// validation against the atom count is left to `Topology::set_bonds`.
pub(crate) fn bond_storage(
    pairs: &PyArrayLike2<'_, i64, AllowTypeChange>,
    orders: Option<Vec<String>>,
) -> PyResult<BondStorage> {
    let shape = pairs.shape();
    if shape[1] != 2 && shape[0] != 0 {
        return Err(PyValueError::new_err(format!(
            "bond pairs must have shape (n_bonds, 2), not {shape:?}"
        )));
    }
    let n = shape[0];
    let orders = match orders {
        Some(names) if names.len() != n => {
            return Err(PyValueError::new_err(format!(
                "{} bond orders for {n} bonds",
                names.len()
            )));
        }
        Some(names) => names.iter().map(|s| order_from_name(s)).collect::<PyResult<Vec<_>>>()?,
        None => vec![BondOrder::Unspecified; n],
    };
    let view = pairs.as_array();
    let mut bonds = BondStorage::default();
    bonds.reserve(n);
    for (k, order) in orders.into_iter().enumerate() {
        let (i, j) = (view[[k, 0]], view[[k, 1]]);
        if i < 0 || j < 0 {
            return Err(PyValueError::new_err(format!("negative atom index in bond ({i}, {j})")));
        }
        bonds.push(&Bond::with_order(i as usize, j as usize, order));
    }
    Ok(bonds)
}
/// Coordinate frame with simulation time and periodic box.
///
/// Stores coordinates, optional periodic box, and time for one frame.
///
/// **Example**
///
/// .. code-block:: python
///
///    import pymolar
///    fh = pymolar.FileHandler("traj.xtc", "r")
///    st = fh.read_state()
///    print(st.time)    # simulation time in ps
///    print(st.box)     # PeriodicBox
///    sys.state = st    # hot-swap coordinates into system

#[pyclass(name = "State", frozen)]
pub struct StatePy {
    pub(crate) st: UnsafeCell<State>,
    /// Monotonic generation counter, bumped on every in-place coordinate mutation
    /// (`translate` / `apply_transform` / `coords` setter / `Particle` position
    /// setters / …). An embedding viewer (molar_vis) polls it lock-free to re-render
    /// only when coordinates actually change. Unused — and ~free (one relaxed
    /// increment per mutation) — for standalone pymolar.
    coords_version: std::sync::atomic::AtomicU64,
}

unsafe impl Send for StatePy {}
unsafe impl Sync for StatePy {}

impl From<State> for StatePy {
    fn from(value: State) -> Self {
        Self {
            st: UnsafeCell::new(value),
            coords_version: std::sync::atomic::AtomicU64::new(0),
        }
    }
}

impl StatePy {
    pub(crate) fn inner(&self) -> &State {
        unsafe { &*self.st.get() }
    }

    pub(crate) fn inner_mut(&self) -> &mut State {
        unsafe { &mut *self.st.get() }
    }

    /// Bump the coordinate generation counter — call after any in-place coordinate
    /// edit. `Release` so a viewer that observes the new value (with `Acquire`) also
    /// sees the coordinate writes that preceded it.
    pub fn bump_coords_version(&self) {
        self.coords_version
            .fetch_add(1, std::sync::atomic::Ordering::Release);
    }

    /// The coordinate generation atomic, for an embedding viewer to poll lock-free
    /// (no GIL needed for the load).
    pub fn coords_version_atomic(&self) -> &std::sync::atomic::AtomicU64 {
        &self.coords_version
    }

    pub(crate) fn into_py(self) -> Py<StatePy> {
        Python::attach(|py| Py::new(py, self).unwrap())
    }
}

#[pymethods]
impl StatePy {
    /// Number of coordinates in this state.
    ///
    /// :returns: Number of atoms/coordinates.
    /// :rtype: int
    fn __len__(&self) -> usize {
        self.inner().len()
    }

    fn __repr__(&self) -> String {
        format!("State(n={}, t={:.2})", self.inner().len(), self.inner().time)
    }

    /// Periodic box of this state.
    ///
    /// :returns: Periodic box.
    /// :rtype: PeriodicBox
    #[getter]
    fn get_box(&self) -> PyResult<PeriodicBoxPy> {
        Ok(PeriodicBoxPy(
            self.inner()
                .pbox
                .as_ref()
                .ok_or_else(|| PyAttributeError::new_err("No periodic box to get"))?
                .clone(),
        ))
    }

    /// Set periodic box of this state.
    ///
    /// :param val: New periodic box.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter]
    fn set_box(&self, val: Bound<'_, PeriodicBoxPy>) -> PyResult<()> {
        let b = self
            .inner_mut()
            .pbox
            .as_mut()
            .ok_or_else(|| PyAttributeError::new_err("No periodic box to set"))?;
        *b = val.borrow().0.clone();
        Ok(())
    }

    /// Simulation time value.
    ///
    /// :returns: Frame time.
    /// :rtype: float
    #[getter]
    fn get_time(&self) -> Float {
        self.inner().time
    }

    /// Set simulation time value.
    ///
    /// :param t: New frame time.
    /// :returns: ``None``.
    /// :rtype: None
    #[setter]
    fn set_time(&self, t: Float) {
        self.inner_mut().time = t;
    }

    /// Copy periodic box from a ``System`` or ``Sel``.
    ///
    /// :param arg: Source object (``System`` or ``Sel``).
    /// :returns: ``None``.
    /// :rtype: None
    fn set_box_from(&self, arg: Bound<'_, PyAny>) -> PyResult<()> {
        let st_ref = match arg.cast::<SystemPy>() { Ok(sys) => {
            sys.get().r_st()
        } _ => { match arg.cast::<SelPy>() { Ok(sel) => {
            sel.get().r_st()
        } _ => {
            let ty_name = arg.get_type().name()?.to_string();
            return Err(PyTypeError::new_err(format!(
                "Invalid argument type {ty_name} in set_box_from()"
            )));
        }}}};
        self.inner_mut().pbox = st_ref.pbox.clone();
        Ok(())
    }
}

impl SaveState for StatePy {
    fn iter_pos_dyn<'a>(&'a self) -> Box<dyn ExactSizeIterator<Item = &'a Pos> + 'a> {
        self.inner().iter_pos_dyn()
    }
}

impl LenProvider for StatePy {
    fn len(&self) -> usize {
        self.__len__()
    }
}

impl TimeProvider for StatePy {
    fn get_time(&self) -> Float {
        self.inner().time
    }
}

impl BoxProvider for StatePy {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.inner().pbox.as_ref()
    }
}

impl IndexProvider for StatePy {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        i
    }
}

impl PosProvider for StatePy {
    unsafe fn coords_ptr(&self) -> *const Pos {
        self.inner().coords.as_ptr()
    }
}

//----------------------------------------------------------------
/// Molecular topology container (atoms and connectivity).
///
/// **Example**
///
/// .. code-block:: python
///
///    top = sys.topology
///    print(len(top))    # number of atoms

#[pyclass(name = "Topology", frozen)]
pub struct TopologyPy(pub(crate) UnsafeCell<Topology>);

unsafe impl Send for TopologyPy {}
unsafe impl Sync for TopologyPy {}

impl From<Topology> for TopologyPy {
    fn from(value: Topology) -> Self {
        Self(UnsafeCell::new(value))
    }
}

impl TopologyPy {
    pub(crate) fn inner(&self) -> &Topology {
        unsafe { &*self.0.get() }
    }

    pub(crate) fn inner_mut(&self) -> &mut Topology {
        unsafe { &mut *self.0.get() }
    }

    
    pub(crate) fn into_py(self) -> Py<TopologyPy> {
        Python::attach(|py| Py::new(py, self).unwrap())
    }
}

impl LenProvider for TopologyPy {
    fn len(&self) -> usize {
        self.inner().len()
    }
}

#[pymethods]
impl TopologyPy {
    /// Assign atom types to all topology atoms without coordinates.
    ///
    /// :param ff: FFType.Gaff, FFType.Gaff2, or a string "gaff" / "gaff2".
    /// :raises FFError: If bond orders or chemistry are invalid.
    ///
    /// Changes shared atom type_name fields. Preparation requires a System.
    #[pyo3(signature = (ff=FFArg::Name("gaff".to_owned())), text_signature = "($self, ff='gaff')")]
    fn apply_ff(&self, ff: FFArg) -> PyResult<()> {
        self.inner_mut().apply_ff(ff.resolve()?).map_err(to_ff_error)
    }

    /// Predict partial charges for all topology atoms without coordinates.
    ///
    /// :param model: ChargeModel.Espaloma or "espaloma".
    /// :raises ChargeError: If bond orders or elements are invalid.
    ///
    /// Changes shared partial charges; preserves the total formal charge.
    #[pyo3(signature = (model=ChargeArg::Name("espaloma".to_owned())), text_signature = "($self, model='espaloma')")]
    fn apply_charges(&self, model: ChargeArg) -> PyResult<()> {
        self.inner_mut().apply_charges(model.resolve()?).map_err(to_charge_error)
    }

    /// Number of atoms in this topology.
    fn __len__(&self) -> usize {
        self.inner().len()
    }

    fn __repr__(&self) -> String {
        format!("Topology(n={} atoms)", self.inner().len())
    }

    /// Bonds as a new ``(n_bonds, 2)`` integer array of atom indices.
    ///
    /// :returns: Copy of the bond table; editing it does not change the topology.
    /// :rtype: numpy.ndarray
    #[getter("bonds")]
    fn get_bonds<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<usize>> {
        bonds_array(py, &self.inner().bonds)
    }

    /// Bond orders, one per bond: ``"unspecified"``, ``"single"``, ``"double"``,
    /// ``"triple"`` or ``"aromatic"``.
    ///
    /// :returns: ``None`` when the source recorded no orders (PDB, GRO, TPR, ...).
    /// :rtype: list[str] | None
    #[getter("bond_orders")]
    fn get_bond_orders(&self) -> Option<Vec<&'static str>> {
        bond_order_names(&self.inner().bonds)
    }

    /// Replace the whole bond table.
    ///
    /// :param pairs: ``(n_bonds, 2)`` atom indices (array or sequence of pairs).
    /// :param orders: Optional order name per bond (see :attr:`bond_orders`).
    /// :raises ValueError: On a malformed table, an index out of range, or a self-bond.
    /// :returns: ``None``.
    /// :rtype: None
    #[pyo3(signature = (pairs, orders=None))]
    fn set_bonds(
        &self,
        pairs: PyArrayLike2<'_, i64, AllowTypeChange>,
        orders: Option<Vec<String>>,
    ) -> PyResult<()> {
        let bonds = bond_storage(&pairs, orders)?;
        self.inner_mut()
            .set_bonds(bonds)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }
}

impl BondProvider for TopologyPy {
    fn num_bonds(&self) -> usize {
        BondProvider::num_bonds(self.inner())
    }

    unsafe fn get_bond_unchecked(&self, i: usize) -> BondRef<'_> { unsafe {
        self.inner().get_bond_unchecked(i)
    }}

    fn iter_bonds(&self) -> impl Iterator<Item = BondRef<'_>> {
        self.inner().iter_bonds()
    }
}

impl SaveTopology for TopologyPy {
    fn iter_atoms_dyn<'a>(&'a self) -> Box<dyn Iterator<Item = AtomRef<'a>> + 'a> {
        self.inner().iter_atoms_dyn()
    }
    fn iter_bonds_dyn<'a>(&'a self) -> Box<dyn Iterator<Item = BondRef<'a>> + 'a> {
        self.inner().iter_bonds_dyn()
    }
    fn num_bonds(&self) -> usize {
        BondProvider::num_bonds(self)
    }
}
