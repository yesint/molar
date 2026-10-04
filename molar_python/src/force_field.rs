//! Python force-field choices, preparation options, and error conversion.

use molar::prelude::BondPerceptionError as NativeBondPerceptionError;
use molar::prelude::*;
use pyo3::{create_exception, exceptions::PyValueError, prelude::*, types::PyDict};

/// Force-field atom typing choice. Use ``FFType.Gaff`` or ``FFType.Gaff2``.
/// Existing string names ``"gaff"`` and ``"gaff2"`` also remain supported.
#[pyclass(name = "FFType", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum FFTypePy {
    /// General Amber Force Field.
    Gaff,
    /// General Amber Force Field 2.
    Gaff2,
}

impl From<FFTypePy> for molar_ff::FFType {
    fn from(value: FFTypePy) -> Self {
        match value {
            FFTypePy::Gaff => Self::Gaff,
            FFTypePy::Gaff2 => Self::Gaff2,
        }
    }
}

/// Partial-charge model. Use ``ChargeModel.Espaloma`` or ``"espaloma"``.
#[pyclass(name = "ChargeModel", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum ChargeModelPy {
    /// Espaloma-charge; preserves the total formal charge of atoms in scope.
    Espaloma,
}

impl From<ChargeModelPy> for molar_ff::ChargeModel {
    fn from(_: ChargeModelPy) -> Self {
        Self::Espaloma
    }
}

#[derive(FromPyObject)]
pub enum FFArg {
    #[pyo3(transparent)]
    Choice(FFTypePy),
    #[pyo3(transparent)]
    Name(String),
}

impl FFArg {
    pub(crate) fn resolve(&self) -> PyResult<molar_ff::FFType> {
        match self {
            Self::Choice(value) => Ok((*value).into()),
            Self::Name(name) => crate::utils::parse_ff(name),
        }
    }
}

#[derive(FromPyObject)]
pub enum ChargeArg {
    #[pyo3(transparent)]
    Choice(ChargeModelPy),
    #[pyo3(transparent)]
    Name(String),
}

impl ChargeArg {
    pub(crate) fn resolve(&self) -> PyResult<molar_ff::ChargeModel> {
        match self {
            Self::Choice(value) => Ok((*value).into()),
            Self::Name(name) => crate::utils::parse_charge_model(name),
        }
    }
}

/// Treatment of known input bond orders during preparation.
#[pyclass(name = "InputOrders", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum InputOrdersPy {
    /// Keep concrete input orders; solve unspecified bonds. Kekulize aromatic bonds.
    PreserveKnown,
    /// Reassign non-aromatic orders. Keep the result of aromatic kekulization.
    ReassignAll,
}

impl From<InputOrdersPy> for InputOrders {
    fn from(value: InputOrdersPy) -> Self {
        match value {
            InputOrdersPy::PreserveKnown => Self::PreserveKnown,
            InputOrdersPy::ReassignAll => Self::ReassignAll,
        }
    }
}

/// Hydrogen treatment during bond-order perception.
#[pyclass(name = "HydrogenPolicy", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum HydrogenPolicyPy {
    /// Require every hydrogen to be an explicit atom.
    AllExplicit,
    /// Infer missing hydrogens from geometry. Requires coordinates.
    InferFromGeometry,
}

impl From<HydrogenPolicyPy> for HydrogenPolicy {
    fn from(value: HydrogenPolicyPy) -> Self {
        match value {
            HydrogenPolicyPy::AllExplicit => Self::AllExplicit,
            HydrogenPolicyPy::InferFromGeometry => Self::InferFromGeometry,
        }
    }
}

/// Distance-based connectivity settings. Distances use nm.
/// Nested option objects are shared Python references; property edits take effect.
#[pyclass(name = "ConnectivityOptions", get_all, set_all)]
#[derive(Clone)]
pub struct ConnectivityOptionsPy {
    /// Added to the covalent-radius sum, in nm; finite and nonnegative.
    pub tolerance: Float,
    /// Minimum accepted bond distance, in nm; finite and nonnegative.
    pub minimum_distance: Float,
    /// Periodic dimensions in x, y, z order. A periodic box is required if enabled.
    pub pbc: [bool; 3],
    /// Remove the longest candidate bonds from over-coordinated atoms.
    pub cleanup_overcoordination: bool,
}

#[pymethods]
impl ConnectivityOptionsPy {
    #[new]
    #[pyo3(signature = (*, tolerance=0.045, minimum_distance=0.04, pbc=None, cleanup_overcoordination=true))]
    /// Create connectivity options. None for pbc disables periodic boundaries.
    fn new(
        tolerance: Float,
        minimum_distance: Float,
        pbc: Option<[bool; 3]>,
        cleanup_overcoordination: bool,
    ) -> Self {
        Self {
            tolerance,
            minimum_distance,
            pbc: pbc.unwrap_or([false; 3]),
            cleanup_overcoordination,
        }
    }
}

impl ConnectivityOptionsPy {
    pub(crate) fn options(&self) -> ConnectivityOptions {
        ConnectivityOptions {
            tolerance: self.tolerance,
            minimum_distance: self.minimum_distance,
            pbc: PbcDims::new(self.pbc[0], self.pbc[1], self.pbc[2]),
            cleanup_overcoordination: self.cleanup_overcoordination,
        }
    }
}

/// Work limits for each bonded fragment during bond-order search.
#[pyclass(name = "SearchLimits", get_all, set_all)]
#[derive(Clone)]
pub struct SearchLimitsPy {
    /// Maximum visited search nodes per fragment. Default: 5,000,000.
    pub max_branches: u64,
}

#[pymethods]
impl SearchLimitsPy {
    #[new]
    #[pyo3(signature = (*, max_branches=5_000_000))]
    /// Create search limits. An exhausted search reports an error or warning.
    fn new(max_branches: u64) -> Self {
        Self { max_branches }
    }
}

/// All bond-order settings used by preparation.
/// The standalone default reassigns orders, as in Rust BondOrderOptions.
/// PrepareOptions defaults instead preserve known orders.
#[pyclass(name = "BondOrderOptions", get_all, set_all)]
pub struct BondOrderOptionsPy {
    /// Preserve known input orders or reassign them.
    pub input_orders: InputOrdersPy,
    /// Require explicit hydrogens or infer them from coordinates.
    pub hydrogens: HydrogenPolicyPy,
    /// Use canonical bond orders for recognized functional groups.
    pub use_functional_groups: bool,
    /// Use canonical bond orders for recognized polymer residues.
    pub use_residue_templates: bool,
    /// Optional net formal charge; requires a single bonded fragment.
    pub total_charge: Option<i32>,
    /// Mutable search limits. Nested edits apply to this option object.
    pub limits: Py<SearchLimitsPy>,
}

#[pymethods]
impl BondOrderOptionsPy {
    #[new]
    #[pyo3(signature = (*, input_orders=InputOrdersPy::ReassignAll, hydrogens=HydrogenPolicyPy::AllExplicit, use_functional_groups=true, use_residue_templates=true, total_charge=None, limits=None))]
    /// Create bond-order options. Omitted limits use a separate default object.
    fn new(
        py: Python<'_>,
        input_orders: InputOrdersPy,
        hydrogens: HydrogenPolicyPy,
        use_functional_groups: bool,
        use_residue_templates: bool,
        total_charge: Option<i32>,
        limits: Option<Py<SearchLimitsPy>>,
    ) -> PyResult<Self> {
        Ok(Self {
            input_orders,
            hydrogens,
            use_functional_groups,
            use_residue_templates,
            total_charge,
            limits: match limits {
                Some(value) => value,
                None => Py::new(
                    py,
                    SearchLimitsPy::new(SearchLimits::default().max_branches),
                )?,
            },
        })
    }
}

impl BondOrderOptionsPy {
    pub(crate) fn options(&self, py: Python<'_>) -> BondOrderOptions {
        BondOrderOptions {
            input_orders: self.input_orders.into(),
            hydrogens: self.hydrogens.into(),
            use_functional_groups: self.use_functional_groups,
            use_residue_templates: self.use_residue_templates,
            total_charge: self.total_charge,
            limits: SearchLimits {
                max_branches: self.limits.borrow(py).max_branches,
            },
        }
    }
}

/// Options for explicit hydrogen addition during preparation.
#[pyclass(name = "HydrogenOptions", get_all, set_all)]
#[derive(Clone)]
pub struct HydrogenOptionsPy {
    /// Give new hydrogens zero velocities and forces if dynamics are present.
    /// If false, adding hydrogens to such a state raises BondPerceptionError.
    pub zero_fill_dynamics: bool,
}

#[pymethods]
impl HydrogenOptionsPy {
    #[new]
    #[pyo3(signature = (*, zero_fill_dynamics=true))]
    /// Create hydrogen-addition options.
    fn new(zero_fill_dynamics: bool) -> Self {
        Self { zero_fill_dynamics }
    }
}

/// Complete options for System.prepare_for_ff. All Rust PrepareOptions fields
/// are available. Set add_hydrogens to HydrogenOptions to enable addition.
/// Nested option properties are references, so nested edits take effect.
#[pyclass(name = "PrepareOptions", get_all, set_all)]
pub struct PrepareOptionsPy {
    /// Connectivity perception settings; used only when no bonds are present.
    pub connectivity: Py<ConnectivityOptionsPy>,
    /// Bond-order settings. Defaults preserve known input orders.
    pub bond_orders: Py<BondOrderOptionsPy>,
    /// Hydrogen addition settings, or None to leave hydrogens implicit.
    pub add_hydrogens: Option<Py<HydrogenOptionsPy>>,
}

#[pymethods]
impl PrepareOptionsPy {
    #[new]
    #[pyo3(signature = (*, connectivity=None, bond_orders=None, add_hydrogens=None))]
    /// Create complete preparation options. Default nested objects are independent
    /// for each construction. Explicit nested objects are shared, not copied.
    fn new(
        py: Python<'_>,
        connectivity: Option<Py<ConnectivityOptionsPy>>,
        bond_orders: Option<Py<BondOrderOptionsPy>>,
        add_hydrogens: Option<Py<HydrogenOptionsPy>>,
    ) -> PyResult<Self> {
        Ok(Self {
            connectivity: match connectivity {
                Some(value) => value,
                None => Py::new(py, ConnectivityOptionsPy::new(0.045, 0.04, None, true))?,
            },
            bond_orders: match bond_orders {
                Some(value) => value,
                None => Py::new(
                    py,
                    BondOrderOptionsPy::new(
                        py,
                        InputOrdersPy::PreserveKnown,
                        HydrogenPolicyPy::AllExplicit,
                        true,
                        true,
                        None,
                        None,
                    )?,
                )?,
            },
            add_hydrogens,
        })
    }
}

impl PrepareOptionsPy {
    pub(crate) fn options(&self, py: Python<'_>) -> molar_ff::PrepareOptions {
        molar_ff::PrepareOptions {
            connectivity: self.connectivity.borrow(py).options(),
            bond_orders: self.bond_orders.borrow(py).options(py),
            add_hydrogens: self.add_hydrogens.as_ref().map(|value| HydrogenOptions {
                zero_fill_dynamics: value.borrow(py).zero_fill_dynamics,
            }),
        }
    }
}

create_exception!(
    molar,
    FFError,
    PyValueError,
    "Atom typing failed. kind identifies the Rust error variant; details holds its fields."
);
create_exception!(
    molar,
    ChargeError,
    PyValueError,
    "Charge prediction failed. kind identifies the Rust error variant; details holds its fields."
);
create_exception!(
    molar,
    BondPerceptionError,
    PyValueError,
    "Preparation failed. kind identifies the Rust perception error variant; details holds its fields."
);

pub(crate) fn to_ff_error(error: molar_ff::FFError) -> PyErr {
    Python::attach(|py| -> PyResult<PyErr> {
        let details = PyDict::new(py);
        let kind = match &error {
            molar_ff::FFError::MissingBondOrders(i, j) => {
                details.set_item("atoms", (*i, *j))?;
                "MissingBondOrders"
            }
            molar_ff::FFError::OpenSelection { global, neighbor } => {
                details.set_item("atom", *global)?;
                details.set_item("neighbor", *neighbor)?;
                "OpenSelection"
            }
            molar_ff::FFError::InvalidAromatic(reason) => {
                details.set_item("reason", reason.to_string())?;
                "InvalidAromatic"
            }
            molar_ff::FFError::UntypedAtom { ff, local, z } => {
                details.set_item("ff", format!("{ff:?}").to_lowercase())?;
                details.set_item("local", *local)?;
                details.set_item("atomic_number", *z)?;
                "UntypedAtom"
            }
        };
        Ok(attach_error(
            py,
            FFError::new_err(error.to_string()),
            kind,
            details,
        ))
    })
    .unwrap_or_else(|error| error)
}

pub(crate) fn to_charge_error(error: molar_ff::ChargeError) -> PyErr {
    Python::attach(|py| -> PyResult<PyErr> {
        let details = PyDict::new(py);
        let kind = match &error {
            molar_ff::ChargeError::MissingBondOrders(i, j) => {
                details.set_item("atoms", (*i, *j))?;
                "MissingBondOrders"
            }
            molar_ff::ChargeError::OpenSelection { global, neighbor } => {
                details.set_item("atom", *global)?;
                details.set_item("neighbor", *neighbor)?;
                "OpenSelection"
            }
            molar_ff::ChargeError::Kekulize(reason) => {
                details.set_item("reason", reason.to_string())?;
                "Kekulize"
            }
            molar_ff::ChargeError::UnsupportedElement(z, model) => {
                details.set_item("atomic_number", *z)?;
                details.set_item("model", format!("{model:?}").to_lowercase())?;
                "UnsupportedElement"
            }
            molar_ff::ChargeError::Inference(reason) => {
                details.set_item("reason", reason)?;
                "Inference"
            }
        };
        Ok(attach_error(
            py,
            ChargeError::new_err(error.to_string()),
            kind,
            details,
        ))
    })
    .unwrap_or_else(|error| error)
}

pub(crate) fn to_perception_error(error: NativeBondPerceptionError) -> PyErr {
    Python::attach(|py| -> PyResult<PyErr> {
        let details = PyDict::new(py);
        let kind = match &error {
            NativeBondPerceptionError::DistanceSearch(reason) => {
                details.set_item("reason", reason.to_string())?;
                "DistanceSearch"
            }
            NativeBondPerceptionError::MissingPeriodicBox => "MissingPeriodicBox",
            NativeBondPerceptionError::UnsupportedElement {
                atom,
                atomic_number,
            } => {
                details.set_item("atom", *atom)?;
                details.set_item("atomic_number", *atomic_number)?;
                "UnsupportedElement"
            }
            NativeBondPerceptionError::InvalidTolerance(value) => {
                details.set_item("value", *value)?;
                "InvalidTolerance"
            }
            NativeBondPerceptionError::InvalidMinimumDistance(value) => {
                details.set_item("value", *value)?;
                "InvalidMinimumDistance"
            }
            NativeBondPerceptionError::AssignmentLength {
                field,
                expected,
                actual,
            } => {
                details.set_item("field", *field)?;
                details.set_item("expected", *expected)?;
                details.set_item("actual", *actual)?;
                "AssignmentLength"
            }
            NativeBondPerceptionError::AtomCountChanged { expected, actual } => {
                details.set_item("expected", *expected)?;
                details.set_item("actual", *actual)?;
                "AtomCountChanged"
            }
            NativeBondPerceptionError::AtomTableChanged {
                atom,
                expected,
                actual,
            } => {
                details.set_item("atom", *atom)?;
                details.set_item("expected", *expected)?;
                details.set_item("actual", *actual)?;
                "AtomTableChanged"
            }
            NativeBondPerceptionError::BondTableChanged {
                bond,
                expected,
                actual,
            } => {
                details.set_item("bond", *bond)?;
                details.set_item("expected", expected)?;
                details.set_item("actual", actual)?;
                "BondTableChanged"
            }
            NativeBondPerceptionError::Kekulization(reason) => {
                details.set_item("reason", reason)?;
                "Kekulization"
            }
            NativeBondPerceptionError::NoValidAssignment { atom } => {
                details.set_item("atom", *atom)?;
                "NoValidAssignment"
            }
            NativeBondPerceptionError::SearchLimitExceeded { atom } => {
                details.set_item("atom", *atom)?;
                "SearchLimitExceeded"
            }
            NativeBondPerceptionError::TotalChargeWithMultipleComponents => {
                "TotalChargeWithMultipleComponents"
            }
            NativeBondPerceptionError::MissingCoordinates => "MissingCoordinates",
            NativeBondPerceptionError::CoordinateCountMismatch => "CoordinateCountMismatch",
            NativeBondPerceptionError::DynamicsPresent => "DynamicsPresent",
        };
        Ok(attach_error(
            py,
            BondPerceptionError::new_err(error.to_string()),
            kind,
            details,
        ))
    })
    .unwrap_or_else(|error| error)
}

fn attach_error(py: Python<'_>, error: PyErr, kind: &str, details: Bound<'_, PyDict>) -> PyErr {
    // Python exception instances have dictionaries; these assignments cannot fail
    // except for interpreter allocation errors, which must propagate.
    let value = error.value(py);
    if let Err(error) = value
        .setattr("kind", kind)
        .and_then(|_| value.setattr("details", details))
    {
        return error;
    }
    error
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<FFTypePy>()?;
    m.add_class::<ChargeModelPy>()?;
    m.add_class::<InputOrdersPy>()?;
    m.add_class::<HydrogenPolicyPy>()?;
    m.add_class::<ConnectivityOptionsPy>()?;
    m.add_class::<SearchLimitsPy>()?;
    m.add_class::<BondOrderOptionsPy>()?;
    m.add_class::<HydrogenOptionsPy>()?;
    m.add_class::<PrepareOptionsPy>()?;
    m.add("FFError", m.py().get_type::<FFError>())?;
    m.add("ChargeError", m.py().get_type::<ChargeError>())?;
    m.add(
        "BondPerceptionError",
        m.py().get_type::<BondPerceptionError>(),
    )?;
    Ok(())
}
