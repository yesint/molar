//! Opt-in preparation for force-field typing and charging (delivery step 8 of the
//! bond-perception plan).
//!
//! [`ApplyFF`](crate::ApplyFF) and [`ApplyCharges`](crate::ApplyCharges) are deliberately
//! strict: they require a molecule that already carries bond orders, and reject an order-less
//! input (PDB/GRO/XYZ). [`PrepareForFF::prepare_for_ff`] is the opt-in bridge that runs MolAR's
//! perception pipeline so such an input can be typed and charged:
//!
//! 1. perceive connectivity from coordinates when the topology has no bonds;
//! 2. perceive bond orders and formal charges;
//! 3. optionally add the missing hydrogens.
//!
//! The perception logic itself lives in `molar`; this wrapper only orders the steps for the
//! force-field workflow and leaves the low-level typing/charging operations untouched.
//!
//! ```no_run
//! use molar::prelude::*;
//! use molar_ff::{ApplyFF, FFType, PrepareForFF, PrepareOptions};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let mut sys = System::from_file("ligand.pdb")?; // no bond records
//! sys.prepare_for_ff(&PrepareOptions::default())?; // perceive bonds/orders/charges
//! sys.apply_ff(FFType::Gaff)?;                     // now typeable
//! # Ok(())
//! # }
//! ```

use molar::prelude::*;

/// How to prepare a molecule for force-field typing and charging.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrepareOptions {
    /// Distance-based connectivity perception, used only when the topology has no bonds.
    pub connectivity: ConnectivityOptions,
    /// Bond-order and formal-charge perception. Defaults to [`InputOrders::PreserveKnown`] so any
    /// orders already on the input are kept and only the unspecified bonds are solved.
    pub bond_orders: BondOrderOptions,
    /// When set, add the perceived implicit hydrogens as explicit atoms with these options.
    pub add_hydrogens: Option<HydrogenOptions>,
}

impl Default for PrepareOptions {
    fn default() -> Self {
        Self {
            connectivity: ConnectivityOptions::default(),
            bond_orders: BondOrderOptions {
                input_orders: InputOrders::PreserveKnown,
                ..BondOrderOptions::default()
            },
            add_hydrogens: None,
        }
    }
}

/// Prepare a molecule for the strict [`ApplyFF`](crate::ApplyFF) /
/// [`ApplyCharges`](crate::ApplyCharges) operations (see the module docs).
pub trait PrepareForFF {
    fn prepare_for_ff(&mut self, options: &PrepareOptions) -> Result<(), BondPerceptionError>;
}

impl PrepareForFF for System {
    fn prepare_for_ff(&mut self, options: &PrepareOptions) -> Result<(), BondPerceptionError> {
        if self.topology().bonds.is_empty() {
            self.perceive_connectivity(&options.connectivity)?;
        }
        let assignment = self.assign_bond_orders(&options.bond_orders)?;
        self.apply_bond_assignment(&assignment)?;
        if let Some(hydrogen_options) = options.add_hydrogens {
            let plan = plan_hydrogen_addition(self, &hydrogen_options);
            self.add_hydrogens(&plan)?;
        }
        Ok(())
    }
}
