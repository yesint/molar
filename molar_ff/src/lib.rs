//! Force-field atom typing for MolAR.
//!
//! This crate assigns GAFF and GAFF2 atom types to a molecule. Typing is
//! exposed through the [`ApplyFF`] trait, which is blanket-implemented for every MolAR
//! object that exposes atoms and bonds — [`System`], [`Topology`], and the bound
//! selection types.
//!
//! ```no_run
//! use molar::prelude::*;
//! use molar_ff::{ApplyFF, FFType};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let mut sys = System::from_file("ligand.sdf")?;
//! sys.apply_ff(FFType::Gaff)?;                                  // whole system
//! sys.select_bound_mut("resname LIG")?.apply_ff(FFType::Gaff)?; // just a selection
//! # Ok(())
//! # }
//! ```
//!
//! # Input requirement
//! Molecules must already carry bond orders (SDF/mol2). Inputs without bond orders
//! (PDB/GRO) return [`FFError::MissingBondOrders`] — this crate does not perceive bond
//! orders itself.

use std::collections::HashMap;

use molar::prelude::*;

mod gaff;

#[cfg(feature = "espaloma")]
pub mod charge;

/// The force field whose atom types should be assigned.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FFType {
    /// General Amber Force Field (GAFF, `gaff.dat` / `ATOMTYPE_GFF.DEF`).
    Gaff,
    /// General Amber Force Field 2 (GAFF2, `gaff2.dat` / `ATOMTYPE_GFF2.DEF`).
    Gaff2,
}

/// Errors returned by [`ApplyFF::apply_ff`].
#[derive(Debug, thiserror::Error)]
pub enum FFError {
    /// A bond in scope has no known order. Force-field typing requires known bond
    /// orders; supply an SDF/mol2 input (this crate does not perceive bond orders).
    #[error(
        "bond {0}-{1} has no known order; force-field typing requires known bond \
         orders (use an SDF/mol2 input)"
    )]
    MissingBondOrders(usize, usize),

    /// The selection is not bond-complete: a selected atom is bonded to an atom
    /// outside the selection, so it cannot be typed as a whole molecule.
    #[error(
        "selection is not bond-complete: atom {global} is bonded to non-selected \
         atom {neighbor}"
    )]
    OpenSelection { global: usize, neighbor: usize },

    /// No matching rule assigned a type to this atom.
    #[error("could not assign a {ff:?} atom type to atom {local} (element Z={z})")]
    UntypedAtom { ff: FFType, local: usize, z: u8 },
}

/// Assign force-field atom types, writing the result into each atom's `type_name`.
///
/// Blanket-implemented for every type that can both read atoms + bonds and write atoms
/// ([`AtomMutProvider`] + [`BondProvider`]): [`System`], [`Topology`], and the mutable
/// bound selection types. A bare detached [`Sel`] does not qualify — bind it first
/// (e.g. `sys.try_bind_mut(&sel)`).
///
/// The object is treated as **the molecule**: perception (rings/aromaticity/conjugation)
/// runs over the atoms in scope using only bonds whose both endpoints are in scope. For
/// a whole [`System`] this is the entire topology; for a selection it is the selected
/// atoms (which should form complete molecule(s) — see [`FFError::OpenSelection`]).
pub trait ApplyFF {
    fn apply_ff(&mut self, ff: FFType) -> Result<(), FFError>;
}

impl<T: AtomMutProvider + BondProvider> ApplyFF for T {
    fn apply_ff(&mut self, ff: FFType) -> Result<(), FFError> {
        // 1. Build a global -> local index map. For `System`/`Topology` this is the
        //    identity (0..n); for a bound selection it maps the sorted global indices.
        let global: Vec<usize> = self.iter_index().collect();
        let g2l: HashMap<usize, usize> =
            global.iter().enumerate().map(|(l, &g)| (g, l)).collect();

        // 2. Local atoms (atomic numbers) in local (== ascending-global) order.
        let z: Vec<u8> = self.iter_atoms().map(|a| a.get_atomic_number()).collect();

        // 3. Local bonds: keep only bonds whose both endpoints are in scope, remapped
        //    to local indices; error on unknown orders and boundary-crossing bonds.
        let mut bonds: Vec<gaff::LocalBond> = Vec::new();
        for b in self.iter_bonds() {
            match (g2l.get(&b.i1).copied(), g2l.get(&b.i2).copied()) {
                (Some(i), Some(j)) => {
                    let order = bond_order_code(b.order)
                        .ok_or(FFError::MissingBondOrders(b.i1, b.i2))?;
                    bonds.push(gaff::LocalBond { i, j, order });
                }
                (Some(_), None) => {
                    return Err(FFError::OpenSelection { global: b.i1, neighbor: b.i2 })
                }
                (None, Some(_)) => {
                    return Err(FFError::OpenSelection { global: b.i2, neighbor: b.i1 })
                }
                (None, None) => {} // bond belongs to another molecule; ignore
            }
        }

        // 4. Compute the types on the local subgraph.
        let types = gaff::gaff_types(&z, &bonds, ff)?;

        // 5. Write results back in local order.
        for (mut a, t) in self.iter_atoms_mut().zip(types) {
            a.set_type_name(&t);
        }
        Ok(())
    }
}

/// Map a MolAR [`BondOrder`] to the integer bond-type code used by the typing engine
/// (1 single, 2 double, 3 triple, 4 aromatic). `Unspecified` yields `None`
/// (→ [`FFError::MissingBondOrders`]).
fn bond_order_code(o: BondOrder) -> Option<u8> {
    match o {
        BondOrder::Single => Some(1),
        BondOrder::Double => Some(2),
        BondOrder::Triple => Some(3),
        BondOrder::Aromatic => Some(4),
        BondOrder::Unspecified => None,
    }
}

/// The partial-charge model used by [`ApplyCharges::apply_charges`].
#[cfg(feature = "espaloma")]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChargeModel {
    /// espaloma-charge GNN (a bundled ONNX). Charges are equilibrated to sum to zero over the
    /// atoms in scope (espaloma's whole-graph convention — it does not split by fragment).
    Espaloma,
}

/// Errors returned by [`ApplyCharges::apply_charges`].
#[cfg(feature = "espaloma")]
#[derive(Debug, thiserror::Error)]
pub enum ChargeError {
    /// A bond in scope lacks an explicit Kekulé order. Charge prediction needs single/double/
    /// triple bonds (supply an SDF/mol2 input); aromatic-order or unspecified bonds are rejected.
    #[error(
        "bond {0}-{1} has no explicit Kekulé order; charge prediction requires single/double/\
         triple bonds (use an SDF/mol2 input)"
    )]
    MissingBondOrders(usize, usize),

    /// A selected atom is bonded to an atom outside the selection, so the scope is not a
    /// complete molecule and cannot be charged as one.
    #[error(
        "selection is not bond-complete: atom {global} is bonded to non-selected atom {neighbor}"
    )]
    OpenSelection { global: usize, neighbor: usize },

    /// The molecule contains an element the charge model was not trained on.
    #[error("element Z={0} is not supported by the {1:?} charge model")]
    UnsupportedElement(u8, ChargeModel),

    /// The underlying model failed to run.
    #[error("charge model inference failed: {0}")]
    Inference(String),
}

/// Assign partial charges, writing the result into each atom's `charge`.
///
/// Blanket-implemented for every type that can read atoms + bonds and write atoms
/// ([`AtomMutProvider`] + [`BondProvider`]): [`System`], [`Topology`], and the mutable bound
/// selection types — exactly like [`ApplyFF`]. The object in scope is treated as **the
/// molecule** (perception runs over the atoms in scope using only bonds whose endpoints are both
/// in scope); a selection should be bond-complete (see [`ChargeError::OpenSelection`]).
///
/// Any integer formal charge already on an atom's `charge` (e.g. read from an SDF `M  CHG`
/// record) is consumed as input to the featurization and then **overwritten** with the predicted
/// partial charge.
///
/// ```no_run
/// use molar::prelude::*;
/// use molar_ff::{ApplyCharges, ChargeModel};
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let mut sys = System::from_file("ligand.sdf")?;
/// sys.apply_charges(ChargeModel::Espaloma)?;
/// # Ok(())
/// # }
/// ```
#[cfg(feature = "espaloma")]
pub trait ApplyCharges {
    fn apply_charges(&mut self, model: ChargeModel) -> Result<(), ChargeError>;
}

#[cfg(feature = "espaloma")]
impl<T: AtomMutProvider + BondProvider> ApplyCharges for T {
    fn apply_charges(&mut self, model: ChargeModel) -> Result<(), ChargeError> {
        match model {
            ChargeModel::Espaloma => {}
        }

        // 1. global -> local index map (identity for System/Topology; sorted globals for a sel).
        let global: Vec<usize> = self.iter_index().collect();
        let g2l: HashMap<usize, usize> =
            global.iter().enumerate().map(|(l, &g)| (g, l)).collect();

        // 2. Atomic numbers + integer formal charges in local order.
        let z: Vec<u8> = self.iter_atoms().map(|a| a.get_atomic_number()).collect();
        if let Some(&bad) = z.iter().find(|&&zi| {
            !matches!(zi, 1 | 6 | 7 | 8 | 9 | 15 | 16 | 17 | 35 | 53)
        }) {
            return Err(ChargeError::UnsupportedElement(bad, model));
        }
        let fc: Vec<i32> = self.iter_atoms().map(|a| a.get_formal_charge().unwrap_or(0)).collect();

        // 3. Local Kekulé bonds; error on boundary-crossing or non-Kekulé bonds.
        let mut bonds: Vec<gaff::LocalBond> = Vec::new();
        for b in self.iter_bonds() {
            match (g2l.get(&b.i1).copied(), g2l.get(&b.i2).copied()) {
                (Some(i), Some(j)) => {
                    let order = match b.order {
                        BondOrder::Single => 1,
                        BondOrder::Double => 2,
                        BondOrder::Triple => 3,
                        _ => return Err(ChargeError::MissingBondOrders(b.i1, b.i2)),
                    };
                    bonds.push(gaff::LocalBond { i, j, order });
                }
                (Some(_), None) => {
                    return Err(ChargeError::OpenSelection { global: b.i1, neighbor: b.i2 })
                }
                (None, Some(_)) => {
                    return Err(ChargeError::OpenSelection { global: b.i2, neighbor: b.i1 })
                }
                (None, None) => {}
            }
        }

        // 4. Predict, then write the partial charges back in local order.
        let q = charge::espaloma_charges(&z, &fc, &bonds)
            .map_err(|e| ChargeError::Inference(e.to_string()))?;
        for (mut a, c) in self.iter_atoms_mut().zip(q) {
            a.set_charge(c as Float);
        }
        Ok(())
    }
}
