//! Force-field atom typing for MolAR.
//!
//! This crate assigns GAFF (and, later, GAFF2) atom types to a molecule. Typing is
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

/// The force field whose atom types should be assigned.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FFType {
    /// General Amber Force Field (GAFF, `gaff.dat` / `ATOMTYPE_GFF.DEF`).
    Gaff,
    /// General Amber Force Field 2 (GAFF2). Not implemented yet.
    Gaff2,
}

/// Errors returned by [`ApplyFF::apply_ff`].
#[derive(Debug, thiserror::Error)]
pub enum FFError {
    /// GAFF2 typing is not implemented yet.
    #[error("GAFF2 typing is not implemented yet")]
    Gaff2NotImplemented,

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
        // GAFF2 not yet supported — fail before doing any work.
        if ff == FFType::Gaff2 {
            return Err(FFError::Gaff2NotImplemented);
        }

        // 1. Build a global -> local index map. For `System`/`Topology` this is the
        //    identity (0..n); for a bound selection it maps the sorted global indices.
        let global: Vec<usize> = self.iter_index().collect();
        let g2l: HashMap<usize, usize> =
            global.iter().enumerate().map(|(l, &g)| (g, l)).collect();

        // 2. Local atoms (atomic numbers) in local (== ascending-global) order.
        let z: Vec<u8> = self.iter_atoms().map(|a| a.atomic_number).collect();

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
        let types = match ff {
            FFType::Gaff => gaff::gaff_types(&z, &bonds)?,
            FFType::Gaff2 => unreachable!("handled above"),
        };

        // 5. Write results back in local order.
        for (a, t) in self.iter_atoms_mut().zip(types) {
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
