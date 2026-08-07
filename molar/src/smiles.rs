//! **SMILES output** for any atom scope — a whole [`System`](crate::System), a
//! [`Topology`](crate::Topology), or a bound selection.
//!
//! ```no_run
//! use molar::prelude::*;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let sys = System::from_file("ligand.sdf")?;
//! println!("{}", sys.to_smiles()?);
//! # Ok(())
//! # }
//! ```
//!
//! # What this writes
//! A valid **Kekulé** SMILES: all-uppercase atoms with explicit alternating single/double bonds,
//! never the lowercase aromatic form. That is deliberate — it is always valid and never depends
//! on the consumer sharing molar's aromaticity model. Aromatic-order input (an SDF order-4
//! record, or anything that has been through [`System::perceive`](crate::System::perceive)) is
//! resolved by [`kekulize`](crate::perception::kekulize) on the way out.
//!
//! Hydrogens are folded into the atom they hang off, so an all-explicit-H MD structure comes out
//! as `CCO` rather than a bracketed mess. A hydrogen that cannot be folded — H₂, a lone H, a
//! bridging or charged H — is written `[H]`.
//!
//! # Deliberate limits
//! - **Not canonical.** Output is deterministic (traversal is driven by ascending atom index), so
//!   the same scope always yields the same string, but two different orderings of the same
//!   molecule need not agree. Do not use it as an identity key.
//! - **No stereochemistry.** Neither `@`/`@@` nor `/`,`\` are emitted, so the string is lossy for
//!   a chiral molecule. Take stereo from the 3D coordinates instead.
//! - **Bond orders must already be known.** PDB/GRO/XTC/TPR record no orders, so those inputs are
//!   rejected rather than guessed at — this module does not perceive bond orders.

mod valence;
mod write;

use std::collections::HashMap;
use std::fmt;
use thiserror::Error;

use crate::perception::KekulizeError;
use crate::prelude::*;

/// A written SMILES string plus the mapping back to molar atom indices.
#[derive(Debug, Clone)]
pub struct Smiles {
    string: String,
    atom_order: Vec<usize>,
}

impl Smiles {
    pub fn as_str(&self) -> &str {
        &self.string
    }

    /// Global atom index of each atom in the string, in written order — hydrogens folded into a
    /// neighbour are absent. This is what lets a per-atom result computed elsewhere on the SMILES
    /// (an ML prediction, say) be mapped back onto the atoms it came from.
    pub fn atom_order(&self) -> &[usize] {
        &self.atom_order
    }

    pub fn into_string(self) -> String {
        self.string
    }
}

impl fmt::Display for Smiles {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.string)
    }
}

impl AsRef<str> for Smiles {
    fn as_ref(&self) -> &str {
        &self.string
    }
}

/// Why a scope could not be written as SMILES. Atom indices are **global** throughout, as
/// elsewhere in molar.
#[derive(Debug, Error)]
pub enum SmilesError {
    #[error("nothing to write: no atoms in scope")]
    Empty,

    #[error(
        "bond {0}-{1} has no known order; SMILES needs explicit bond orders (use an SDF/mol \
         input, or install connectivity that has them via System::set_bonds)"
    )]
    UnspecifiedBondOrder(usize, usize),

    #[error(
        "selection is not bond-complete: atom {global} is bonded to non-selected atom \
         {neighbor}, and SMILES cannot express a dangling bond"
    )]
    OpenSelection { global: usize, neighbor: usize },

    #[error("atom {global} has no recognised element (atomic number {z})")]
    UnknownElement { global: usize, z: u8 },

    #[error("more than 99 ring-closure bonds are open at once, which SMILES cannot express")]
    TooManyRingClosures,

    #[error(transparent)]
    Kekulize(KekulizeError),
}

/// The molecule reduced to a local `0..n` subgraph, which is all the writer works on.
///
/// Private to this module; visible to the `write` child. Every array is indexed by **local**
/// atom (or bond) index, with `global` the way back out.
struct MolView<'a> {
    z: &'a [u8],
    fc: &'a [i32],
    /// Kekulé orders — no `Aromatic`, no `Unspecified`.
    orders: &'a [BondOrder],
    adj: &'a BondAdjacency,
    /// Hydrogens folded into each atom.
    h_count: &'a [u8],
    /// Atoms folded away and therefore not written.
    suppressed: &'a [bool],
    global: &'a [usize],
    /// The single input deciding traversal order. v1: the local index.
    rank: &'a [u32],
}

/// Write the atoms in scope as a SMILES string.
///
/// Blanket-implemented for everything that can read atoms and bonds — [`System`](crate::System),
/// [`Topology`](crate::Topology), and the bound selection types — exactly like
/// `molar_ff`'s `ApplyFF`. A bare detached [`Sel`] does not qualify (it carries no atom storage);
/// bind it first with [`System::bind`](crate::System::bind) or
/// [`select_bound`](crate::System::select_bound).
///
/// The scope is treated as **the molecule**: only bonds with both endpoints in scope are used, and
/// a bond leaving the scope is an error ([`SmilesError::OpenSelection`]) rather than a silently
/// dropped valence. Disconnected scopes are written as `.`-separated fragments.
pub trait ToSmiles {
    fn to_smiles(&self) -> Result<Smiles, SmilesError>;
}

impl<T: AtomProvider + BondProvider> ToSmiles for T {
    fn to_smiles(&self) -> Result<Smiles, SmilesError> {
        // 1. Map the scope onto a local 0..n subgraph. Same shape as `molar_ff`'s ApplyFF.
        let global: Vec<usize> = self.iter_index().collect();
        if global.is_empty() {
            return Err(SmilesError::Empty);
        }
        let n = global.len();
        let g2l: HashMap<usize, usize> =
            global.iter().enumerate().map(|(l, &g)| (g, l)).collect();

        let z: Vec<u8> = self.iter_atoms().map(|a| a.get_atomic_number()).collect();
        let fc: Vec<i32> =
            self.iter_atoms().map(|a| a.get_formal_charge().unwrap_or(0)).collect();

        // 2. Local bonds. On a bound selection `iter_bonds` yields the *whole* system's bonds, so
        //    filtering through `g2l` is what restricts them to the scope.
        let mut pairs: Vec<[usize; 2]> = Vec::new();
        let mut orders: Vec<BondOrder> = Vec::new();
        for b in self.iter_bonds() {
            let ([g1, g2], order) = (b.pair(), b.order());
            match (g2l.get(&g1).copied(), g2l.get(&g2).copied()) {
                (Some(i), Some(j)) => {
                    if order == BondOrder::Unspecified {
                        return Err(SmilesError::UnspecifiedBondOrder(g1, g2));
                    }
                    pairs.push([i, j]);
                    orders.push(order);
                }
                (Some(_), None) => {
                    return Err(SmilesError::OpenSelection { global: g1, neighbor: g2 })
                }
                (None, Some(_)) => {
                    return Err(SmilesError::OpenSelection { global: g2, neighbor: g1 })
                }
                (None, None) => {} // belongs to another molecule entirely
            }
        }

        let adj = BondAdjacency::build(n, pairs.iter().copied());

        // 3. Resolve aromatic orders into a Kekulé structure (a no-op if there are none).
        let orders = crate::perception::kekulize(&z, &fc, &orders, &adj)
            .map_err(|e| SmilesError::Kekulize(e.remap_atoms(&global)))?;
        debug_assert!(!orders.contains(&BondOrder::Aromatic));

        // 4. Fold hydrogens into their heavy neighbour. Anything that cannot be folded without
        //    losing information — H2, a lone H, a bridging H, a charged or multiply-bonded H —
        //    stays a written atom.
        let mut suppressed = vec![false; n];
        let mut h_count = vec![0u8; n];
        for i in 0..n {
            if z[i] != 1 || fc[i] != 0 {
                continue;
            }
            let nb = adj.neighbors(i);
            if nb.len() != 1 {
                continue;
            }
            let (partner, bond) = (nb[0].atom(), nb[0].bond());
            if z[partner] == 1 || orders[bond] != BondOrder::Single {
                continue;
            }
            suppressed[i] = true;
            h_count[partner] += 1;
        }

        // 5. Emit. Ranking by local index is what makes the output deterministic; canonical
        //    output later replaces this array and nothing else.
        let rank: Vec<u32> = (0..n as u32).collect();
        let (string, atom_order) = write::write(&MolView {
            z: &z,
            fc: &fc,
            orders: &orders,
            adj: &adj,
            h_count: &h_count,
            suppressed: &suppressed,
            global: &global,
            rank: &rank,
        })?;

        Ok(Smiles { string, atom_order })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use BondOrder::{Double as D, Single as S, Triple as T};

    /// A topology from atomic numbers + `(i, j, order)` bonds. `Topology` implements both
    /// provider traits, so the writer can be exercised without any coordinates.
    fn mol(z: &[u8], bonds: &[(usize, usize, BondOrder)]) -> Topology {
        let mut t = Topology::default();
        for &n in z {
            t.atoms.push(&Atom::new().with_atomic_number(n).with_resid(1));
        }
        for &(i, j, o) in bonds {
            t.bonds.push(&Bond::with_order(i, j, o));
        }
        t
    }

    /// Append `count` hydrogens bonded to heavy atom `to`, as a real molar structure carries them.
    fn add_h(z: &mut Vec<u8>, bonds: &mut Vec<(usize, usize, BondOrder)>, to: usize, count: usize) {
        for _ in 0..count {
            z.push(1);
            bonds.push((to, z.len() - 1, S));
        }
    }

    fn smiles_of(t: &Topology) -> String {
        t.to_smiles().expect("should write").into_string()
    }

    // ---- the organic subset written bare ----

    #[test]
    fn simple_acyclic_molecules() {
        // Methane: four explicit H all fold away.
        let mut z = vec![6];
        let mut b = vec![];
        add_h(&mut z, &mut b, 0, 4);
        assert_eq!(smiles_of(&mol(&z, &b)), "C");

        // Ethene.
        let mut z = vec![6, 6];
        let mut b = vec![(0, 1, D)];
        add_h(&mut z, &mut b, 0, 2);
        add_h(&mut z, &mut b, 1, 2);
        assert_eq!(smiles_of(&mol(&z, &b)), "C=C");

        // Acetylene, exercising the triple-bond symbol.
        let mut z = vec![6, 6];
        let mut b = vec![(0, 1, T)];
        add_h(&mut z, &mut b, 0, 1);
        add_h(&mut z, &mut b, 1, 1);
        assert_eq!(smiles_of(&mol(&z, &b)), "C#C");

        // Ethanol: a heteroatom in the chain, and an O-H that folds away too.
        let mut z = vec![6, 6, 8];
        let mut b = vec![(0, 1, S), (1, 2, S)];
        add_h(&mut z, &mut b, 0, 3);
        add_h(&mut z, &mut b, 1, 2);
        add_h(&mut z, &mut b, 2, 1);
        assert_eq!(smiles_of(&mol(&z, &b)), "CCO");
    }

    /// A branch: the carboxyl carbon has two children, so the first is parenthesised and the
    /// second continues the chain.
    #[test]
    fn branches_are_parenthesised() {
        let mut z = vec![6, 6, 8, 8];
        let mut b = vec![(0, 1, S), (1, 2, D), (1, 3, S)];
        add_h(&mut z, &mut b, 0, 3);
        add_h(&mut z, &mut b, 3, 1);
        assert_eq!(smiles_of(&mol(&z, &b)), "CC(=O)O", "acetic acid");
    }

    // ---- brackets: charge, and hydrogen counts the notation cannot imply ----

    #[test]
    fn charged_atoms_use_the_bracket_form() {
        // Acetate: the anionic O carries a charge, so it brackets; the carbonyl O does not.
        let mut z = vec![6, 6, 8, 8];
        let mut b = vec![(0, 1, S), (1, 2, D), (1, 3, S)];
        add_h(&mut z, &mut b, 0, 3);
        let mut t = mol(&z, &b);
        t.atoms.get_mut(3).unwrap().set_formal_charge(-1);
        assert_eq!(smiles_of(&t), "CC(=O)[O-]");

        // Ammonium: charge plus a hydrogen count, both stated explicitly.
        let mut z = vec![7];
        let mut b = vec![];
        add_h(&mut z, &mut b, 0, 4);
        let mut t = mol(&z, &b);
        t.atoms.get_mut(0).unwrap().set_formal_charge(1);
        assert_eq!(smiles_of(&t), "[NH4+]");
    }

    /// A radical: three hydrogens where the notation would imply four, so the count has to be
    /// spelled out rather than left to the default valence.
    #[test]
    fn a_mismatched_hydrogen_count_forces_brackets() {
        let mut z = vec![6];
        let mut b = vec![];
        add_h(&mut z, &mut b, 0, 3);
        assert_eq!(smiles_of(&mol(&z, &b)), "[CH3]");
    }

    /// Hydrogens that cannot be folded into a heavy neighbour stay written.
    #[test]
    fn unfoldable_hydrogens_are_written_as_bracket_atoms() {
        assert_eq!(smiles_of(&mol(&[1, 1], &[(0, 1, S)])), "[H][H]", "H2 has no heavy atom");
        assert_eq!(smiles_of(&mol(&[1], &[])), "[H]", "a lone hydrogen");
    }

    /// Elements outside the organic subset always bracket, hydrogens or not.
    #[test]
    fn non_subset_elements_bracket() {
        let mut t = mol(&[11], &[]); // Na+
        t.atoms.get_mut(0).unwrap().set_formal_charge(1);
        assert_eq!(smiles_of(&t), "[Na+]");

        let mut z = vec![14]; // silane
        let mut b = vec![];
        add_h(&mut z, &mut b, 0, 4);
        assert_eq!(smiles_of(&mol(&z, &b)), "[SiH4]");
    }

    // ---- rings ----

    fn benzene_kekule() -> Topology {
        let mut z = vec![6; 6];
        let mut b = vec![(0, 1, D), (1, 2, S), (2, 3, D), (3, 4, S), (4, 5, D), (5, 0, S)];
        for i in 0..6 {
            add_h(&mut z, &mut b, i, 1);
        }
        mol(&z, &b)
    }

    #[test]
    fn ring_closure_numbering() {
        assert_eq!(smiles_of(&benzene_kekule()), "C1=CC=CC=C1");

        // Cyclohexane: every ring bond single, so no bond symbols at all.
        let mut z = vec![6; 6];
        let mut b = vec![(0, 1, S), (1, 2, S), (2, 3, S), (3, 4, S), (4, 5, S), (5, 0, S)];
        for i in 0..6 {
            add_h(&mut z, &mut b, i, 2);
        }
        assert_eq!(smiles_of(&mol(&z, &b)), "C1CCCCC1");
    }

    /// Aromatic-order input has no Kekulé structure of its own, so the writer has to build one.
    /// Which of benzene's two resonance forms comes out is arbitrary — that it is alternating is
    /// not.
    #[test]
    fn aromatic_input_is_kekulized_on_the_way_out() {
        let mut t = benzene_kekule();
        crate::perception::perceive(&mut t); // aromatizes every ring bond
        let s = smiles_of(&t);
        assert_eq!(s.matches('=').count(), 3, "benzene needs three double bonds: {s}");
        assert!(!s.contains('c'), "v1 never writes lowercase aromatic atoms: {s}");
        // Same shape as the Kekulé input, up to which bonds carry the doubles.
        assert_eq!(s.len(), "C1=CC=CC=C1".len(), "{s}");
    }

    /// Ring-closure numbers past 9 need the `%nn` form. Ten bonds from atom 0 to the far end of a
    /// chain all open at atom 0 and close one at a time, so ten are open at once.
    #[test]
    fn ring_closure_numbers_past_nine_use_the_percent_form() {
        let n = 15;
        let mut b: Vec<(usize, usize, BondOrder)> = (0..n - 1).map(|i| (i, i + 1, S)).collect();
        for far in 5..n {
            b.push((0, far, S));
        }
        let s = smiles_of(&mol(&vec![6; n], &b));
        assert!(s.contains("%10"), "a tenth simultaneous ring closure needs %10: {s}");
        // Every number opened is closed again: each digit run appears exactly twice.
        for k in 1..=9u32 {
            let c = s.matches(&k.to_string()).count();
            assert!(c % 2 == 0, "ring number {k} appears {c} times in {s}");
        }
    }

    // ---- fragments ----

    #[test]
    fn disconnected_scopes_are_dot_separated() {
        // Methane and water, unbonded.
        let mut z = vec![6];
        let mut b = vec![];
        add_h(&mut z, &mut b, 0, 4);
        z.push(8);
        let o = z.len() - 1;
        add_h(&mut z, &mut b, o, 2);
        assert_eq!(smiles_of(&mol(&z, &b)), "C.O");
    }

    // ---- the mapping back to molar indices ----

    #[test]
    fn atom_order_maps_written_atoms_back_to_global_indices() {
        let mut z = vec![6, 6, 8];
        let mut b = vec![(0, 1, S), (1, 2, S)];
        add_h(&mut z, &mut b, 0, 3);
        add_h(&mut z, &mut b, 1, 2);
        add_h(&mut z, &mut b, 2, 1);
        let s = mol(&z, &b).to_smiles().unwrap();
        assert_eq!(s.as_str(), "CCO");
        assert_eq!(s.atom_order(), &[0, 1, 2], "folded hydrogens are absent");
    }

    /// On a selection the mapping must carry *global* indices, not selection-local ones.
    #[test]
    fn atom_order_is_global_for_a_selection() {
        // Water (0,1,2) then ethanol-ish fragment; select only the second fragment.
        let mut z = vec![8];
        let mut b = vec![];
        add_h(&mut z, &mut b, 0, 2); // atoms 1,2
        let c0 = z.len();
        z.push(6);
        add_h(&mut z, &mut b, c0, 4); // methane at 3, H 4..7
        let mut t = mol(&z, &b);
        t.assign_resindex();
        let n = t.atoms.len();
        let sys =
            System::new(t, State { coords: vec![Pos::default(); n], ..Default::default() }).unwrap();

        let sel = sys.select_bound(c0..n).unwrap();
        let s = sel.to_smiles().unwrap();
        assert_eq!(s.as_str(), "C");
        assert_eq!(s.atom_order(), &[c0], "the methane carbon's global index");
    }

    // ---- errors ----

    #[test]
    fn unspecified_bond_orders_are_rejected_not_guessed() {
        let t = mol(&[6, 6], &[(0, 1, BondOrder::Unspecified)]);
        assert!(matches!(
            t.to_smiles(),
            Err(SmilesError::UnspecifiedBondOrder(0, 1))
        ));
    }

    /// A selection that cuts through a bond has a dangling valence, which SMILES cannot express.
    #[test]
    fn a_selection_cutting_a_bond_is_rejected() {
        let mut z = vec![6, 6, 8];
        let mut b = vec![(0, 1, S), (1, 2, S)];
        add_h(&mut z, &mut b, 0, 3);
        add_h(&mut z, &mut b, 1, 2);
        add_h(&mut z, &mut b, 2, 1);
        let mut t = mol(&z, &b);
        t.assign_resindex();
        let n = t.atoms.len();
        let sys =
            System::new(t, State { coords: vec![Pos::default(); n], ..Default::default() }).unwrap();

        // Just the two carbons: the C-O bond leaves the scope.
        let sel = sys.select_bound(0..2).unwrap();
        assert!(matches!(
            sel.to_smiles(),
            Err(SmilesError::OpenSelection { global: 1, neighbor: 2 })
        ));
    }

    #[test]
    fn an_empty_topology_has_nothing_to_write() {
        assert!(matches!(Topology::default().to_smiles(), Err(SmilesError::Empty)));
    }

    #[test]
    fn an_unknown_element_is_reported() {
        let mut z = vec![0]; // never assigned an element
        let mut b = vec![];
        add_h(&mut z, &mut b, 0, 1);
        assert!(matches!(
            mol(&z, &b).to_smiles(),
            Err(SmilesError::UnknownElement { global: 0, z: 0 })
        ));
    }
}
