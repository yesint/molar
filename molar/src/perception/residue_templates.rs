//! Residue templates (delivery step 6 of the bond-perception plan).
//!
//! A standard biopolymer residue has a known chemical graph, so its non-single bonds can be
//! pinned directly instead of searched. This matters at polymer scale: a peptide's repeating
//! amide backbone gives the general search many equal-cost charge-separated resonance forms, and
//! it can settle on a wrong one. Pinning each residue's backbone carbonyl (and its side-chain
//! double bonds) removes that ambiguity, and the remaining bonds — the peptide links between
//! residues, and any ligand — are left to the general solver.
//!
//! Templates match on residue name plus PDB atom names, and pin **orders** only; formal charges
//! then follow from the valence model, exactly as for [`functional_groups`](super::functional_groups).
//! A residue name that is not a standard amino acid matches nothing. Aromatic side-chain rings
//! (Phe/Tyr/Trp/His) are intentionally left to the general search, which resolves them from the
//! explicit ring hydrogens without a tautomer assumption.

use crate::prelude::*;

/// Canonical non-single bond orders for standard-residue bonds, indexed by bond. `None` means no
/// template applies. Only bonds whose two atoms share a residue are considered, so inter-residue
/// peptide/nucleotide links stay single via the general solver.
pub(super) fn residue_template_orders(top: &Topology, pairs: &[[usize; 2]]) -> Vec<Option<BondOrder>> {
    let mut orders: Vec<Option<BondOrder>> = vec![None; pairs.len()];
    let atoms = &top.atoms;
    for (b, &[i, j]) in pairs.iter().enumerate() {
        let ai = atoms.get(i).unwrap();
        let aj = atoms.get(j).unwrap();
        if ai.get_resindex() != aj.get_resindex() {
            continue; // an inter-residue link
        }
        if let Some(o) = residue_bond_order(ai.get_resname(), ai.get_name(), aj.get_name()) {
            orders[b] = Some(o);
        }
    }
    orders
}

/// The standard amino acids, including the common protonation/tautomer and terminal variants
/// that share the same heavy-atom graph.
fn is_amino_acid(resname: &str) -> bool {
    matches!(
        resname,
        "ALA" | "ARG" | "ASN" | "ASP" | "CYS" | "CYX" | "GLN" | "GLU" | "GLY" | "HIS" | "HID"
            | "HIE" | "HIP" | "HSD" | "HSE" | "HSP" | "ILE" | "LEU" | "LYS" | "MET" | "PHE"
            | "PRO" | "SER" | "THR" | "TRP" | "TYR" | "VAL"
    )
}

/// The canonical order of the bond between two named atoms of a residue, or `None` if the
/// residue/pair is not templated (then the bond is single or solved).
fn residue_bond_order(resname: &str, name_a: &str, name_b: &str) -> Option<BondOrder> {
    let is = |x: &str, y: &str| (name_a == x && name_b == y) || (name_a == y && name_b == x);

    // Backbone carbonyl, common to every amino acid.
    if is_amino_acid(resname) && is("C", "O") {
        return Some(BondOrder::Double);
    }
    // Side-chain double bonds. Carboxylate/carboxyl oxygens are also caught by the general
    // carboxyl group pass, but naming them here keeps the residue chemistry explicit.
    let double = match resname {
        "ASP" => is("CG", "OD1"),
        "GLU" => is("CD", "OE1"),
        "ASN" => is("CG", "OD1"), // side-chain amide
        "GLN" => is("CD", "OE1"), // side-chain amide
        "ARG" => is("CZ", "NH1"), // guanidinium
        _ => false,
    };
    double.then_some(BondOrder::Double)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn backbone_carbonyl_is_double_for_amino_acids() {
        assert_eq!(residue_bond_order("ALA", "C", "O"), Some(BondOrder::Double));
        assert_eq!(residue_bond_order("GLY", "O", "C"), Some(BondOrder::Double));
        // Not an amino acid: no template.
        assert_eq!(residue_bond_order("LIG", "C", "O"), None);
        // A different backbone bond is not the carbonyl.
        assert_eq!(residue_bond_order("ALA", "CA", "C"), None);
    }

    #[test]
    fn side_chain_double_bonds() {
        assert_eq!(residue_bond_order("ASP", "CG", "OD1"), Some(BondOrder::Double));
        assert_eq!(residue_bond_order("ASN", "OD1", "CG"), Some(BondOrder::Double));
        assert_eq!(residue_bond_order("ARG", "CZ", "NH1"), Some(BondOrder::Double));
        // The other carboxylate oxygen stays single.
        assert_eq!(residue_bond_order("ASP", "CG", "OD2"), None);
    }

    #[test]
    fn template_pins_a_dipeptide_backbone() {
        // Two glycine residues: N0 CA1 C2 O3 (res 0) - N4 CA5 C6 O7 (res 1), peptide C2-N4.
        let mut top = Topology::default();
        let names = ["N", "CA", "C", "O", "N", "CA", "C", "O"];
        let residues = [0, 0, 0, 0, 1, 1, 1, 1];
        let z = [7u8, 6, 6, 8, 7, 6, 6, 8];
        for k in 0..8 {
            top.atoms.push(
                &Atom::new()
                    .with_atomic_number(z[k])
                    .with_name(names[k])
                    .with_resname("GLY")
                    .with_resindex(residues[k]),
            );
        }
        let bonds = [[0, 1], [1, 2], [2, 3], [2, 4], [4, 5], [5, 6], [6, 7]];
        for &[i, j] in &bonds {
            top.bonds.push(&Bond::new(i, j));
        }
        let pairs: Vec<[usize; 2]> = top.bonds.iter_pairs().collect();
        let orders = residue_template_orders(&top, &pairs);
        // The two backbone carbonyls (C2=O3, C6=O7) are pinned double; everything else, including
        // the peptide bond C2-N4, is left alone.
        assert_eq!(orders[2], Some(BondOrder::Double), "res0 C=O");
        assert_eq!(orders[6], Some(BondOrder::Double), "res1 C=O");
        assert_eq!(orders[3], None, "peptide bond C-N not pinned");
        assert!(orders[0].is_none() && orders[1].is_none());
    }
}
