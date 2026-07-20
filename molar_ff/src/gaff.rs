//! GAFF atom-type perception.
//!
//! The typing core operates purely on a LOCAL (0..n) index space: an array of atomic
//! numbers plus a list of [`LocalBond`]s. It is deliberately free of MolAR selection
//! types so it can be unit-tested directly. [`crate::ApplyFF`] is responsible for
//! extracting this local view from a `System`/selection and writing results back.
//!
//! molar convention: `gaff.rs` + a `gaff/` directory, never `mod.rs`.

use crate::FFError;

mod aromatic;
mod conjugation;
mod matcher;
mod props;
mod ring;
mod tables;

/// A bond in the local (0..n) index space, with an integer bond-type code
/// (1 single, 2 double, 3 triple, 4 aromatic).
#[derive(Debug, Clone, Copy)]
pub struct LocalBond {
    pub i: usize,
    pub j: usize,
    pub order: u8,
}

/// Per-atom perceived context (rings, ring-size counts, aromaticity, EW flag).
pub(crate) struct Perceived {
    /// Neighbour lists in input-bond order.
    pub con: Vec<Vec<usize>>,
    /// The perceived ring set (each ring = sorted member atoms).
    #[allow(dead_code)]
    pub rings: Vec<Vec<usize>>,
    /// Per-atom ring-size counts: `rg[a][0]` total, `rg[a][k]` = #k-membered rings.
    pub rg: Vec<[u16; 11]>,
    /// AR1..AR5 counts, `ewd`, `nr`.
    pub arom: aromatic::Arom,
}

/// Build neighbour lists (`con`) from the local bond list, in bond order.
pub(crate) fn build_con(n: usize, bonds: &[LocalBond]) -> Vec<Vec<usize>> {
    let mut con = vec![Vec::new(); n];
    for b in bonds {
        con[b.i].push(b.j);
        con[b.j].push(b.i);
    }
    con
}

/// Run the coordinate-free perception pipeline: ring detection → ring properties →
/// aromaticity/EW classification.
pub(crate) fn perceive(z: &[u8], bonds: &[LocalBond]) -> Perceived {
    let n = z.len();
    let con = build_con(n, bonds);
    let rings = ring::detect_rings(z, &con);
    let rg = ring::ring_property(n, &rings);
    let arom = aromatic::aromatic(z, &con, bonds, &rings, &rg);
    Perceived { con, rings, rg, arom }
}

/// Assign a GAFF atom type to every atom, given atomic numbers `z` and the local bond
/// list. Returns one owned type name per atom, in the same order as `z`.
///
/// This runs the full pipeline: perception (rings/aromaticity/EW) → per-atom property
/// precompute → the rule-matching loop over the embedded DEF → the conjugation parity
/// split.
pub fn gaff_types(z: &[u8], bonds: &[LocalBond]) -> Result<Vec<String>, FFError> {
    let p = perceive(z, bonds);
    let pr = props::compute(z, &p.con, bonds, &p.arom.ewd);
    let ctx = matcher::Ctx::new(
        z, &p.con, bonds, &pr, &p.rg, &p.arom.ar, &p.arom.nr, &p.arom.ewd,
    );
    let mut types = matcher::jat(&ctx)?;
    conjugation::adjust(&mut types, bonds);
    Ok(types)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: build (z, bonds) and perceive.
    fn perceive_mol(z: &[u8], bonds: &[(usize, usize, u8)]) -> Perceived {
        let b: Vec<LocalBond> = bonds.iter().map(|&(i, j, o)| LocalBond { i, j, order: o }).collect();
        perceive(z, &b)
    }

    #[test]
    fn benzene_all_ar1() {
        // 0..5 ring C, 6..11 H on each C.
        let z = [6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1];
        let bonds = [
            (0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1),
            (0, 6, 1), (1, 7, 1), (2, 8, 1), (3, 9, 1), (4, 10, 1), (5, 11, 1),
        ];
        let p = perceive_mol(&z, &bonds);
        assert_eq!(p.rings.len(), 1);
        assert_eq!(p.rings[0].len(), 6);
        for c in 0..6 {
            assert_eq!(p.arom.ar[c][1], 1, "ring C{c} must be AR1");
            assert_eq!(p.rg[c][6], 1, "ring C{c} must be in one 6-ring");
            assert!(!p.arom.nr[c]);
        }
        for h in 6..12 {
            assert!(p.arom.nr[h], "H{h} must be non-ring");
        }
    }

    #[test]
    fn pyridine_ar1_with_n() {
        // 0 = N (connum 2, no H), 1..5 = C ring, H on each C (6..10).
        let z = [7, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1];
        let bonds = [
            (0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1),
            (1, 6, 1), (2, 7, 1), (3, 8, 1), (4, 9, 1), (5, 10, 1),
        ];
        let p = perceive_mol(&z, &bonds);
        assert_eq!(p.rings.len(), 1);
        for a in 0..6 {
            assert_eq!(p.arom.ar[a][1], 1, "atom {a} must be AR1");
        }
    }

    #[test]
    fn cyclohexane_all_ar5() {
        // 6 sp3 C (each connum 4: 2 ring C + 2 H), 12 H.
        let mut z = vec![6u8; 6];
        z.extend(std::iter::repeat(1u8).take(12));
        let mut bonds = vec![
            (0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1),
        ];
        // two H per ring carbon
        for c in 0..6 {
            bonds.push((c, 6 + 2 * c, 1));
            bonds.push((c, 7 + 2 * c, 1));
        }
        let p = perceive_mol(&z, &bonds);
        assert_eq!(p.rings.len(), 1);
        for c in 0..6 {
            assert_eq!(p.arom.ar[c][5], 1, "ring C{c} must be AR5 (aliphatic)");
            assert_eq!(p.arom.ar[c][1], 0);
        }
    }

    #[test]
    fn imidazole_is_ar2() {
        // c1cnc[nH]1 : ring atoms 0=C,1=C,2=N(pyridine,connum2),3=C,4=N(H,connum3)
        // H: 5,6 on C0,C1 (aromatic CH); 7 on C3; 8 on N4.
        let z = [6, 6, 7, 6, 7, 1, 1, 1, 1];
        let bonds = [
            (0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 0, 1),
            (0, 5, 1), (1, 6, 1), (3, 7, 1), (4, 8, 1),
        ];
        let p = perceive_mol(&z, &bonds);
        assert_eq!(p.rings.len(), 1);
        assert_eq!(p.rings[0].len(), 5);
        for a in 0..5 {
            assert_eq!(p.arom.ar[a][2], 1, "ring atom {a} must be AR2");
            assert_eq!(p.arom.ar[a][1], 0, "5-ring atom {a} must not be AR1");
        }
    }

    #[test]
    fn naphthalene_two_six_rings() {
        // 10 C fused bicyclic; fusion atoms 4 and 9 (connum 3, no H); others have 1 H.
        let mut z = vec![6u8; 10];
        z.extend(std::iter::repeat(1u8).take(8));
        let mut bonds = vec![
            (0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 1),
            (5, 6, 2), (6, 7, 1), (7, 8, 2), (8, 9, 1), (9, 0, 1), (4, 9, 2),
        ];
        // H on peripheral carbons 0,1,2,3,5,6,7,8
        let periph = [0, 1, 2, 3, 5, 6, 7, 8];
        for (k, &c) in periph.iter().enumerate() {
            bonds.push((c, 10 + k, 1));
        }
        let p = perceive_mol(&z, &bonds);
        assert_eq!(p.rings.len(), 2, "naphthalene must yield two 6-rings (no perimeter)");
        for r in &p.rings {
            assert_eq!(r.len(), 6);
        }
        // fusion atoms are in two rings; peripheral in one
        assert_eq!(p.rg[4][6], 2);
        assert_eq!(p.rg[9][6], 2);
        assert_eq!(p.rg[0][6], 1);
        // all ring carbons are AR1 (pure aromatic)
        for c in 0..10 {
            assert!(p.arom.ar[c][1] >= 1, "ring C{c} must be AR1");
        }
    }
}
