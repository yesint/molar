//! Molecular perception from the connection table (no coordinates needed):
//! ring detection (SSSR), aromaticity, and valence / implicit-hydrogen counts.
//!
//! Results live in the **canonical per-element storage**, not parallel arrays:
//! - bond aromaticity → `bond.order == BondOrder::Aromatic` (perception writes it);
//! - per-atom *in-ring* / *aromatic* → the optional `AtomFlags` column
//!   ([`AtomLike::is_in_ring`] / [`AtomLike::is_aromatic`]).
//!
//! The only things with nowhere canonical to live — the ring list and the net charge
//! — are returned in a small [`Perception`] result. Implicit-hydrogen counts are a
//! pure function ([`implicit_hydrogens`]), not stored.

use std::collections::VecDeque;

use crate::prelude::*;

/// Outcome of [`perceive`]: the SSSR ring list and the net molecular charge. (Bond
/// aromaticity is recorded in `bond.order`; per-atom flags in the `AtomFlags` column —
/// see the module docs.)
pub struct Perception {
    rings: Vec<Vec<usize>>,
    /// Parallel to `rings`: whether each ring was perceived aromatic.
    aromatic: Vec<bool>,
    total_charge: Float,
}

impl Perception {
    /// The smallest set of smallest rings; each ring is its atom indices in cycle order.
    pub fn rings(&self) -> &[Vec<usize>] {
        &self.rings
    }

    /// Net molecular charge (Σ per-atom partial charge), summed before perception.
    pub fn total_charge(&self) -> Float {
        self.total_charge
    }

    /// The rings perceived aromatic.
    pub fn aromatic_rings(&self) -> impl Iterator<Item = &Vec<usize>> {
        self.rings
            .iter()
            .zip(&self.aromatic)
            .filter_map(|(r, &a)| a.then_some(r))
    }
}

/// Perceive rings + aromaticity for a topology, annotating it **in place**: sets
/// `BondOrder::Aromatic` on the bonds of every aromatic ring, and the in-ring /
/// aromatic flags on the atoms. Sums the formal charges up front. Returns the
/// [`Perception`] (rings + net charge).
///
/// Destructive of Kekulé structure (an aromatic ring's bonds all become `Aromatic`);
/// idempotent. See also [`System::perceive`](crate::System::perceive).
pub fn perceive(top: &mut Topology) -> Perception {
    let n = top.atoms.len();
    let total_charge: Float = top
        .atoms
        .iter()
        .map(|a| a.get_formal_charge().unwrap_or(0) as Float)
        .sum();

    let rings = sssr(n, &top.bonds);
    let adj = adjacency(n, &top.bonds);
    let z: Vec<u8> = top.atoms.iter().map(|a| a.get_atomic_number()).collect();

    // Global ring membership (any SSSR ring) — lets the aromaticity test tell a π double
    // bond shared with a *fused* ring (both atoms in rings) from a genuinely exocyclic one
    // (e.g. a carbonyl O outside any ring).
    let mut in_ring = vec![false; n];
    for r in &rings {
        for &a in &r.atoms {
            in_ring[a] = true;
        }
    }

    // Decide aromaticity for every ring against the *original* (Kekulé) bond orders
    // first, so the result doesn't depend on the order rings are processed in (a shared
    // bond of a fused system would otherwise be aromatized before its second ring is
    // tested).
    let aromatic: Vec<bool> = rings
        .iter()
        .map(|r| ring_is_aromatic(r, &top.bonds, &adj, &z, &in_ring))
        .collect();

    // Now annotate. Every ring atom is in a ring; aromatic rings additionally set the
    // aromatic flag on their atoms and `Aromatic` on their bonds.
    for r in &rings {
        for &a in &r.atoms {
            top.atoms.get_mut_unchecked(a).set_in_ring(true);
        }
    }
    for (r, &is_arom) in rings.iter().zip(&aromatic) {
        if is_arom {
            for &bi in &r.bonds {
                top.bonds[bi].order = BondOrder::Aromatic;
            }
            for &a in &r.atoms {
                top.atoms.get_mut_unchecked(a).set_aromatic(true);
            }
        }
    }

    Perception {
        rings: rings.into_iter().map(|r| r.atoms).collect(),
        aromatic,
        total_charge,
    }
}

/// Implicit hydrogens per atom: `round(target_valence − Σ incident bond orders)`,
/// clamped to ≥ 0. `target_valence` is the element's neutral valence adjusted by the
/// atom's formal charge (`Atom.formal_charge`) — so a protonated amine N⁺ targets 4.
///
/// **Exact on a Kekulé structure** (Single/Double/Triple bonds). For `Aromatic`-typed
/// bonds (an order-4 SDF, or a structure already run through [`perceive`]) it uses a
/// ring-size heuristic for the aromatic valence (correct for benzene / pyridine /
/// pyrrole / furan / thiophene; approximate for poly-heteroatom azoles like imidazole,
/// where the per-N π role isn't recoverable without the Kekulé form).
pub fn implicit_hydrogens(sel: &(impl AtomProvider + BondProvider + LenProvider)) -> Vec<u8> {
    let n = sel.len();
    let bonds: Vec<Bond> = sel.iter_bonds().copied().collect();
    let z: Vec<u8> = sel.iter_atoms().map(|a| a.get_atomic_number()).collect();
    let formal_charge: Vec<i32> =
        sel.iter_atoms().map(|a| a.get_formal_charge().unwrap_or(0)).collect();

    // Ring size per atom is only needed to weight aromatic bonds; skip the SSSR work
    // entirely for a plain Kekulé molecule.
    let has_aromatic = bonds.iter().any(|b| b.order == BondOrder::Aromatic);
    let ring_size = if has_aromatic {
        let mut rs = vec![0usize; n];
        for r in sssr(n, &bonds) {
            let sz = r.atoms.len();
            for a in r.atoms {
                if rs[a] == 0 || sz < rs[a] {
                    rs[a] = sz;
                }
            }
        }
        rs
    } else {
        vec![0; n]
    };

    let mut explicit = vec![0.0f32; n];
    for b in &bonds {
        if b.i1 >= n || b.i2 >= n {
            continue;
        }
        explicit[b.i1] += bond_valence(b.order, z[b.i1], ring_size[b.i1]);
        explicit[b.i2] += bond_valence(b.order, z[b.i2], ring_size[b.i2]);
    }

    (0..n)
        .map(|i| {
            let target = target_valence(z[i], formal_charge[i]);
            let h = (target as f32 - explicit[i]).round();
            h.max(0.0) as u8
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Valence
// ---------------------------------------------------------------------------

/// Neutral valence of the common (organic) elements; 0 = unhandled (→ no implicit H).
fn base_valence(z: u8) -> i32 {
    match z {
        1 => 1,                       // H
        5 => 3,                       // B
        6 => 4,                       // C
        7 => 3,                       // N
        8 => 2,                       // O
        9 | 17 | 35 | 53 => 1,        // F, Cl, Br, I
        15 => 3,                      // P
        16 => 2,                      // S
        _ => 0,
    }
}

/// Target valence adjusted for the atom's formal charge `fc`.
fn target_valence(z: u8, fc: i32) -> i32 {
    let base = base_valence(z);
    if base == 0 {
        return 0;
    }
    match z {
        6 => (base - fc.abs()).max(0),      // carbocation/carbanion both → 3
        7 | 15 | 8 | 16 => base + fc,       // N⁺→4, O⁺→3, O⁻→1, …
        _ => (base + fc).max(0),
    }
}

/// Valence contributed by one incident bond, given the bonded atom's element + the size
/// of the smallest ring it sits in (used only for `Aromatic` bonds).
fn bond_valence(order: BondOrder, z: u8, ring_size: usize) -> f32 {
    match order {
        BondOrder::Single | BondOrder::Unspecified => 1.0,
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        // pyrrole-N / furan-O / thiophene-S are σ-bonded lone-pair donors (the aromatic
        // bond is really order 1 for valence); everything else (C, pyridine-type N in a
        // 6-ring) averages to 1.5.
        BondOrder::Aromatic => match z {
            7 if ring_size == 5 => 1.0,
            8 | 16 => 1.0,
            _ => 1.5,
        },
    }
}

// ---------------------------------------------------------------------------
// Graph + SSSR
// ---------------------------------------------------------------------------

/// Non-mutating **SSSR** ring perception: the smallest set of smallest rings, each as its
/// atom indices in cycle order. Unlike [`perceive`] this does not touch the topology (no
/// Kekulé destruction), so callers that need the original bond orders — e.g. force-field
/// featurization — can get the ring set without side effects.
pub fn sssr_rings(n: usize, bonds: &[Bond]) -> Vec<Vec<usize>> {
    sssr(n, bonds).into_iter().map(|r| r.atoms).collect()
}

/// One ring: its atom indices (cycle order) and the bond indices forming it.
struct RingData {
    atoms: Vec<usize>,
    bonds: Vec<usize>,
}

/// Per-atom neighbour list of `(neighbor_atom, bond_index)`, skipping self/out-of-range bonds.
fn adjacency(n: usize, bonds: &[Bond]) -> Vec<Vec<(usize, usize)>> {
    let mut adj = vec![Vec::new(); n];
    for (bi, b) in bonds.iter().enumerate() {
        if b.i1 < n && b.i2 < n && b.i1 != b.i2 {
            adj[b.i1].push((b.i2, bi));
            adj[b.i2].push((b.i1, bi));
        }
    }
    adj
}

fn connected_components(n: usize, adj: &[Vec<(usize, usize)>]) -> usize {
    let mut seen = vec![false; n];
    let mut count = 0;
    for s in 0..n {
        if seen[s] {
            continue;
        }
        count += 1;
        let mut q = VecDeque::from([s]);
        seen[s] = true;
        while let Some(x) = q.pop_front() {
            for &(y, _) in &adj[x] {
                if !seen[y] {
                    seen[y] = true;
                    q.push_back(y);
                }
            }
        }
    }
    count
}

/// Smallest ring through bond `(u,v)` (the bond `excl` is the closing edge): BFS the
/// shortest `u→v` path that doesn't use `excl`.
fn shortest_cycle(
    adj: &[Vec<(usize, usize)>],
    n: usize,
    u: usize,
    v: usize,
    excl: usize,
) -> Option<RingData> {
    let mut prev = vec![usize::MAX; n];
    let mut prev_bond = vec![usize::MAX; n];
    let mut visited = vec![false; n];
    let mut q = VecDeque::from([u]);
    visited[u] = true;
    while let Some(x) = q.pop_front() {
        if x == v {
            break;
        }
        for &(y, bi) in &adj[x] {
            if bi == excl || visited[y] {
                continue;
            }
            visited[y] = true;
            prev[y] = x;
            prev_bond[y] = bi;
            q.push_back(y);
        }
    }
    if !visited[v] {
        return None;
    }
    let mut atoms = Vec::new();
    let mut bonds = vec![excl];
    let mut cur = v;
    while cur != u {
        atoms.push(cur);
        bonds.push(prev_bond[cur]);
        cur = prev[cur];
        if cur == usize::MAX {
            return None; // shouldn't happen once `v` is reached
        }
    }
    atoms.push(u);
    atoms.reverse();
    Some(RingData { atoms, bonds })
}

/// Smallest set of smallest rings, via "smallest ring per bond" + GF(2) independence.
fn sssr(n: usize, bonds: &[Bond]) -> Vec<RingData> {
    if n == 0 || bonds.is_empty() {
        return Vec::new();
    }
    let adj = adjacency(n, bonds);
    let e = bonds.len();
    let comps = connected_components(n, &adj);
    let mu = (e as isize - n as isize + comps as isize).max(0) as usize; // cyclomatic number
    if mu == 0 {
        return Vec::new();
    }

    let mut cands: Vec<RingData> = Vec::new();
    for (bi, b) in bonds.iter().enumerate() {
        if b.i1 >= n || b.i2 >= n || b.i1 == b.i2 {
            continue;
        }
        if let Some(r) = shortest_cycle(&adj, n, b.i1, b.i2, bi) {
            cands.push(r);
        }
    }
    cands.sort_by_key(|r| r.bonds.len());

    // Keep candidates that are linearly independent over GF(2) on the edge set, until we
    // have `mu` of them. Each ring is a bitvector over bond indices; Gaussian elimination.
    let words = e.div_ceil(64);
    let mut basis: Vec<(usize, Vec<u64>)> = Vec::new(); // (pivot bit, reduced row)
    let mut chosen = Vec::new();
    for cand in cands {
        if chosen.len() == mu {
            break;
        }
        let mut bits = vec![0u64; words];
        for &bi in &cand.bonds {
            bits[bi / 64] |= 1u64 << (bi % 64);
        }
        for (piv, row) in &basis {
            if bits[piv / 64] & (1u64 << (piv % 64)) != 0 {
                for (d, s) in bits.iter_mut().zip(row) {
                    *d ^= *s;
                }
            }
        }
        if let Some(piv) = lowest_set_bit(&bits) {
            basis.push((piv, bits));
            chosen.push(cand);
        }
    }
    chosen
}

fn lowest_set_bit(v: &[u64]) -> Option<usize> {
    for (wi, &w) in v.iter().enumerate() {
        if w != 0 {
            return Some(wi * 64 + w.trailing_zeros() as usize);
        }
    }
    None
}

// ---------------------------------------------------------------------------
// Aromaticity (Hückel)
// ---------------------------------------------------------------------------

/// Whether a ring is aromatic: all ring bonds already `Aromatic` (trust the input), or a
/// Hückel 4n+2 π-count over sp2 ring atoms. v1 handles 5- and 6-membered rings; an
/// exocyclic double bond (e.g. carbonyl) or an sp3 ring atom breaks aromaticity.
fn ring_is_aromatic(
    ring: &RingData,
    bonds: &[Bond],
    adj: &[Vec<(usize, usize)>],
    z: &[u8],
    in_ring: &[bool],
) -> bool {
    let sz = ring.atoms.len();
    if !(5..=6).contains(&sz) {
        return false;
    }
    if ring
        .bonds
        .iter()
        .all(|&bi| bonds[bi].order == BondOrder::Aromatic)
    {
        return true; // already aromatized / SDF order-4
    }

    let mut pi = 0i32;
    for &a in &ring.atoms {
        // A double bond to a ring neighbour (possibly in a fused ring) puts a π electron
        // on `a`; a double bond to a non-ring atom (carbonyl, imine to a substituent) is
        // exocyclic and breaks aromaticity.
        let mut ring_double = false;
        for &(nb, bi) in &adj[a] {
            if bonds[bi].order == BondOrder::Double {
                if in_ring[nb] {
                    ring_double = true;
                } else {
                    return false;
                }
            }
        }
        match z[a] {
            6 => {
                if ring_double {
                    pi += 1;
                } else {
                    return false; // sp3 carbon
                }
            }
            7 => pi += if ring_double { 1 } else { 2 }, // pyridine vs pyrrole
            8 | 16 => {
                if ring_double {
                    return false;
                } else {
                    pi += 2; // furan-O / thiophene-S lone pair
                }
            }
            _ => return false,
        }
    }
    matches!(pi, 2 | 6 | 10)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a topology from atomic numbers + (i, j, order) bonds.
    fn topo(z: &[u8], bonds: &[(usize, usize, BondOrder)]) -> Topology {
        let mut t = Topology::default();
        for &n in z {
            t.atoms.push_row(&Atom::new().with_atomic_number(n));
        }
        for &(i, j, o) in bonds {
            t.bonds.push(Bond::with_order(i, j, o));
        }
        t
    }

    use BondOrder::{Double as D, Single as S};

    fn benzene() -> Topology {
        topo(
            &[6, 6, 6, 6, 6, 6],
            &[(0, 1, D), (1, 2, S), (2, 3, D), (3, 4, S), (4, 5, D), (5, 0, S)],
        )
    }

    #[test]
    fn benzene_aromatic_one_h_each() {
        let mut t = benzene();
        let p = perceive(&mut t);
        assert_eq!(p.rings().len(), 1);
        assert_eq!(p.aromatic_rings().count(), 1);
        assert!(t.bonds.iter().all(|b| b.order == BondOrder::Aromatic));
        assert!(t.atoms.iter().all(|a| a.is_aromatic() && a.is_in_ring()));
        for h in implicit_hydrogens(&t) {
            assert_eq!(h, 1);
        }
    }

    #[test]
    fn pyridine_n_no_h() {
        // atom 0 = N, ring otherwise carbons; Kekulé alternating.
        let mut t = topo(
            &[7, 6, 6, 6, 6, 6],
            &[(0, 1, D), (1, 2, S), (2, 3, D), (3, 4, S), (4, 5, D), (5, 0, S)],
        );
        perceive(&mut t);
        let h = implicit_hydrogens(&t);
        assert_eq!(h[0], 0, "pyridine N has no H");
        assert_eq!(h[1], 1, "ring C has 1 H");
    }

    #[test]
    fn pyrrole_n_one_h() {
        // 5-ring: N(0)-C(1)=C(2)-C(3)=C(4)-N; N single-bonded both sides.
        let mut t = topo(
            &[7, 6, 6, 6, 6],
            &[(0, 1, S), (1, 2, D), (2, 3, S), (3, 4, D), (4, 0, S)],
        );
        let p = perceive(&mut t);
        assert_eq!(p.aromatic_rings().count(), 1);
        let h = implicit_hydrogens(&t);
        assert_eq!(h[0], 1, "pyrrole N-H");
        assert_eq!(h[1], 1, "ring C-H");
    }

    #[test]
    fn furan_o_no_h() {
        let mut t = topo(
            &[8, 6, 6, 6, 6],
            &[(0, 1, S), (1, 2, D), (2, 3, S), (3, 4, D), (4, 0, S)],
        );
        let p = perceive(&mut t);
        assert_eq!(p.aromatic_rings().count(), 1);
        assert_eq!(implicit_hydrogens(&t)[0], 0, "furan O has no H");
    }

    #[test]
    fn cyclohexane_not_aromatic_two_h() {
        let mut t = topo(
            &[6, 6, 6, 6, 6, 6],
            &[(0, 1, S), (1, 2, S), (2, 3, S), (3, 4, S), (4, 5, S), (5, 0, S)],
        );
        let p = perceive(&mut t);
        assert_eq!(p.rings().len(), 1);
        assert_eq!(p.aromatic_rings().count(), 0);
        assert!(t.bonds.iter().all(|b| b.order == BondOrder::Single));
        assert!(t.atoms.iter().all(|a| a.is_in_ring() && !a.is_aromatic()));
        for h in implicit_hydrogens(&t) {
            assert_eq!(h, 2);
        }
    }

    #[test]
    fn cyclohexanone_not_aromatic() {
        // 6 ring carbons (all single) + exocyclic O double-bonded to C0.
        let mut t = topo(
            &[6, 6, 6, 6, 6, 6, 8],
            &[
                (0, 1, S), (1, 2, S), (2, 3, S), (3, 4, S), (4, 5, S), (5, 0, S),
                (0, 6, D),
            ],
        );
        let p = perceive(&mut t);
        assert_eq!(p.aromatic_rings().count(), 0, "carbonyl breaks aromaticity");
    }

    #[test]
    fn naphthalene_two_aromatic_rings() {
        // Two fused 6-rings sharing the 0–1 bond. A valid Kekulé form.
        let mut t = topo(
            &[6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            &[
                (0, 1, S),               // shared bond
                (1, 2, D), (2, 3, S), (3, 4, D), (4, 5, S), (5, 0, D), // ring A
                (1, 6, S), (6, 7, D), (7, 8, S), (8, 9, D), (9, 0, S), // ring B
            ],
        );
        let p = perceive(&mut t);
        assert_eq!(p.rings().len(), 2);
        assert_eq!(p.aromatic_rings().count(), 2);
        assert!(t.bonds.iter().all(|b| b.order == BondOrder::Aromatic));
    }

    #[test]
    fn biphenyl_link_bond_not_aromatic() {
        // Two benzene rings joined by a single bond (atom 0 – atom 6).
        let mut bonds = vec![
            (0, 1, D), (1, 2, S), (2, 3, D), (3, 4, S), (4, 5, D), (5, 0, S),
            (6, 7, D), (7, 8, S), (8, 9, D), (9, 10, S), (10, 11, D), (11, 6, S),
            (0, 6, S), // inter-ring link
        ];
        let link = bonds.len() - 1;
        let mut t = topo(&[6; 12], &bonds);
        let p = perceive(&mut t);
        assert_eq!(p.aromatic_rings().count(), 2);
        assert_eq!(t.bonds[link].order, BondOrder::Single, "link bond stays single");
        let _ = &mut bonds;
    }

    #[test]
    fn ammonium_charge_adjusts_valence() {
        // CH3–NH3⁺: a C–N single bond, N carries +1; H are implicit.
        let mut t = topo(&[6, 7], &[(0, 1, S)]);
        t.atoms.get_mut_unchecked(1).set_formal_charge(1);
        let p = perceive(&mut t);
        assert_eq!(p.total_charge(), 1.0);
        let h = implicit_hydrogens(&t);
        assert_eq!(h[0], 3, "methyl C → 3 H");
        assert_eq!(h[1], 3, "ammonium N⁺ (valence 4) → 3 H");
    }

    #[test]
    fn carboxylate_oxygen_minus() {
        // A lone O⁻ with one single bond → valence 1 → 0 implicit H.
        let mut t = topo(&[8, 6], &[(0, 1, S)]);
        t.atoms.get_mut_unchecked(0).set_formal_charge(-1);
        let _ = perceive(&mut t);
        assert_eq!(implicit_hydrogens(&t)[0], 0);
    }

    #[test]
    fn acyclic_implicit_h() {
        // Ethene C=C → 2 H each.
        let t = topo(&[6, 6], &[(0, 1, D)]);
        for h in implicit_hydrogens(&t) {
            assert_eq!(h, 2);
        }
        // Methane (lone C) → 4 H.
        let t = topo(&[6], &[]);
        assert_eq!(implicit_hydrogens(&t)[0], 4);
    }

    #[test]
    fn perceive_preserves_type_id() {
        let mut t = benzene();
        for i in 0..t.atoms.len() {
            t.atoms.get_mut_unchecked(i).set_type_id(42); // a "real" force-field type id
        }
        perceive(&mut t);
        for a in t.atoms.iter() {
            assert!(a.is_aromatic() && a.is_in_ring());
            assert_eq!(a.get_type_id(), Some(42), "perceive leaves type_id untouched");
        }
    }

    #[test]
    fn sdf_order4_ring_stays_aromatic() {
        // A 6-ring whose bonds are already Aromatic (as an SDF order-4 record yields).
        let mut t = topo(
            &[6, 6, 6, 6, 6, 6],
            &[
                (0, 1, BondOrder::Aromatic), (1, 2, BondOrder::Aromatic),
                (2, 3, BondOrder::Aromatic), (3, 4, BondOrder::Aromatic),
                (4, 5, BondOrder::Aromatic), (5, 0, BondOrder::Aromatic),
            ],
        );
        let p = perceive(&mut t);
        assert_eq!(p.aromatic_rings().count(), 1);
        for h in implicit_hydrogens(&t) {
            assert_eq!(h, 1);
        }
    }
}
