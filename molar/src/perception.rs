//! Molecular perception from the connection table (no coordinates needed):
//! ring detection (SSSR), aromaticity, and valence / implicit-hydrogen counts.
//!
//! Results live in the **canonical per-element storage**, not parallel arrays:
//! - bond aromaticity → `bond.order == BondOrder::Aromatic` (perception writes it);
//! - per-atom *in-ring* / *aromatic* → the optional `AtomFlags` column
//!   ([`AtomLike::is_in_ring`] / [`AtomLike::is_aromatic`]).
//!
//! The only things with nowhere canonical to live — the ring list and the net charge
//! — are returned in a small [`Perception`] result.
//!
//! Every graph routine here takes a prebuilt [`BondAdjacency`]; none of them builds a
//! throwaway index. For a [`Topology`] that means calling
//! [`BondStorage::ensure_adjacency`] first (the cache then survives the order writes
//! [`perceive`] makes); for a remapped local subgraph — `molar_ff`'s case — it means
//! [`BondAdjacency::build`] over the local pairs.

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

    let z: Vec<u8> = top.atoms.iter().map(|a| a.get_atomic_number()).collect();

    // Build the adjacency into the storage. Discarding the returned reference and re-borrowing
    // below is deliberate: it releases the `&mut` so `ring_is_aromatic` can also read the
    // order column.
    top.bonds.ensure_adjacency(n);

    // The adjacency lives *inside* `top.bonds`, which the annotation pass mutates. `rings` and
    // `aromatic` are owned, so scoping the read here ends the borrow before that.
    let (rings, aromatic) = {
        let adj = top.bonds.get_adjacency().unwrap();
        rings_with_aromaticity(&*top, adj, &z)
    };

    // Now annotate. Every ring atom is in a ring; aromatic rings additionally set the
    // aromatic flag on their atoms and `Aromatic` on their bonds. Writing orders does not
    // invalidate the adjacency, so the cache built above survives this.
    for r in &rings {
        for &a in &r.atoms {
            top.atoms.get_mut(a).unwrap().set_in_ring(true);
        }
    }
    for (r, &is_arom) in rings.iter().zip(&aromatic) {
        if is_arom {
            for &bi in &r.bonds {
                top.bonds.set_order(bi, BondOrder::Aromatic);
            }
            for &a in &r.atoms {
                top.atoms.get_mut(a).unwrap().set_aromatic(true);
            }
        }
    }

    Perception {
        rings: rings.into_iter().map(|r| r.atoms).collect(),
        aromatic,
        total_charge,
    }
}

/// Non-mutating **aromatic**-ring perception: of the SSSR rings, the ones perceived
/// aromatic, each as its atom indices in cycle order.
///
/// The read-only counterpart of [`perceive`], in the same sense [`sssr_rings`] is: it does
/// not aromatize bond orders, so the input's Kekulé structure survives. That matters for
/// callers that need the rings *and* the original orders — force-field typing, charge
/// assignment (espaloma rejects `Aromatic` bonds outright), or a viewer that draws ring
/// circles without editing the molecule it is drawing.
///
/// Takes a prebuilt [`BondAdjacency`] like every other graph routine here.
pub fn aromatic_rings(
    mol: &(impl AtomProvider + BondProvider),
    adj: &BondAdjacency,
) -> Vec<Vec<usize>> {
    let z: Vec<u8> = mol.iter_atoms().map(|a| a.get_atomic_number()).collect();
    let (rings, aromatic) = rings_with_aromaticity(mol, adj, &z);
    rings
        .into_iter()
        .zip(aromatic)
        .filter_map(|(r, is_arom)| is_arom.then_some(r.atoms))
        .collect()
}

/// The SSSR rings plus, parallel to them, whether each is aromatic. Shared by [`perceive`]
/// (which then annotates) and [`aromatic_rings`] (which just filters).
///
/// Aromaticity is decided for **every** ring against the *original* (Kekulé) bond orders
/// before anything is written, so the result never depends on the order rings are processed
/// in — a bond shared by a fused system would otherwise be aromatized before its second
/// ring is tested.
fn rings_with_aromaticity(
    mol: &impl BondProvider,
    adj: &BondAdjacency,
    z: &[u8],
) -> (Vec<RingData>, Vec<bool>) {
    let rings = sssr(adj);

    // Global ring membership (any SSSR ring) — lets the aromaticity test tell a π double
    // bond shared with a *fused* ring (both atoms in rings) from a genuinely exocyclic one
    // (e.g. a carbonyl O outside any ring).
    let mut in_ring = vec![false; adj.n_atoms()];
    for r in &rings {
        for &a in &r.atoms {
            in_ring[a] = true;
        }
    }

    let aromatic = rings.iter().map(|r| ring_is_aromatic(r, mol, adj, z, &in_ring)).collect();
    (rings, aromatic)
}

/// Implicit hydrogens per atom: `round(target_valence − Σ incident bond orders)`,
/// clamped to ≥ 0. `target_valence` is the element's neutral valence adjusted by the
/// atom's formal charge — so a protonated amine N⁺ targets 4.
///
/// **Exact on a Kekulé structure** (Single/Double/Triple bonds). For `Aromatic`-typed
/// bonds (an order-4 SDF, or a structure already run through [`perceive`]) it uses a
/// ring-size heuristic for the aromatic valence (correct for benzene / pyridine /
/// pyrrole / furan / thiophene; approximate for poly-heteroatom azoles like imidazole,
/// where the per-N π role isn't recoverable without the Kekulé form).
///
/// Takes a prebuilt [`BondAdjacency`] like every other graph routine here — it supplies
/// both the atom count and the incident-bond walk, so no throwaway index is built. `mol`
/// must be indexed consistently with it: for a whole [`Topology`], call
/// [`BondStorage::ensure_adjacency`] and pass the topology; for a remapped local subgraph,
/// a [`BondAdjacency::build`] over the local pairs.
///
/// Unspecified-order bonds count as single, so a distance-guessed connectivity table
/// yields the valence-completion counts a sketching/protonation step wants.
pub fn implicit_hydrogens(
    mol: &(impl AtomProvider + BondProvider),
    adj: &BondAdjacency,
) -> Vec<u8> {
    let n = adj.n_atoms();
    let z: Vec<u8> = mol.iter_atoms().map(|a| a.get_atomic_number()).collect();
    let formal_charge: Vec<i32> =
        mol.iter_atoms().map(|a| a.get_formal_charge().unwrap_or(0)).collect();
    let order = |bi: usize| mol.get_bond(bi).map(|b| b.order()).unwrap_or_default();

    // Ring size per atom is only needed to weight aromatic bonds; skip the SSSR work
    // entirely for a plain Kekulé molecule.
    let has_aromatic = mol.iter_bonds().any(|b| b.order() == BondOrder::Aromatic);
    let ring_size = if has_aromatic {
        let mut rs = vec![0usize; n];
        for r in sssr(adj) {
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

    (0..n)
        .map(|i| {
            // Σ over this atom's incident bonds. The valence a bond contributes depends on
            // *this* atom's element/ring (an aromatic pyrrole-N counts 1, a ring carbon 1.5),
            // which is why it's summed per atom rather than once per bond.
            let explicit: f32 = adj
                .neighbors(i)
                .iter()
                .map(|nb| bond_valence(order(nb.bond()), z[i], ring_size[i]))
                .sum();
            let target = target_valence(z[i], formal_charge[i]);
            (target as f32 - explicit).round().max(0.0) as u8
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Valence
// ---------------------------------------------------------------------------

/// Neutral valence of the common (organic) elements; 0 = unhandled (→ no implicit H).
fn base_valence(z: u8) -> i32 {
    match z {
        1 => 1,                // H
        5 => 3,                // B
        6 => 4,                // C
        7 => 3,                // N
        8 => 2,                // O
        9 | 17 | 35 | 53 => 1, // F, Cl, Br, I
        15 => 3,               // P
        16 => 2,               // S
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
        6 => (base - fc.abs()).max(0), // carbocation/carbanion both → 3
        7 | 15 | 8 | 16 => base + fc,  // N⁺→4, O⁺→3, O⁻→1, …
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
///
/// Takes a prebuilt [`BondAdjacency`], which carries its own atom and bond counts; bond order
/// is irrelevant to ring finding, so connectivity is all this needs.
pub fn sssr_rings(adj: &BondAdjacency) -> Vec<Vec<usize>> {
    sssr(adj).into_iter().map(|r| r.atoms).collect()
}

/// One ring: its atom indices (cycle order) and the bond indices forming it.
struct RingData {
    atoms: Vec<usize>,
    bonds: Vec<usize>,
}

fn connected_components(adj: &BondAdjacency) -> usize {
    let n = adj.n_atoms();
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
            for nb in adj.neighbors(x) {
                let y = nb.atom();
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
fn shortest_cycle(adj: &BondAdjacency, u: usize, v: usize, excl: usize) -> Option<RingData> {
    let n = adj.n_atoms();
    let mut prev = vec![usize::MAX; n];
    let mut prev_bond = vec![usize::MAX; n];
    let mut visited = vec![false; n];
    let mut q = VecDeque::from([u]);
    visited[u] = true;
    while let Some(x) = q.pop_front() {
        if x == v {
            break;
        }
        for nb in adj.neighbors(x) {
            let (y, bi) = (nb.atom(), nb.bond());
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
fn sssr(adj: &BondAdjacency) -> Vec<RingData> {
    let n = adj.n_atoms();
    let e = adj.n_bonds();
    if n == 0 || e == 0 {
        return Vec::new();
    }
    let comps = connected_components(adj);
    let mu = (e as isize - n as isize + comps as isize).max(0) as usize; // cyclomatic number
    if mu == 0 {
        return Vec::new();
    }

    // Candidates are generated in ascending bond index — the same order the old bond-slice
    // loop used — because `sort_by_key` below is stable and ties among equal-size rings would
    // otherwise resolve differently, yielding a different (still valid) SSSR set.
    let mut cands: Vec<RingData> = Vec::new();
    for (bi, &[u, v]) in adj.bond_endpoints().iter().enumerate() {
        if u == usize::MAX || v == usize::MAX {
            continue; // self-bond or out-of-range: absent from the adjacency
        }
        if let Some(r) = shortest_cycle(adj, u, v, bi) {
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
    bonds: &impl BondProvider,
    adj: &BondAdjacency,
    z: &[u8],
    in_ring: &[bool],
) -> bool {
    let order = |bi: usize| bonds.get_bond(bi).expect("ring bond index out of range").order();

    let sz = ring.atoms.len();
    if !(5..=6).contains(&sz) {
        return false;
    }
    if ring.bonds.iter().all(|&bi| order(bi) == BondOrder::Aromatic) {
        return true; // already aromatized / SDF order-4
    }

    let mut pi = 0i32;
    for &a in &ring.atoms {
        // A double bond to a ring neighbour (possibly in a fused ring) puts a π electron
        // on `a`; a double bond to a non-ring atom (carbonyl, imine to a substituent) is
        // exocyclic and breaks aromaticity.
        let mut ring_double = false;
        for n in adj.neighbors(a) {
            if order(n.bond()) == BondOrder::Double {
                if in_ring[n.atom()] {
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
            t.atoms.push(&Atom::new().with_atomic_number(n));
        }
        for &(i, j, o) in bonds {
            t.bonds.push(&Bond::with_order(i, j, o));
        }
        t
    }

    /// Every bond order, for asserting what `perceive` wrote into the column.
    fn orders(t: &Topology) -> Vec<BondOrder> {
        t.bonds.iter().map(|b| b.order()).collect()
    }

    use BondOrder::{Double as D, Single as S};

    fn benzene() -> Topology {
        topo(
            &[6, 6, 6, 6, 6, 6],
            &[(0, 1, D), (1, 2, S), (2, 3, D), (3, 4, S), (4, 5, D), (5, 0, S)],
        )
    }

    #[test]
    fn benzene_aromatic() {
        let mut t = benzene();
        let p = perceive(&mut t);
        assert_eq!(p.rings().len(), 1);
        assert_eq!(p.aromatic_rings().count(), 1);
        assert!(orders(&t).iter().all(|&o| o == BondOrder::Aromatic));
        assert!(t.atoms.iter().all(|a| a.is_aromatic() && a.is_in_ring()));
    }

    #[test]
    fn pyridine_aromatic() {
        // atom 0 = N, ring otherwise carbons; Kekulé alternating.
        let mut t = topo(
            &[7, 6, 6, 6, 6, 6],
            &[(0, 1, D), (1, 2, S), (2, 3, D), (3, 4, S), (4, 5, D), (5, 0, S)],
        );
        let p = perceive(&mut t);
        assert_eq!(p.aromatic_rings().count(), 1, "pyridine-type N contributes one π electron");
    }

    #[test]
    fn pyrrole_aromatic() {
        // 5-ring: N(0)-C(1)=C(2)-C(3)=C(4)-N; N single-bonded both sides.
        let mut t = topo(
            &[7, 6, 6, 6, 6],
            &[(0, 1, S), (1, 2, D), (2, 3, S), (3, 4, D), (4, 0, S)],
        );
        let p = perceive(&mut t);
        assert_eq!(p.aromatic_rings().count(), 1, "pyrrole N donates its lone pair");
    }

    #[test]
    fn furan_aromatic() {
        let mut t = topo(
            &[8, 6, 6, 6, 6],
            &[(0, 1, S), (1, 2, D), (2, 3, S), (3, 4, D), (4, 0, S)],
        );
        let p = perceive(&mut t);
        assert_eq!(p.aromatic_rings().count(), 1, "furan O donates its lone pair");
    }

    #[test]
    fn cyclohexane_not_aromatic() {
        let mut t = topo(
            &[6, 6, 6, 6, 6, 6],
            &[(0, 1, S), (1, 2, S), (2, 3, S), (3, 4, S), (4, 5, S), (5, 0, S)],
        );
        let p = perceive(&mut t);
        assert_eq!(p.rings().len(), 1);
        assert_eq!(p.aromatic_rings().count(), 0);
        assert!(orders(&t).iter().all(|&o| o == BondOrder::Single));
        assert!(t.atoms.iter().all(|a| a.is_in_ring() && !a.is_aromatic()));
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
        assert!(orders(&t).iter().all(|&o| o == BondOrder::Aromatic));
    }

    #[test]
    fn biphenyl_link_bond_not_aromatic() {
        // Two benzene rings joined by a single bond (atom 0 – atom 6).
        let bonds = vec![
            (0, 1, D), (1, 2, S), (2, 3, D), (3, 4, S), (4, 5, D), (5, 0, S),
            (6, 7, D), (7, 8, S), (8, 9, D), (9, 10, S), (10, 11, D), (11, 6, S),
            (0, 6, S), // inter-ring link
        ];
        let link = bonds.len() - 1;
        let mut t = topo(&[6; 12], &bonds);
        let p = perceive(&mut t);
        assert_eq!(p.aromatic_rings().count(), 2);
        assert_eq!(orders(&t)[link], BondOrder::Single, "link bond stays single");
    }

    #[test]
    fn formal_charges_sum_into_perception() {
        // CH3–NH3⁺ carries net +1; flipping the N to −1 gives net −1.
        let mut t = topo(&[6, 7], &[(0, 1, S)]);
        t.atoms.get_mut(1).unwrap().set_formal_charge(1);
        assert_eq!(perceive(&mut t).total_charge(), 1.0);
        t.atoms.get_mut(1).unwrap().set_formal_charge(-1);
        assert_eq!(perceive(&mut t).total_charge(), -1.0);
    }

    /// `perceive` writes only into the order column, so the adjacency it built survives — the
    /// whole reason the pair and order columns are stored separately.
    #[test]
    fn perceive_keeps_its_adjacency() {
        let mut t = benzene();
        perceive(&mut t);
        assert!(t.bonds.get_adjacency().is_some());
    }

    #[test]
    fn perceive_preserves_type_id() {
        let mut t = benzene();
        for i in 0..t.atoms.len() {
            t.atoms.get_mut(i).unwrap().set_type_id(42); // a "real" force-field type id
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
        assert!(orders(&t).iter().all(|&o| o == BondOrder::Aromatic));
    }

    /// `aromatic_rings` finds what `perceive` would, without touching the molecule —
    /// the Kekulé orders and the atom flags must all survive.
    #[test]
    fn aromatic_rings_is_non_mutating() {
        let mut t = benzene();
        let before = orders(&t);
        t.bonds.ensure_adjacency(t.atoms.len());
        let adj = t.bonds.get_adjacency().unwrap();
        let rings = aromatic_rings(&t, adj);

        assert_eq!(rings.len(), 1);
        assert_eq!(rings[0].len(), 6);
        assert_eq!(orders(&t), before, "Kekulé orders must survive");
        assert!(t.atoms.iter().all(|a| !a.is_aromatic() && !a.is_in_ring()));

        // ...and it agrees with the mutating pass.
        let perceived: Vec<Vec<usize>> = perceive(&mut t).aromatic_rings().cloned().collect();
        assert_eq!(rings, perceived);
    }

    /// Only aromatic rings come back — an sp3 ring is a ring but not aromatic.
    #[test]
    fn aromatic_rings_skips_aliphatic_rings() {
        let mut t = topo(
            &[6, 6, 6, 6, 6, 6],
            &[(0, 1, S), (1, 2, S), (2, 3, S), (3, 4, S), (4, 5, S), (5, 0, S)],
        );
        t.bonds.ensure_adjacency(t.atoms.len());
        let adj = t.bonds.get_adjacency().unwrap();
        assert!(aromatic_rings(&t, adj).is_empty(), "cyclohexane is not aromatic");
    }

    // -----------------------------------------------------------------------
    // Implicit hydrogens
    // -----------------------------------------------------------------------

    /// `implicit_hydrogens` over a whole topology, building the adjacency it needs.
    fn imp_h(t: &mut Topology) -> Vec<u8> {
        let n = t.atoms.len();
        t.bonds.ensure_adjacency(n);
        let adj = t.bonds.get_adjacency().unwrap();
        implicit_hydrogens(t, adj)
    }

    #[test]
    fn acyclic_implicit_h() {
        // Ethene C=C → 2 H each.
        let mut t = topo(&[6, 6], &[(0, 1, D)]);
        assert_eq!(imp_h(&mut t), vec![2, 2]);
        // Methane (lone C) → 4 H.
        let mut t = topo(&[6], &[]);
        assert_eq!(imp_h(&mut t)[0], 4);
        // Acetylene C#C → 1 H each.
        let mut t = topo(&[6, 6], &[(0, 1, BondOrder::Triple)]);
        assert_eq!(imp_h(&mut t), vec![1, 1]);
    }

    #[test]
    fn unspecified_order_counts_as_single() {
        // A distance-guessed table (PDB/GRO) carries no orders: methyl-ish C with one
        // neighbour still completes to 3 H, which is what a sketching tool wants.
        let mut t = topo(&[6, 6], &[(0, 1, BondOrder::Unspecified)]);
        assert_eq!(imp_h(&mut t), vec![3, 3]);
    }

    #[test]
    fn benzene_aromatic_one_h_each() {
        let mut t = benzene();
        perceive(&mut t);
        for h in imp_h(&mut t) {
            assert_eq!(h, 1, "aromatic ring CH");
        }
    }

    #[test]
    fn pyridine_n_no_h() {
        let mut t = topo(
            &[7, 6, 6, 6, 6, 6],
            &[(0, 1, D), (1, 2, S), (2, 3, D), (3, 4, S), (4, 5, D), (5, 0, S)],
        );
        perceive(&mut t);
        let h = imp_h(&mut t);
        assert_eq!(h[0], 0, "pyridine N has no H");
        assert_eq!(h[1], 1, "ring C has 1 H");
    }

    #[test]
    fn pyrrole_n_one_h() {
        let mut t = topo(&[7, 6, 6, 6, 6], &[(0, 1, S), (1, 2, D), (2, 3, S), (3, 4, D), (4, 0, S)]);
        perceive(&mut t);
        let h = imp_h(&mut t);
        assert_eq!(h[0], 1, "pyrrole N-H");
        assert_eq!(h[1], 1, "ring C-H");
    }

    #[test]
    fn furan_o_no_h() {
        let mut t = topo(&[8, 6, 6, 6, 6], &[(0, 1, S), (1, 2, D), (2, 3, S), (3, 4, D), (4, 0, S)]);
        perceive(&mut t);
        assert_eq!(imp_h(&mut t)[0], 0, "furan O has no H");
    }

    #[test]
    fn cyclohexane_two_h_each() {
        let mut t = topo(
            &[6, 6, 6, 6, 6, 6],
            &[(0, 1, S), (1, 2, S), (2, 3, S), (3, 4, S), (4, 5, S), (5, 0, S)],
        );
        perceive(&mut t);
        for h in imp_h(&mut t) {
            assert_eq!(h, 2, "sp3 ring CH2");
        }
    }

    #[test]
    fn formal_charge_adjusts_valence() {
        // CH3–NH3⁺: a C–N single bond, N carries +1; H are implicit.
        let mut t = topo(&[6, 7], &[(0, 1, S)]);
        t.atoms.get_mut(1).unwrap().set_formal_charge(1);
        let h = imp_h(&mut t);
        assert_eq!(h[0], 3, "methyl C → 3 H");
        assert_eq!(h[1], 3, "ammonium N⁺ (valence 4) → 3 H");

        // A carboxylate-style O⁻ with one single bond → valence 1 → no implicit H.
        let mut t = topo(&[8, 6], &[(0, 1, S)]);
        t.atoms.get_mut(0).unwrap().set_formal_charge(-1);
        assert_eq!(imp_h(&mut t)[0], 0);
    }
}
