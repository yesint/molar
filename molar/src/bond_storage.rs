//! Column-wise storage for bonds, plus the derived bonded-adjacency index.
//!
//! [`BondStorage`] replaces `Vec<Bond>` inside [`Topology`](crate::Topology). It holds:
//!
//! - the **pair column**, always present, one `[u32; 2]` per bond;
//! - the **order column**, `Option<Vec<BondOrder>>` — `None` costs nothing and stays `None` for
//!   sources that record no orders (PDB CONECT, TPR, GRO); when `Some` it is full-length;
//! - a cached [`BondAdjacency`], invalidated by structural change but **surviving order
//!   writes** — which is what lets [`perceive`](crate::perceive) aromatize rings
//!   without throwing its adjacency away.
//!
//! Bonds are read through the borrowed [`BondRef`] proxy; the owned [`Bond`] row remains the
//! construction/interchange type.
//!
//! Index width is `u32` **internally only** — every public boundary speaks `usize`. This caps a
//! system at `u32::MAX` atoms, which [`to_idx`] enforces with a panic rather than a silent
//! truncation.

use crate::bond::{Bond, BondOrder};

/// Internal index width for the bond columns. Never appears in a public signature.
type BondIdx = u32;

/// Narrow an atom or bond index into the storage width, panicking rather than truncating.
#[inline]
fn to_idx(i: usize) -> BondIdx {
    BondIdx::try_from(i).expect("index exceeds the bond-storage index width (u32)")
}

/// Column-wise storage for bonds. See the [module docs](self).
#[derive(Debug, Default)]
pub struct BondStorage {
    /// Always present, `len() == n_bonds`.
    pairs: Vec<[BondIdx; 2]>,
    /// `None` => no source recorded an order; `Some` => full-length.
    orders: Option<Vec<BondOrder>>,
    /// Derived index. A plain field, deliberately not a `OnceCell`: it must be invalidatable,
    /// and interior mutability would cost `Topology` its `Sync`-through-`&` sharing across rayon.
    adj: Option<BondAdjacency>,
}

/// A derived cache is never worth copying — a clone rebuilds it on demand.
impl Clone for BondStorage {
    fn clone(&self) -> Self {
        Self { pairs: self.pairs.clone(), orders: self.orders.clone(), adj: None }
    }
}

/// Push into the order column, materializing it (backfilling `Unspecified`) on the first bond
/// that carries a concrete order. `new_len` is the pair-column length *after* the push.
#[inline]
fn push_order(col: &mut Option<Vec<BondOrder>>, val: BondOrder, new_len: usize) {
    match (col.as_mut(), val) {
        (Some(v), x) => v.push(x),
        // A bond with no order carries no information, so it does not force the column.
        (None, BondOrder::Unspecified) => {}
        (None, x) => {
            let mut v = vec![BondOrder::Unspecified; new_len - 1];
            v.push(x);
            *col = Some(v);
        }
    }
}

impl BondStorage {
    /// Number of bonds.
    pub fn len(&self) -> usize {
        self.pairs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    /// Whether the order column is allocated. `false` for connectivity-only sources, in which
    /// case every [`BondRef::order`] reads `Unspecified`.
    pub fn has_orders(&self) -> bool {
        self.orders.is_some()
    }

    /// Append a bond, scattering it into the columns. The order column materializes on the
    /// first bond carrying a concrete order. Invalidates the adjacency.
    pub fn push(&mut self, b: &Bond) {
        self.pairs.push([to_idx(b.i1), to_idx(b.i2)]);
        let n = self.pairs.len();
        push_order(&mut self.orders, b.order, n);
        self.adj = None;
        debug_assert!(self.invariant_holds());
    }

    /// Reserve capacity for `n` more bonds in every present column.
    pub fn reserve(&mut self, n: usize) {
        self.pairs.reserve(n);
        if let Some(c) = self.orders.as_mut() {
            c.reserve(n);
        }
    }

    /// A read-only proxy for bond `i`.
    ///
    /// # Safety
    /// `i` must be `< len()`.
    #[inline]
    pub unsafe fn get_unchecked(&self, i: usize) -> BondRef<'_> {
        BondRef { st: self, idx: i }
    }

    /// A read-only proxy for bond `i`, or `None` if out of bounds.
    #[inline]
    pub fn get(&self, i: usize) -> Option<BondRef<'_>> {
        (i < self.len()).then_some(BondRef { st: self, idx: i })
    }

    /// Iterate read-only proxies over all bonds.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = BondRef<'_>> {
        // SAFETY: `i` ranges over `0..len()`.
        (0..self.len()).map(move |i| unsafe { self.get_unchecked(i) })
    }

    /// Set the order of bond `i`, materializing the column if absent.
    ///
    /// Deliberately does **not** invalidate the adjacency — orders don't affect connectivity.
    pub fn set_order(&mut self, i: usize, o: BondOrder) {
        let n = self.pairs.len();
        assert!(i < n, "bond index {i} out of bounds (n_bonds = {n})");
        match (self.orders.as_mut(), o) {
            (Some(v), x) => v[i] = x,
            (None, BondOrder::Unspecified) => {}
            (None, x) => {
                let mut v = vec![BondOrder::Unspecified; n];
                v[i] = x;
                self.orders = Some(v);
            }
        }
        debug_assert!(self.invariant_holds());
    }

    /// Drop every bond with a removed or out-of-range endpoint and renumber the survivors into
    /// the post-removal atom index space. `removed_atoms` and `n_atoms` are the *pre-removal*
    /// atom indices to drop and the pre-removal atom count — i.e. exactly the inputs
    /// `AtomStorage::retain_by_index` is given.
    pub fn remove_by_index(&mut self, removed_atoms: &[usize], n_atoms: usize) {
        let mut keep = vec![true; n_atoms];
        for &i in removed_atoms {
            if i < n_atoms {
                keep[i] = false;
            }
        }
        // Old atom index -> new atom index (only meaningful where `keep`).
        let mut remap = vec![0 as BondIdx; n_atoms];
        let mut next: BondIdx = 0;
        for i in 0..n_atoms {
            if keep[i] {
                remap[i] = next;
                next += 1;
            }
        }

        // Compact in place: `w` trails `r`, so the order-column copy never reads a stale slot.
        let mut orders = self.orders.take();
        let mut w = 0usize;
        for r in 0..self.pairs.len() {
            let [a, b] = self.pairs[r];
            let (a, b) = (a as usize, b as usize);
            if a >= n_atoms || b >= n_atoms || !keep[a] || !keep[b] {
                continue;
            }
            self.pairs[w] = [remap[a], remap[b]];
            if let Some(v) = orders.as_mut() {
                v[w] = v[r];
            }
            w += 1;
        }
        self.pairs.truncate(w);
        if let Some(v) = orders.as_mut() {
            v.truncate(w);
        }
        self.orders = orders;
        self.adj = None;
        debug_assert!(self.invariant_holds());
    }

    /// The cached bonded adjacency, or `None` if it was never built or has been invalidated.
    ///
    /// Cheap and parallel-safe: no locking, no interior mutability. Build it with
    /// [`ensure_adjacency`](Self::ensure_adjacency) *before* entering a parallel region — the
    /// same rule the optional atom columns follow.
    pub fn get_adjacency(&self) -> Option<&BondAdjacency> {
        self.adj.as_ref()
    }

    /// Build the adjacency if absent or built for a different atom count, then return it.
    pub fn ensure_adjacency(&mut self, n_atoms: usize) -> &BondAdjacency {
        let stale = self.adj.as_ref().is_none_or(|a| a.n_atoms() != n_atoms);
        if stale {
            self.adj = Some(BondAdjacency::build(n_atoms, self.iter_pairs()));
        }
        self.adj.as_ref().unwrap()
    }

    /// Discard the cached adjacency. Callers that change the **atom count** must invoke this:
    /// `offsets` is sized `n_atoms + 1`, so a stale index would be short.
    pub fn invalidate_adjacency(&mut self) {
        self.adj = None;
    }

    /// Endpoint pairs in bond-index order. Feeds [`BondAdjacency::build`] and the SSSR
    /// candidate loop, both of which need bonds in their original index order.
    pub fn iter_pairs(&self) -> impl ExactSizeIterator<Item = [usize; 2]> + Clone + '_ {
        self.pairs.iter().map(|&[a, b]| [a as usize, b as usize])
    }

    /// Debug check: a present order column matches the pair-column length.
    fn invariant_holds(&self) -> bool {
        self.orders.as_ref().is_none_or(|c| c.len() == self.pairs.len())
    }
}

impl Extend<Bond> for BondStorage {
    fn extend<I: IntoIterator<Item = Bond>>(&mut self, iter: I) {
        let it = iter.into_iter();
        self.reserve(it.size_hint().0);
        for b in it {
            self.push(&b);
        }
    }
}

impl FromIterator<Bond> for BondStorage {
    fn from_iter<I: IntoIterator<Item = Bond>>(iter: I) -> Self {
        let mut s = BondStorage::default();
        s.extend(iter);
        s
    }
}

//============================================================================
// Proxy
//============================================================================

/// A borrowed read-only view of one stored bond — a two-word `{storage, index}` handle.
#[derive(Debug, Clone, Copy)]
pub struct BondRef<'a> {
    st: &'a BondStorage,
    idx: usize,
}

impl BondRef<'_> {
    #[inline]
    pub fn i1(&self) -> usize {
        // SAFETY: `idx` was bounds-checked when the proxy was handed out.
        unsafe { self.st.pairs.get_unchecked(self.idx)[0] as usize }
    }

    #[inline]
    pub fn i2(&self) -> usize {
        unsafe { self.st.pairs.get_unchecked(self.idx)[1] as usize }
    }

    /// The two atom indices as a pair.
    #[inline]
    pub fn pair(&self) -> [usize; 2] {
        let p = unsafe { self.st.pairs.get_unchecked(self.idx) };
        [p[0] as usize, p[1] as usize]
    }

    /// The chemical order, or [`BondOrder::Unspecified`] when the column is absent.
    #[inline]
    pub fn order(&self) -> BondOrder {
        match self.st.orders.as_ref() {
            Some(v) => unsafe { *v.get_unchecked(self.idx) },
            None => BondOrder::Unspecified,
        }
    }

    /// Whether `idx` is one of this bond's endpoints.
    #[inline]
    pub fn contains(&self, idx: usize) -> bool {
        self.i1() == idx || self.i2() == idx
    }
}

impl From<&BondRef<'_>> for Bond {
    fn from(b: &BondRef<'_>) -> Self {
        Bond { i1: b.i1(), i2: b.i2(), order: b.order() }
    }
}

//============================================================================
// Adjacency
//============================================================================

/// One bonded neighbor: the atom at the other end, and which bond connects them.
///
/// The two fields index **different** spaces, hence the named accessors rather than a tuple.
/// They are private `u32` on purpose — public `usize` fields would double the entry to 16 bytes,
/// and there are two entries per bond.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BondNeighbor {
    atom: BondIdx,
    bond: BondIdx,
}

impl BondNeighbor {
    #[inline]
    pub fn atom(&self) -> usize {
        self.atom as usize
    }

    #[inline]
    pub fn bond(&self) -> usize {
        self.bond as usize
    }
}

/// Per-atom bonded neighbors in a compressed row layout: `offsets[i]..offsets[i+1]` slices
/// `entries` to atom `i`'s neighbors.
///
/// Self-describing — the atom count it was built for is `offsets.len() - 1`, so staleness is
/// detectable without a generation counter.
///
/// Usable **standalone**, outside any [`BondStorage`]: [`build`](Self::build) takes a bare
/// iterator of endpoint pairs, which is how `molar_ff` indexes a remapped local (`0..n`)
/// subgraph with no topology behind it.
#[derive(Debug, Clone)]
pub struct BondAdjacency {
    offsets: Vec<BondIdx>,
    entries: Vec<BondNeighbor>,
    /// Size of the bond-index space. Stored because filtering makes it unrecoverable from
    /// `entries.len()`, and because bond indices in `entries` refer to the *unfiltered* space.
    n_bonds: usize,
}

impl BondAdjacency {
    /// Index `pairs` (endpoint pairs in bond-index order) over `n_atoms` atoms.
    ///
    /// Self-bonds and out-of-range endpoints are skipped. Two passes — count degrees, then
    /// fill — hence the `Clone` bound.
    ///
    /// **Neighbor order is guaranteed**: pass 2 walks `pairs` in bond order and appends to each
    /// endpoint's run, so `neighbors(a)` comes out in ascending bond index — identical to
    /// pushing onto per-atom `Vec`s in the same loop. `molar_ff`'s GAFF port depends on this,
    /// since it indexes neighbors positionally and truncates to the first 4 or 6.
    pub fn build(n_atoms: usize, pairs: impl Iterator<Item = [usize; 2]> + Clone) -> Self {
        // Pass 1: degrees, shifted by one so the prefix sum lands directly in `offsets`.
        let mut offsets = vec![0 as BondIdx; n_atoms + 1];
        let mut n_bonds = 0usize;
        for [a, b] in pairs.clone() {
            // Counted before the filter: bond indices stay in the caller's space.
            n_bonds += 1;
            if a >= n_atoms || b >= n_atoms || a == b {
                continue;
            }
            offsets[a + 1] += 1;
            offsets[b + 1] += 1;
        }
        for i in 0..n_atoms {
            offsets[i + 1] += offsets[i];
        }

        // Pass 2: fill each atom's run in bond order.
        let total = offsets[n_atoms] as usize;
        let mut entries = vec![BondNeighbor { atom: 0, bond: 0 }; total];
        let mut fill = vec![0 as BondIdx; n_atoms];
        for (bi, [a, b]) in pairs.enumerate() {
            if a >= n_atoms || b >= n_atoms || a == b {
                continue;
            }
            let bond = to_idx(bi);
            let pa = (offsets[a] + fill[a]) as usize;
            entries[pa] = BondNeighbor { atom: to_idx(b), bond };
            fill[a] += 1;
            let pb = (offsets[b] + fill[b]) as usize;
            entries[pb] = BondNeighbor { atom: to_idx(a), bond };
            fill[b] += 1;
        }

        Self { offsets, entries, n_bonds }
    }

    /// The atom count this index was built for.
    pub fn n_atoms(&self) -> usize {
        self.offsets.len() - 1
    }

    /// The size of the bond-index space (all input pairs, including any that were filtered).
    pub fn n_bonds(&self) -> usize {
        self.n_bonds
    }

    /// Endpoints of each bond, by bond index; `[usize::MAX; 2]` for bonds that were filtered
    /// out (self-bonds, out-of-range). Recovered from `entries`, where every kept bond appears
    /// once per endpoint.
    ///
    /// Endpoint *orientation* is normalized (ascending atom index), not the caller's original
    /// `(i1, i2)`. Ring perception only needs the pair, and a reversed cycle is still a cycle.
    pub fn bond_endpoints(&self) -> Vec<[usize; 2]> {
        let mut ends = vec![[usize::MAX; 2]; self.n_bonds];
        for a in 0..self.n_atoms() {
            for nb in self.neighbors(a) {
                let e = &mut ends[nb.bond()];
                if e[0] == usize::MAX {
                    e[0] = a;
                } else {
                    e[1] = a;
                }
            }
        }
        ends
    }

    /// Atom `atom`'s bonded neighbors, in ascending bond-index order.
    ///
    /// Panics if `atom >= n_atoms()` — which is also how a stale (too-short) index surfaces.
    #[inline]
    pub fn neighbors(&self, atom: usize) -> &[BondNeighbor] {
        let s = self.offsets[atom] as usize;
        let e = self.offsets[atom + 1] as usize;
        &self.entries[s..e]
    }
}

//============================================================================
#[cfg(test)]
mod tests {
    use super::*;

    // The whole point of the u32 columns: 8 bytes per bond, 8 per neighbor entry.
    const _: () = assert!(size_of::<[BondIdx; 2]>() == 8);
    const _: () = assert!(size_of::<BondNeighbor>() == 8);

    fn storage(bonds: &[(usize, usize, BondOrder)]) -> BondStorage {
        bonds.iter().map(|&(i, j, o)| Bond::with_order(i, j, o)).collect()
    }

    fn nbr_atoms(adj: &BondAdjacency, a: usize) -> Vec<usize> {
        adj.neighbors(a).iter().map(BondNeighbor::atom).collect()
    }

    fn nbr_bonds(adj: &BondAdjacency, a: usize) -> Vec<usize> {
        adj.neighbors(a).iter().map(BondNeighbor::bond).collect()
    }

    #[test]
    fn order_column_stays_absent_without_concrete_orders() {
        let s = storage(&[(0, 1, BondOrder::Unspecified), (1, 2, BondOrder::Unspecified)]);
        assert!(!s.has_orders(), "connectivity-only source must allocate no order column");
        assert_eq!(s.get(0).unwrap().order(), BondOrder::Unspecified);
    }

    #[test]
    fn order_column_materializes_and_backfills() {
        // Third bond is the first with a concrete order; the first two backfill.
        let s = storage(&[
            (0, 1, BondOrder::Unspecified),
            (1, 2, BondOrder::Unspecified),
            (2, 3, BondOrder::Double),
        ]);
        assert!(s.has_orders());
        assert_eq!(s.get(0).unwrap().order(), BondOrder::Unspecified);
        assert_eq!(s.get(2).unwrap().order(), BondOrder::Double);
        assert_eq!(s.len(), 3);
    }

    #[test]
    fn set_order_materializes_a_missing_column() {
        let mut s = storage(&[(0, 1, BondOrder::Unspecified), (1, 2, BondOrder::Unspecified)]);
        assert!(!s.has_orders());
        s.set_order(1, BondOrder::Aromatic);
        assert!(s.has_orders());
        assert_eq!(s.get(0).unwrap().order(), BondOrder::Unspecified);
        assert_eq!(s.get(1).unwrap().order(), BondOrder::Aromatic);
    }

    #[test]
    fn proxy_reads_endpoints_and_pair() {
        let s = storage(&[(3, 7, BondOrder::Triple)]);
        let b = s.get(0).unwrap();
        assert_eq!((b.i1(), b.i2()), (3, 7));
        assert_eq!(b.pair(), [3, 7]);
        assert!(b.contains(3) && b.contains(7) && !b.contains(5));
        assert_eq!(Bond::from(&b), Bond::with_order(3, 7, BondOrder::Triple));
    }

    /// The invariant `molar_ff`'s GAFF port depends on: per-atom neighbors appear in the order
    /// their bonds were pushed, exactly as `Vec::push` into per-atom lists would give.
    #[test]
    fn neighbor_order_matches_bond_insertion_order() {
        // Atom 0 is bonded to 3, then 1, then 2 — deliberately not ascending by atom index.
        let s = storage(&[
            (0, 3, BondOrder::Single),
            (0, 1, BondOrder::Single),
            (0, 2, BondOrder::Single),
        ]);
        let adj = BondAdjacency::build(4, s.iter_pairs());
        assert_eq!(nbr_atoms(&adj, 0), vec![3, 1, 2], "neighbors follow bond order, not atom order");
        assert_eq!(nbr_bonds(&adj, 0), vec![0, 1, 2], "and carry ascending bond indices");
        // Every neighbor list is a reference back to the same bond.
        assert_eq!(nbr_atoms(&adj, 3), vec![0]);
        assert_eq!(nbr_bonds(&adj, 3), vec![0]);
    }

    #[test]
    fn build_skips_self_bonds_and_out_of_range() {
        // bond 1 is a self-bond, bond 2 points past the atom count.
        let s = storage(&[
            (0, 1, BondOrder::Single),
            (2, 2, BondOrder::Single),
            (0, 9, BondOrder::Single),
        ]);
        let adj = BondAdjacency::build(3, s.iter_pairs());
        assert_eq!(nbr_atoms(&adj, 0), vec![1], "only the valid bond survives");
        assert_eq!(nbr_atoms(&adj, 2), Vec::<usize>::new());
        // Bond indices are preserved from the original space, not renumbered by filtering.
        assert_eq!(nbr_bonds(&adj, 1), vec![0]);
    }

    #[test]
    fn adjacency_is_empty_for_no_atoms() {
        let s = storage(&[]);
        let adj = BondAdjacency::build(0, s.iter_pairs());
        assert_eq!(adj.n_atoms(), 0);
    }

    // ---- invalidation matrix (one test per row of the lifecycle table) ----

    #[test]
    fn set_order_preserves_adjacency() {
        let mut s = storage(&[(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)]);
        s.ensure_adjacency(3);
        s.set_order(0, BondOrder::Aromatic);
        assert!(
            s.get_adjacency().is_some(),
            "an order write must not cost the adjacency — this is why the columns are split"
        );
    }

    #[test]
    fn push_invalidates_adjacency() {
        let mut s = storage(&[(0, 1, BondOrder::Single)]);
        s.ensure_adjacency(3);
        s.push(&Bond::new(1, 2));
        assert!(s.get_adjacency().is_none());
    }

    #[test]
    fn remove_by_index_invalidates_adjacency() {
        let mut s = storage(&[(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)]);
        s.ensure_adjacency(3);
        s.remove_by_index(&[2], 3);
        assert!(s.get_adjacency().is_none());
    }

    #[test]
    fn invalidate_adjacency_clears_the_cache() {
        let mut s = storage(&[(0, 1, BondOrder::Single)]);
        s.ensure_adjacency(2);
        s.invalidate_adjacency();
        assert!(s.get_adjacency().is_none());
    }

    #[test]
    fn clone_drops_the_cache_but_keeps_the_columns() {
        let mut s = storage(&[(0, 1, BondOrder::Double)]);
        s.ensure_adjacency(2);
        let c = s.clone();
        assert!(c.get_adjacency().is_none(), "a derived cache is never worth copying");
        assert_eq!(c.len(), 1);
        assert_eq!(c.get(0).unwrap().order(), BondOrder::Double);
    }

    #[test]
    fn ensure_adjacency_rebuilds_when_the_atom_count_changed() {
        let mut s = storage(&[(0, 1, BondOrder::Single)]);
        assert_eq!(s.ensure_adjacency(2).n_atoms(), 2);
        assert_eq!(s.ensure_adjacency(5).n_atoms(), 5, "a different atom count is stale");
    }

    // ---- removal ----

    #[test]
    fn remove_by_index_drops_incident_bonds_and_renumbers() {
        // 5 atoms in a chain 0-1-2-3-4; remove atom 2. Bonds 1-2 and 2-3 must go, and the
        // survivors renumber onto the 4-atom space {0,1,2,3} = old {0,1,3,4}.
        let mut s = storage(&[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Double),
            (2, 3, BondOrder::Triple),
            (3, 4, BondOrder::Aromatic),
        ]);
        s.remove_by_index(&[2], 5);
        assert_eq!(s.len(), 2);
        assert_eq!(s.get(0).unwrap().pair(), [0, 1]);
        assert_eq!(s.get(0).unwrap().order(), BondOrder::Single);
        // old 3-4 becomes new 2-3, and keeps its own order (not a neighbor's).
        assert_eq!(s.get(1).unwrap().pair(), [2, 3]);
        assert_eq!(s.get(1).unwrap().order(), BondOrder::Aromatic);
    }

    #[test]
    fn remove_by_index_handles_removing_everything() {
        let mut s = storage(&[(0, 1, BondOrder::Single)]);
        s.remove_by_index(&[0, 1], 2);
        assert_eq!(s.len(), 0);
        assert!(s.is_empty());
    }
}
