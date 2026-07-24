//! Struct-of-Arrays storage for atoms.
//!
//! [`AtomStorage`] replaces `Vec<Atom>` inside [`Topology`](crate::Topology): each atom
//! property is a separate column. The ten *core* properties are always-present `Vec`s of
//! length `n`; the four *optional* (force-field / chemistry) properties are
//! `Option<Vec<T>>` — `None` costs nothing, and when `Some` the column is full-length `n`
//! (mirroring the empty-`Vec`/null-pointer idiom already used for `State` velocities/forces).
//!
//! Atoms are read/written through the borrowed proxies [`AtomRef`] / [`AtomRefMut`], which
//! implement [`AtomLike`] / [`AtomLikeMut`] respectively. The owned [`Atom`] row is retained
//! as the detached construction/interchange type; [`AtomStorage::push_row`] scatters it into
//! the columns.
//!
//! ## Optional-column semantics
//! Presence is tracked **per column**, not per atom. Materializing a column (first `Some`
//! write, or [`AtomStorage::push_row`] of an atom whose optional field is `Some`) backfills
//! the other rows with the property's default. Consequently an optional getter returns
//! `None` iff the whole column is absent; once any atom carries the property, un-set atoms
//! read `Some(default)`. For every optional property the default equals its "unset" meaning
//! (empty type name, type id 0, formal charge 0, no flags), so this is information-preserving.

use crate::atom::{
    Atom, AtomFlags, AtomLike, AtomLikeMut, AtomStr, ATOM_NAME_EXPECT, ATOM_RESNAME_EXPECT,
    ATOM_TYPE_NAME_EXPECT,
};
use crate::aliases::Float;
use std::marker::PhantomData;

/// Column-wise (Struct-of-Arrays) storage for atoms. See the [module docs](self).
#[derive(Debug, Default, Clone)]
pub struct AtomStorage {
    // --- core columns: always allocated, len == n ---
    name: Vec<AtomStr>,
    resname: Vec<AtomStr>,
    resid: Vec<i32>,
    resindex: Vec<usize>,
    atomic_number: Vec<u8>,
    mass: Vec<Float>,
    charge: Vec<Float>,
    chain: Vec<char>,
    bfactor: Vec<Float>,
    occupancy: Vec<Float>,
    // --- optional columns: None => absent (zero alloc); Some => len == n ---
    type_name: Option<Vec<AtomStr>>,
    type_id: Option<Vec<u32>>,
    formal_charge: Option<Vec<i32>>,
    flags: Option<Vec<AtomFlags>>,
}

#[inline]
fn empty_atomstr() -> AtomStr {
    AtomStr::try_from_str("").unwrap()
}

/// Push into an optional column, materializing it (backfilling defaults) on the first
/// `Some`. `new_len` is the length the core columns have *after* this atom was pushed.
#[inline]
fn push_opt<T: Copy>(col: &mut Option<Vec<T>>, val: Option<T>, new_len: usize, default: T) {
    match (col.as_mut(), val) {
        (Some(v), Some(x)) => v.push(x),
        (Some(v), None) => v.push(default),
        (None, Some(x)) => {
            let mut v = vec![default; new_len - 1];
            v.push(x);
            *col = Some(v);
        }
        (None, None) => {}
    }
}

/// Retain the elements of `v` for which `keep[i]` is true. `keep.len() == v.len()`.
#[inline]
fn retain_mask<T>(v: &mut Vec<T>, keep: &[bool]) {
    let mut i = 0;
    v.retain(|_| {
        let k = keep[i];
        i += 1;
        k
    });
}

impl AtomStorage {
    /// Number of atoms (the core-column length — the single source of truth).
    pub fn len(&self) -> usize {
        self.name.len()
    }

    pub fn is_empty(&self) -> bool {
        self.name.is_empty()
    }

    /// Reserve capacity for `n` more atoms in every present column.
    pub fn reserve(&mut self, n: usize) {
        self.name.reserve(n);
        self.resname.reserve(n);
        self.resid.reserve(n);
        self.resindex.reserve(n);
        self.atomic_number.reserve(n);
        self.mass.reserve(n);
        self.charge.reserve(n);
        self.chain.reserve(n);
        self.bfactor.reserve(n);
        self.occupancy.reserve(n);
        if let Some(c) = self.type_name.as_mut() {
            c.reserve(n);
        }
        if let Some(c) = self.type_id.as_mut() {
            c.reserve(n);
        }
        if let Some(c) = self.formal_charge.as_mut() {
            c.reserve(n);
        }
        if let Some(c) = self.flags.as_mut() {
            c.reserve(n);
        }
    }

    /// Append an atom, scattering its fields into the columns. Optional columns materialize
    /// on the first atom that carries the corresponding property (see [module docs](self)).
    pub fn push_row(&mut self, a: &Atom) {
        self.name.push(a.name);
        self.resname.push(a.resname);
        self.resid.push(a.resid);
        self.resindex.push(a.resindex);
        self.atomic_number.push(a.atomic_number);
        self.mass.push(a.mass);
        self.charge.push(a.charge);
        self.chain.push(a.chain);
        self.bfactor.push(a.bfactor);
        self.occupancy.push(a.occupancy);
        let n = self.name.len();
        push_opt(&mut self.type_name, a.type_name, n, empty_atomstr());
        push_opt(&mut self.type_id, a.type_id, n, 0);
        push_opt(&mut self.formal_charge, a.formal_charge, n, 0);
        push_opt(&mut self.flags, a.flags, n, AtomFlags::default());
        debug_assert!(self.invariant_holds());
    }

    /// Overwrite the atom at `i` from an owned [`Atom`]. Optional columns materialize if the
    /// incoming atom carries a property the column doesn't yet hold; an incoming `None` on a
    /// present column writes the property's default.
    pub fn set_row(&mut self, i: usize, a: &Atom) {
        let n = self.len();
        self.name[i] = a.name;
        self.resname[i] = a.resname;
        self.resid[i] = a.resid;
        self.resindex[i] = a.resindex;
        self.atomic_number[i] = a.atomic_number;
        self.mass[i] = a.mass;
        self.charge[i] = a.charge;
        self.chain[i] = a.chain;
        self.bfactor[i] = a.bfactor;
        self.occupancy[i] = a.occupancy;
        set_opt(&mut self.type_name, i, a.type_name, n, empty_atomstr());
        set_opt(&mut self.type_id, i, a.type_id, n, 0);
        set_opt(&mut self.formal_charge, i, a.formal_charge, n, 0);
        set_opt(&mut self.flags, i, a.flags, n, AtomFlags::default());
        debug_assert!(self.invariant_holds());
    }

    /// Remove the atoms at the given (sorted, unique) global indices from every column.
    pub fn retain_by_index(&mut self, removed: &[usize]) {
        let n = self.len();
        let mut keep = vec![true; n];
        for &i in removed {
            if i < n {
                keep[i] = false;
            }
        }
        retain_mask(&mut self.name, &keep);
        retain_mask(&mut self.resname, &keep);
        retain_mask(&mut self.resid, &keep);
        retain_mask(&mut self.resindex, &keep);
        retain_mask(&mut self.atomic_number, &keep);
        retain_mask(&mut self.mass, &keep);
        retain_mask(&mut self.charge, &keep);
        retain_mask(&mut self.chain, &keep);
        retain_mask(&mut self.bfactor, &keep);
        retain_mask(&mut self.occupancy, &keep);
        if let Some(c) = self.type_name.as_mut() {
            retain_mask(c, &keep);
        }
        if let Some(c) = self.type_id.as_mut() {
            retain_mask(c, &keep);
        }
        if let Some(c) = self.formal_charge.as_mut() {
            retain_mask(c, &keep);
        }
        if let Some(c) = self.flags.as_mut() {
            retain_mask(c, &keep);
        }
        debug_assert!(self.invariant_holds());
    }

    /// A read-only proxy for atom `i`.
    ///
    /// # Safety
    /// `i` must be `< len()`: the returned proxy indexes the columns unchecked, so *using* it
    /// with an out-of-bounds `i` is undefined behavior.
    #[inline]
    pub unsafe fn get_unchecked(&self, i: usize) -> AtomRef<'_> {
        AtomRef { st: self, idx: i }
    }

    /// A read-only proxy for atom `i`, or `None` if out of bounds.
    #[inline]
    pub fn get(&self, i: usize) -> Option<AtomRef<'_>> {
        (i < self.len()).then_some(AtomRef { st: self, idx: i })
    }

    /// A mutable proxy for atom `i`.
    ///
    /// # Safety
    /// `i` must be `< len()` (see [`get_unchecked`](Self::get_unchecked)).
    #[inline]
    pub unsafe fn get_mut_unchecked(&mut self, i: usize) -> AtomRefMut<'_> {
        AtomRefMut {
            st: self as *mut AtomStorage,
            idx: i,
            _pd: PhantomData,
        }
    }

    /// A mutable proxy for atom `i`, or `None` if out of bounds.
    #[inline]
    pub fn get_mut(&mut self, i: usize) -> Option<AtomRefMut<'_>> {
        if i < self.len() {
            Some(unsafe { self.get_mut_unchecked(i) })
        } else {
            None
        }
    }

    /// Reconstruct an owned [`Atom`] from row `i`. Panics if `i` is out of bounds.
    pub fn to_atom(&self, i: usize) -> Atom {
        Atom::from(&self.get(i).expect("atom index out of bounds"))
    }

    /// Iterate read-only proxies over all atoms.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = AtomRef<'_>> {
        // SAFETY: `i` ranges over `0..len()`.
        (0..self.len()).map(move |i| unsafe { self.get_unchecked(i) })
    }

    // --- bulk column slices (read-only) ---
    // Direct access to a whole always-present column. Used by the selection hot path to scan
    // one contiguous column instead of materializing a proxy per atom (a cache-locality win),
    // and generally useful for vectorized read access.
    pub fn names(&self) -> &[AtomStr] {
        &self.name
    }
    pub fn resnames(&self) -> &[AtomStr] {
        &self.resname
    }
    pub fn resids(&self) -> &[i32] {
        &self.resid
    }
    pub fn resindices(&self) -> &[usize] {
        &self.resindex
    }
    pub fn atomic_numbers(&self) -> &[u8] {
        &self.atomic_number
    }
    pub fn masses(&self) -> &[Float] {
        &self.mass
    }
    pub fn charges(&self) -> &[Float] {
        &self.charge
    }
    pub fn chains(&self) -> &[char] {
        &self.chain
    }
    pub fn bfactors(&self) -> &[Float] {
        &self.bfactor
    }
    pub fn occupancies(&self) -> &[Float] {
        &self.occupancy
    }

    // --- optional-column materializers (no-op if already present) ---
    fn ensure_type_name(&mut self) -> &mut Vec<AtomStr> {
        let n = self.name.len();
        self.type_name.get_or_insert_with(|| vec![empty_atomstr(); n])
    }
    fn ensure_type_id(&mut self) -> &mut Vec<u32> {
        let n = self.name.len();
        self.type_id.get_or_insert_with(|| vec![0; n])
    }
    fn ensure_formal_charge(&mut self) -> &mut Vec<i32> {
        let n = self.name.len();
        self.formal_charge.get_or_insert_with(|| vec![0; n])
    }
    fn ensure_flags(&mut self) -> &mut Vec<AtomFlags> {
        let n = self.name.len();
        self.flags.get_or_insert_with(|| vec![AtomFlags::default(); n])
    }

    /// Debug check: every present optional column matches the core length.
    fn invariant_holds(&self) -> bool {
        let n = self.len();
        self.type_name.as_ref().is_none_or(|c| c.len() == n)
            && self.type_id.as_ref().is_none_or(|c| c.len() == n)
            && self.formal_charge.as_ref().is_none_or(|c| c.len() == n)
            && self.flags.as_ref().is_none_or(|c| c.len() == n)
    }
}

/// Write into an optional column at index `i`, materializing on a `Some` write.
#[inline]
fn set_opt<T: Copy>(col: &mut Option<Vec<T>>, i: usize, val: Option<T>, n: usize, default: T) {
    match (col.as_mut(), val) {
        (Some(v), Some(x)) => v[i] = x,
        (Some(v), None) => v[i] = default,
        (None, Some(x)) => {
            let mut v = vec![default; n];
            v[i] = x;
            *col = Some(v);
        }
        (None, None) => {}
    }
}

impl Extend<Atom> for AtomStorage {
    fn extend<I: IntoIterator<Item = Atom>>(&mut self, iter: I) {
        let it = iter.into_iter();
        self.reserve(it.size_hint().0);
        for a in it {
            self.push_row(&a);
        }
    }
}

impl<'a> Extend<&'a Atom> for AtomStorage {
    fn extend<I: IntoIterator<Item = &'a Atom>>(&mut self, iter: I) {
        for a in iter {
            self.push_row(a);
        }
    }
}

impl FromIterator<Atom> for AtomStorage {
    fn from_iter<I: IntoIterator<Item = Atom>>(iter: I) -> Self {
        let mut s = AtomStorage::default();
        s.extend(iter);
        s
    }
}

//============================================================================
// Proxies
//============================================================================

/// Read-only view of a single atom stored column-wise. A two-word `{&AtomStorage, index}`
/// handle; implements [`AtomLike`].
#[derive(Debug, Clone, Copy)]
pub struct AtomRef<'a> {
    st: &'a AtomStorage,
    idx: usize,
}

impl<'a> AtomRef<'a> {
    /// Atom name, borrowed from the backing storage (lifetime `'a`, not `&self`).
    #[inline]
    pub fn name(&self) -> &'a str {
        let st: &'a AtomStorage = self.st;
        unsafe { st.name.get_unchecked(self.idx).as_str() }
    }
    /// Residue name, borrowed from the backing storage (lifetime `'a`).
    #[inline]
    pub fn resname(&self) -> &'a str {
        let st: &'a AtomStorage = self.st;
        unsafe { st.resname.get_unchecked(self.idx).as_str() }
    }

    /// The backing storage (used by length-1 providers such as `Particle`).
    #[inline]
    pub(crate) fn storage(&self) -> &'a AtomStorage {
        self.st
    }
}

impl AtomLike for AtomRef<'_> {
    fn get_name(&self) -> &str {
        unsafe { self.st.name.get_unchecked(self.idx).as_str() }
    }
    fn get_resname(&self) -> &str {
        unsafe { self.st.resname.get_unchecked(self.idx).as_str() }
    }
    fn get_resid(&self) -> isize {
        unsafe { *self.st.resid.get_unchecked(self.idx) as isize }
    }
    fn get_resindex(&self) -> usize {
        unsafe { *self.st.resindex.get_unchecked(self.idx) }
    }
    fn get_atomic_number(&self) -> u8 {
        unsafe { *self.st.atomic_number.get_unchecked(self.idx) }
    }
    fn get_mass(&self) -> Float {
        unsafe { *self.st.mass.get_unchecked(self.idx) }
    }
    fn get_charge(&self) -> Float {
        unsafe { *self.st.charge.get_unchecked(self.idx) }
    }
    fn get_chain(&self) -> char {
        unsafe { *self.st.chain.get_unchecked(self.idx) }
    }
    fn get_bfactor(&self) -> Float {
        unsafe { *self.st.bfactor.get_unchecked(self.idx) }
    }
    fn get_occupancy(&self) -> Float {
        unsafe { *self.st.occupancy.get_unchecked(self.idx) }
    }
    fn get_type_name(&self) -> Option<&str> {
        self.st
            .type_name
            .as_ref()
            .map(|c| unsafe { c.get_unchecked(self.idx).as_str() })
    }
    fn get_type_id(&self) -> Option<u32> {
        self.st.type_id.as_ref().map(|c| unsafe { *c.get_unchecked(self.idx) })
    }
    fn get_formal_charge(&self) -> Option<i32> {
        self.st
            .formal_charge
            .as_ref()
            .map(|c| unsafe { *c.get_unchecked(self.idx) })
    }
    fn get_flags(&self) -> Option<AtomFlags> {
        self.st.flags.as_ref().map(|c| unsafe { *c.get_unchecked(self.idx) })
    }
}

/// Mutable view of a single atom stored column-wise. A two-word `{*mut AtomStorage, index}`
/// handle; implements [`AtomLike`] + [`AtomLikeMut`].
///
/// # Safety / parallel use
/// Setters for the *core* columns read the column's `Vec` header through a **shared** borrow
/// (`&AtomStorage`) to obtain the heap-buffer base pointer, then write a single element through
/// a raw pointer — they never form a `&mut [T]`/`&mut Vec` over the whole column. So across a
/// parallel split the `AtomStorage` struct/headers are only ever *read* while disjoint heap
/// elements are written, which is race-free (verified with `cargo miri` under Tree Borrows —
/// see `par_atom_column_write_scoped`). Setters for the *optional* columns may materialize the
/// column (allocate + backfill) and therefore take `&mut AtomStorage`; they must only be used
/// **serially** — materialize optional columns before entering any parallel region.
#[derive(Debug)]
pub struct AtomRefMut<'a> {
    st: *mut AtomStorage,
    idx: usize,
    _pd: PhantomData<&'a mut AtomStorage>,
}

impl<'a> From<AtomRefMut<'a>> for AtomRef<'a> {
    fn from(m: AtomRefMut<'a>) -> Self {
        AtomRef { st: unsafe { &*m.st }, idx: m.idx }
    }
}

impl<'a> AtomRefMut<'a> {
    /// Construct a mutable proxy from a raw storage pointer and global index.
    ///
    /// # Safety
    /// `st` must be valid for `'a`, `idx` must be in bounds, and no other live proxy or
    /// reference may alias the same element mutably (disjoint indices across a parallel
    /// split satisfy this).
    #[inline]
    pub(crate) unsafe fn from_raw(st: *mut AtomStorage, idx: usize) -> Self {
        AtomRefMut { st, idx, _pd: PhantomData }
    }

    #[inline]
    fn st(&self) -> &AtomStorage {
        unsafe { &*self.st }
    }

    /// Raw mutable pointer to the start of a *core* column's heap buffer, obtained by reading
    /// the `Vec` header through a **shared** borrow (never a `&mut [T]`/`&mut Vec` over the
    /// whole column). Writing a disjoint element through the returned pointer is race-free
    /// across threads: the parallel region only writes heap buffers, never the `AtomStorage`
    /// struct/headers (which are only read). See the type-level docs.
    ///
    /// # Safety
    /// `self.idx` must be `< len()`; the caller must not let two proxies write the same element.
    #[inline]
    unsafe fn core_col_mut<T>(&self, proj: impl Fn(&AtomStorage) -> &Vec<T>) -> *mut T {
        proj(self.st()).as_ptr() as *mut T
    }
}

// SAFETY: the raw storage pointer is only dereferenced for element `idx`. Distinct proxies
// over disjoint indices (as produced by a parallel split) never form overlapping references,
// matching the existing `SelParMut` discipline. Optional-column setters that materialize a
// column take `&mut AtomStorage` and must only be used serially (see SoA plan Stage 4).
unsafe impl Send for AtomRefMut<'_> {}
unsafe impl Sync for AtomRefMut<'_> {}

impl AtomLike for AtomRefMut<'_> {
    fn get_name(&self) -> &str {
        unsafe { self.st().name.get_unchecked(self.idx).as_str() }
    }
    fn get_resname(&self) -> &str {
        unsafe { self.st().resname.get_unchecked(self.idx).as_str() }
    }
    fn get_resid(&self) -> isize {
        unsafe { *self.st().resid.get_unchecked(self.idx) as isize }
    }
    fn get_resindex(&self) -> usize {
        unsafe { *self.st().resindex.get_unchecked(self.idx) }
    }
    fn get_atomic_number(&self) -> u8 {
        unsafe { *self.st().atomic_number.get_unchecked(self.idx) }
    }
    fn get_mass(&self) -> Float {
        unsafe { *self.st().mass.get_unchecked(self.idx) }
    }
    fn get_charge(&self) -> Float {
        unsafe { *self.st().charge.get_unchecked(self.idx) }
    }
    fn get_chain(&self) -> char {
        unsafe { *self.st().chain.get_unchecked(self.idx) }
    }
    fn get_bfactor(&self) -> Float {
        unsafe { *self.st().bfactor.get_unchecked(self.idx) }
    }
    fn get_occupancy(&self) -> Float {
        unsafe { *self.st().occupancy.get_unchecked(self.idx) }
    }
    fn get_type_name(&self) -> Option<&str> {
        self.st()
            .type_name
            .as_ref()
            .map(|c| unsafe { c.get_unchecked(self.idx).as_str() })
    }
    fn get_type_id(&self) -> Option<u32> {
        self.st().type_id.as_ref().map(|c| unsafe { *c.get_unchecked(self.idx) })
    }
    fn get_formal_charge(&self) -> Option<i32> {
        self.st()
            .formal_charge
            .as_ref()
            .map(|c| unsafe { *c.get_unchecked(self.idx) })
    }
    fn get_flags(&self) -> Option<AtomFlags> {
        self.st().flags.as_ref().map(|c| unsafe { *c.get_unchecked(self.idx) })
    }
}

impl AtomLikeMut for AtomRefMut<'_> {
    fn set_name(&mut self, name: &str) {
        let s = AtomStr::try_from_str(name).expect(ATOM_NAME_EXPECT);
        unsafe { *self.core_col_mut(|st| &st.name).add(self.idx) = s }
    }
    fn set_resname(&mut self, resname: &str) {
        let s = AtomStr::try_from_str(resname).expect(ATOM_RESNAME_EXPECT);
        unsafe { *self.core_col_mut(|st| &st.resname).add(self.idx) = s }
    }
    fn set_resid(&mut self, resid: isize) {
        unsafe { *self.core_col_mut(|st| &st.resid).add(self.idx) = resid as i32 }
    }
    fn set_resindex(&mut self, resindex: usize) {
        unsafe { *self.core_col_mut(|st| &st.resindex).add(self.idx) = resindex }
    }
    fn set_atomic_number(&mut self, atomic_number: u8) {
        unsafe { *self.core_col_mut(|st| &st.atomic_number).add(self.idx) = atomic_number }
    }
    fn set_mass(&mut self, mass: Float) {
        unsafe { *self.core_col_mut(|st| &st.mass).add(self.idx) = mass }
    }
    fn set_charge(&mut self, charge: Float) {
        unsafe { *self.core_col_mut(|st| &st.charge).add(self.idx) = charge }
    }
    fn set_chain(&mut self, chain: char) {
        unsafe { *self.core_col_mut(|st| &st.chain).add(self.idx) = chain }
    }
    fn set_bfactor(&mut self, bfactor: Float) {
        unsafe { *self.core_col_mut(|st| &st.bfactor).add(self.idx) = bfactor }
    }
    fn set_occupancy(&mut self, occupancy: Float) {
        unsafe { *self.core_col_mut(|st| &st.occupancy).add(self.idx) = occupancy }
    }
    fn set_type_name(&mut self, type_name: &str) {
        let s = AtomStr::try_from_str(type_name).expect(ATOM_TYPE_NAME_EXPECT);
        let idx = self.idx;
        unsafe {
            let col = (&mut *self.st).ensure_type_name();
            *col.get_unchecked_mut(idx) = s;
        }
    }
    fn set_type_id(&mut self, type_id: u32) {
        let idx = self.idx;
        unsafe {
            let col = (&mut *self.st).ensure_type_id();
            *col.get_unchecked_mut(idx) = type_id;
        }
    }
    fn set_formal_charge(&mut self, formal_charge: i32) {
        let idx = self.idx;
        unsafe {
            let col = (&mut *self.st).ensure_formal_charge();
            *col.get_unchecked_mut(idx) = formal_charge;
        }
    }
    fn set_flags(&mut self, flags: AtomFlags) {
        let idx = self.idx;
        unsafe {
            let col = (&mut *self.st).ensure_flags();
            *col.get_unchecked_mut(idx) = flags;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atom(name: &str, resid: i32) -> Atom {
        Atom::new().with_name(name).with_resname("MOL").with_resid(resid)
    }

    #[test]
    fn push_and_read_core() {
        let mut s = AtomStorage::default();
        s.push_row(&atom("CA", 1).with_mass(12.0).with_chain('A'));
        s.push_row(&atom("CB", 2).with_mass(14.0));
        assert_eq!(s.len(), 2);
        assert_eq!(s.get(0).unwrap().get_name(), "CA");
        assert_eq!(s.get(0).unwrap().get_mass(), 12.0);
        assert_eq!(s.get(0).unwrap().get_chain(), 'A');
        assert_eq!(s.get(1).unwrap().get_resid(), 2);
        assert!(s.get(2).is_none());
    }

    #[test]
    fn optional_column_absent_until_set() {
        let mut s = AtomStorage::default();
        s.push_row(&atom("C", 1));
        s.push_row(&atom("N", 1));
        // No atom carried a type id → column absent → getter is None.
        assert_eq!(s.get(0).unwrap().get_type_id(), None);
        assert!(s.type_id.is_none());
    }

    #[test]
    fn optional_materializes_on_first_some_and_backfills() {
        let mut s = AtomStorage::default();
        s.push_row(&atom("C", 1)); // no type id
        s.push_row(&atom("N", 1).with_type_id(7)); // triggers materialization
        s.push_row(&atom("O", 1)); // pushed after column exists → default
        // Column now present and full-length; row 0 backfilled to default 0.
        assert_eq!(s.get(0).unwrap().get_type_id(), Some(0));
        assert_eq!(s.get(1).unwrap().get_type_id(), Some(7));
        assert_eq!(s.get(2).unwrap().get_type_id(), Some(0));
        assert_eq!(s.type_id.as_ref().unwrap().len(), 3);
    }

    #[test]
    fn mut_proxy_sets_core_and_optional() {
        let mut s = AtomStorage::default();
        s.push_row(&atom("C", 1));
        s.push_row(&atom("N", 2));
        {
            let mut m = s.get_mut(1).unwrap();
            m.set_mass(14.5);
            m.set_type_name("n3"); // materializes type_name
            m.set_in_ring(true); // materializes flags via get_flags/set_flags default
        }
        assert_eq!(s.get(1).unwrap().get_mass(), 14.5);
        assert_eq!(s.get(1).unwrap().get_type_name(), Some("n3"));
        assert!(s.get(1).unwrap().is_in_ring());
        // Backfilled sibling: column present, default values.
        assert_eq!(s.get(0).unwrap().get_type_name(), Some(""));
        assert!(!s.get(0).unwrap().is_in_ring());
    }

    #[test]
    fn retain_by_index_keeps_columns_in_sync() {
        let mut s = AtomStorage::default();
        for i in 0..5 {
            s.push_row(&atom("C", i).with_type_id(i as u32 + 10));
        }
        s.retain_by_index(&[1, 3]); // remove atoms 1 and 3
        assert_eq!(s.len(), 3);
        let ids: Vec<u32> = (0..3).map(|i| s.get(i).unwrap().get_type_id().unwrap()).collect();
        assert_eq!(ids, vec![10, 12, 14]);
        let resids: Vec<isize> = (0..3).map(|i| s.get(i).unwrap().get_resid()).collect();
        assert_eq!(resids, vec![0, 2, 4]);
    }

    #[test]
    fn roundtrip_to_atom() {
        let mut s = AtomStorage::default();
        s.push_row(&atom("CA", 7).with_type_name("ca").with_formal_charge(-1));
        let a = s.to_atom(0);
        assert_eq!(a.get_name(), "CA");
        assert_eq!(a.get_resid(), 7);
        assert_eq!(a.get_type_name(), Some("ca"));
        assert_eq!(a.get_formal_charge(), Some(-1));
    }
}
