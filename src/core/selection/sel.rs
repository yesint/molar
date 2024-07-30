use std::{collections::HashMap, marker::PhantomData};
use crate::prelude::*;
use anyhow::{bail, Result, anyhow};
use itertools::Itertools;

#[derive(Default)]
pub(crate) struct System {
    pub topology: Topology,
    pub state: State,
}

//---------------------------------------
// Selection
//---------------------------------------

/// Selection type that acts as a view into given set of indexes from [Topology] and [State].
/// Selections allow to query various properties of the groups of atoms and to
/// change them in different ways.
/// 
/// # Kinds of selections
/// There are three kinds of selections set by the generic marker parameter:
/// ### Serial mutable (`Sel<MutableSerial>`)
/// * May overlap
/// * Mutable
/// * Could only be accessed from the same thread where 
/// they were created (they are neither [Send] nor [Sync]). 
/// * Manipulated directly by the user.
/// ### Parallel immutable (`Sel<ImmutableParallel>`)
/// * May overlap
/// * Immutable
/// * Could be processed in parallel
/// * Not accessible directly. Manipulated by the [SourceParallel] object that created them.
/// ### Parallel mutable (`Sel<MutableParallel>`)
/// * _Can't_ overlap
/// * Mutable
/// * Could be processed in parallel.
/// * Not accessible directly. Manipulated by the [SourceParallel] object that created them.
/// 
/// Immutable selections implement the traits that allow to query properties, while 
/// mutable once also implement traits that allow to modify atoms and coordinates.
/// 
/// # Subselections
/// It is possible to select within existing selection using `select_from_*` methods. 
/// * For [ImmutableParallel] and [MutableSerial] selections subselectons have the same kind.
/// * For [MutableParallel] selections subselectons have [MutableSerial] type instead because otherwise
/// it is impossible to maintain the invariant of non-overlapping mutable access to the underlying data.
/// 
/// # Splitting selections
/// Selections could be split using the custom closire as a criterion. There are two flavours
/// of splitting functions: 
/// * `split_*` produce a number of subselections but leaves the parent selection alive.
/// The parts follow the rules of subselections.
/// * `split_into_*` consume a parent selection and produce the parts, that always have _the same_ 
/// kind as a parent selection.
pub struct Sel<K> {
    system: triomphe::Arc<System>,
    index: Vec<usize>,
    _marker: PhantomData<K>,
}

//-------------------------------------------
// Implementations of index valifdity checking
impl CheckedIndex for Sel<MutableSerial> {
    #[inline(always)]
    fn check_index(&self) -> &Vec<usize> {
        &self.index
    }
}

impl CheckedIndex for Sel<MutableParallel> {
    #[inline(always)]
    fn check_index(&self) -> &Vec<usize> {
        &self.index
    }
}

impl CheckedIndex for Sel<ImmutableParallel> {
    #[inline(always)]
    fn check_index(&self) -> &Vec<usize> {
        &self.index
    }
}

impl CheckedIndex for Sel<BuilderSerial> {
    #[inline(always)]
    fn check_index(&self) -> &Vec<usize> {
        let first = self.index[0];
        let last = self.index[&self.index.len()-1];
        let n = self.system.state.num_coords();
        if first >= n || last >= n {
            panic!(
                "Builder selection indexes [{}:{}] are out of allowed range [0:{}]",
                first,last,n
            );
        }
        &self.index
    }
}
//-------------------------------------------

impl Sel<MutableSerial> {
    // Only visible in selection module
    pub(super) fn from_parallel(sel: Sel<impl ParallelSel>) -> Self {
        Self::new(sel.system,sel.index)
    }
}

impl<K> Sel<K> 
where 
    K: SelectionKind,    
{
    // Only visible in selection module
    pub(super) fn new(system: triomphe::Arc<System>, index: Vec<usize>) -> Self {
        Self {
            system,
            index,
            _marker: PhantomData::default(),
        }
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.num_atoms()
    }

    //===================
    // Subselections
    //===================

    /// Subselection from expression
    pub fn subsel_from_expr(&self, expr: &SelectionExpr) -> Result<Sel<K::SubselKind>> {
        let index = expr.apply_subset(&self.system.topology, &self.system.state, self.index.iter().cloned())?;
        if index.len() > 0 {
            Ok(Sel::new(
                triomphe::Arc::clone(&self.system),
                index,
            ))
        } else {
            bail!("Selection is empty")
        }
    }

    /// Subselection from string
    pub fn subsel_from_str(&self, sel_str: &str) -> Result<Sel<K::SubselKind>> {
        let expr = SelectionExpr::try_from(sel_str)?;
        self.subsel_from_expr(&expr)
    }

    /// Subselection from the range of local selection indexes
    pub fn subsel_from_local_range(
        &self,
        range: std::ops::Range<usize>,
    ) -> Result<Sel<K::SubselKind>> {
        if range.end >= self.index.len() {
            bail!(
                "Invalid local sub-range: {}:{}, valid range: 0:{}",
                range.start,
                range.end,
                self.index.len() - 1
            );
        }

        // Translate range of local indexes to global indexes
        let index: Vec<usize> = self
            .index
            .iter()
            .cloned()
            .skip(range.start)
            .take(range.len())
            .collect();

        Ok(Sel::new(
            triomphe::Arc::clone(&self.system),
            index,
        ))
    }

    /// Subselection from iterator that provides local selection indexes
    pub fn subsel_from_iter(
        &self,
        iter: impl ExactSizeIterator<Item = usize>,
    ) -> Result<Sel<K::SubselKind>> {
        // Remove duplicates if any
        let index = iter
            .sorted()
            .dedup()
            .map(|i| {
                self.index.get(i).cloned().ok_or_else(|| {
                    anyhow!(
                        "Index {} is out of allowed range [0:{}]",
                        i,
                        self.index.len()
                    )
                })
            })
            .collect::<Result<Vec<usize>>>()?;
        // Now it's safe to call
        unsafe { self.subsel_from_vec_unchecked(index) }
    }

    // This method doesn't check if the vector has duplicates and thus unsafe
    pub(super) unsafe fn subsel_from_vec_unchecked<S: SelectionKind>(&self, index: Vec<usize>) -> Result<Sel<S>> {
        if index.len() > 0 {
            Ok(Sel::new(
                triomphe::Arc::clone(&self.system),
                index,
            ))
        } else {
            bail!("Selection is empty")
        }
    }

    /// Get a Particle for i-th selection index.
    /// # Safety
    /// This is an unsafe operation that doesn't check if the index is in bounds.
    pub unsafe fn nth_particle_unchecked(&self, i: usize) -> Particle<'_> {
        let ind = *self.index.get_unchecked(i);
        Particle {
            id: ind,
            atom: self.system.topology.nth_atom_unchecked(ind),
            pos: self.system.state.nth_pos_unchecked(ind),
        }
    }

    /// Get a mutable Particle for i-th selection index with
    /// # Safety
    /// This is an unsafe operation that doesn't check if the index is in bounds.
    pub unsafe fn nth_particle_unchecked_mut(&self, i: usize) -> ParticleMut<'_> {
        let ind = *self.index.get_unchecked(i);
        ParticleMut {
            id: ind,
            atom: self.system.topology.nth_atom_unchecked_mut(ind),
            pos: self.system.state.nth_pos_unchecked_mut(ind),
        }
    }

    //============================
    // Splitting
    //============================

    // Helper splitting function generic over returned selections kind
    fn split_gen<RT, F, C, Kind>(&self, func: F) -> C
    where
        RT: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(Particle) -> RT,
        C: FromIterator<Sel<Kind>> + Default,
        Kind: SelectionKind,
    {
        let mut ids = HashMap::<RT, Vec<usize>>::default();

        for p in self.iter_particle() {
            let i = p.id;
            let id = func(p);
            if let Some(el) = ids.get_mut(&id) {
                el.push(i);
            } else {
                ids.insert(id, vec![i]);
            }
        }

        C::from_iter(
            ids.into_values()
                .map(|ind| unsafe { self.subsel_from_vec_unchecked(ind).unwrap() }),
            // This should never fail because `ind` can't be empty
        )
    }

    /// Splits selection to pieces that could be disjoint
    /// according to the value of function. Parent selection is kept intact.
    /// 
    /// The number of selections correspond to the distinct values returned by `func`.
    /// Selections are stored in a container `C` and has the same kind as subselections.
    pub fn split<RT, F, C>(&self, func: F) -> C
    where
        RT: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(Particle) -> RT,
        C: FromIterator<Sel<K::SubselKind>> + Default,
    {
        self.split_gen(func)
    }

    /// Splits selection to pieces that could be disjoint
    /// according to the value of function. Parent selection is consumed.
    /// 
    /// The number of selections correspond to the distinct values returned by `func`.
    /// Selections are stored in a container `C` and has the same kind as parent selection.
    pub fn split_into<RT, F, C>(self, func: F) -> C
    where
        RT: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(Particle) -> RT,
        C: FromIterator<Sel<K>> + Default,
    {
        self.split_gen(func)
    }

    /// Helper method that splits selection into the parts with distinct resids.
    /// Parent selection is left alive.
    pub fn split_resid<C>(&self) -> C
    where
        C: FromIterator<Sel<K::SubselKind>> + Default,
    {
        self.split_gen(|p| p.atom.resid)
    }

    /// Helper method that splits selection into the parts with distinct resids.
    /// Parent selection is consumed.
    pub fn split_resid_into<C>(self) -> C
    where
        C: FromIterator<Sel<K>> + Default,
    {
        self.split_gen(|p| p.atom.resid)
    }

    /// Computes the Solvet Accessible Surface Area (SASA).
    pub fn sasa(&self) -> (f32, f32) {
        let (areas, volumes) = molar_powersasa::sasa(
            self.len(),
            0.14,
            |i| unsafe {
                let ind = *self.index.get_unchecked(i);
                self.system.state.nth_pos_unchecked_mut(ind).coords.as_mut_ptr()
            },
            |i: usize| self.nth(i).unwrap().atom.vdw(),
        );
        (areas.into_iter().sum(), volumes.into_iter().sum())
    }

    /// Returns a copy of the selection index vector.
    pub fn get_index_vec(&self) -> Vec<usize> {
        self.index.clone()
    }

    /// Saves selection to file. File type is deduced from extension.
    pub fn save(&self, fname: &str) -> Result<()> {
        let mut h = FileHandler::create(fname)?;
        h.write(self)
    }

    pub fn first_index(&self) -> usize {
        self.index[0]
    }

    pub fn last_index(&self) -> usize {
        self.index[self.index.len()-1]
    }

    /// Get a Particle for the first selection index.
    pub fn first(&self) -> Particle {
        unsafe { self.nth_particle_unchecked(0) }
    }

    /// Get a Particle for the first selection index.
    /// Index is bound-checked, an error is returned if it is out of bounds.
    pub fn nth(&self, i: usize) -> Result<Particle> {
        if i > self.len() {
            bail!("Index {} is beyond the allowed range [0:{}]", i, self.len())
        }
        Ok(unsafe { self.nth_particle_unchecked(i) })
    }

    /// Return iterator that splits selection into contigous pieces according to the value of function.
    /// Consumes a selection and returns selections of the same kind.
    /// 
    /// Whenever `func` returns a value different from the previous one, new selection is created.
    /// Selections are computed lazily when iterating.
    pub fn into_split_contig<RT, F>(self, func: F) -> IntoSelectionSplitIterator<RT, F, K>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(usize, &Atom, &Pos) -> RT,
    {
        IntoSelectionSplitIterator::new(self, func)
    }

    /// Return iterator over contigous pieces of selection with distinct contigous resids.
    /// Parent selection is consumed.
    pub fn into_split_contig_resid(
        self,
    ) -> IntoSelectionSplitIterator<i32, fn(usize, &Atom, &Pos) -> i32, K> {
        self.into_split_contig(|_, at, _| at.resid)
    }

    /// Return iterator that splits selection into contigous pieces according to the value of function.
    /// Selection is not consumed. Returned selections are of subselection kind.
    /// 
    /// Whenever `func` returns a value different from the previous one, new selection is created.
    /// Selections are computed lazily when iterating.
    pub fn split_contig<RT, F>(&self, func: F) -> SelectionSplitIterator<'_, RT, F, K>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(Particle) -> RT,
    {
        SelectionSplitIterator::new(self, func)
    }

    /// Return iterator over contigous pieces of selection with distinct contigous resids.
    /// Parent selection is left alive.
    pub fn split_contig_resid(
        &self,
    ) -> SelectionSplitIterator<'_, i32, fn(Particle) -> i32, K> {
        self.split_contig(|p| p.atom.resid)
    }

    /// Tests if two selections are from the same source
    pub fn same_source(&self, other: &Sel<K>) -> bool {
        triomphe::Arc::ptr_eq(&self.system, &other.system)
    }
}

//---------------------------------------------
/// Iterator over the [Particle]s from selection.
pub struct SelectionIterator<'a, K> {
    sel: &'a Sel<K>,
    cur: usize,
}

impl<'a, K: SelectionKind> Iterator for SelectionIterator<'a, K> {
    type Item = Particle<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.cur < self.sel.len() {
            let ret = unsafe { self.sel.nth_particle_unchecked(self.cur) };
            self.cur += 1;
            Some(ret)
        } else {
            None
        }
    }
}

impl<K: SelectionKind> ExactSizeIterator for SelectionIterator<'_, K> {
    fn len(&self) -> usize {
        self.sel.len()
    }
}

impl<K: SelectionKind> ParticleProvider for Sel<K> {
    fn iter_particle(&self) -> impl ExactSizeIterator<Item = Particle<'_>> {
        SelectionIterator { sel: self, cur: 0 }
    }
}

//---------------------------------------------
/// Mutable iterator over the [Particle]s from selection.
pub struct SelectionIteratorMut<'a, K> {
    sel: &'a Sel<K>,
    cur: usize,
}

impl<'a, K: SelectionKind> Iterator for SelectionIteratorMut<'a, K> {
    type Item = ParticleMut<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.cur < self.sel.len() {
            let ret = unsafe { self.sel.nth_particle_unchecked_mut(self.cur) };
            self.cur += 1;
            Some(ret)
        } else {
            None
        }
    }
}

impl<K: SelectionKind> ExactSizeIterator for SelectionIteratorMut<'_, K> {
    fn len(&self) -> usize {
        self.sel.len()
    }
}

impl<K: SelectionKind> ParticleMutProvider for Sel<K> {
    fn iter_particle_mut(&self) -> impl ExactSizeIterator<Item = ParticleMut<'_>> {
        SelectionIteratorMut { sel: self, cur: 0 }
    }
}



//---------------------------------------------
// Implement traits for IO

impl<K> IndexProvider for Sel<K> {
    fn iter_index(&self) -> impl Iterator<Item = usize> {
        self.index.iter().cloned()
    }
}

impl<K> TopologyProvider for Sel<K> {
    fn num_atoms(&self) -> usize {
        self.index.len()
    }
}

impl<K> StateProvider for Sel<K> {
    fn get_time(&self) -> f32 {
        self.system.state.get_time()
    }

    fn num_coords(&self) -> usize {
        self.index.len()
    }
}

//==================================================================
// Implement analysis traits

impl<K> BoxProvider for Sel<K> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.system.state.get_box()
    }
}

impl<T> MeasurePeriodic for Sel<T> {}

impl<K> PosProvider for Sel<K> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        unsafe { self.index.iter().map(|i| self.system.state.nth_pos_unchecked(*i)) }
    }
}

impl<K> MeasurePos for Sel<K> {}

impl<K> AtomsProvider for Sel<K> {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        unsafe {
            self.index
                .iter()
                .map(|i| self.system.topology.nth_atom_unchecked(*i))
        }
    }
}

impl<K> MassesProvider for Sel<K> {
    fn iter_masses(&self) -> impl ExactSizeIterator<Item = f32> {
        unsafe {
            self.index
                .iter()
                .map(|i| self.system.topology.nth_atom_unchecked(*i).mass)
        }
    }
}

impl<K> MeasureMasses for Sel<K> {}

//-------------------------------------------------------
// Mutable analysis traits only for mutable selections

impl<K: MutableSel> PosMutProvider for Sel<K> {
    fn iter_pos_mut(&self) -> impl PosMutIterator<'_> {
        unsafe {
            self.index
                .iter()
                .map(|i| self.system.state.nth_pos_unchecked_mut(*i))
        }
    }
}

impl<K: MutableSel> RandomPosMutProvider for Sel<K> {
    #[inline(always)]
    unsafe fn nth_pos_unchecked_mut(&self, i: usize) -> &mut Pos {
        let ind = *self.index.get_unchecked(i);
        self.system.state.nth_pos_unchecked_mut(ind)
    }
}

impl<K: MutableSel> ModifyPos for Sel<K> {}
impl<K: MutableSel> ModifyPeriodic for Sel<K> {}
impl<K: MutableSel> ModifyRandomAccess for Sel<K> {}

//-------------------------------------------------------
// Analysis traits for builder selections

/*
impl BoxProvider for Sel<BuilderSerial> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.state.get_box()
    }
}

impl<T> MeasurePeriodic for Sel<T> {}

*/