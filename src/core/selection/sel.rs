use crate::prelude::*;

use rayon::iter::IntoParallelIterator;
use sorted_vec::SortedSet;
use std::collections::HashMap;

use super::utils::*;

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

pub struct Sel<K: SelectionKind> {
    topology: Holder<Topology, K>,
    state: Holder<State, K>,
    index_storage: SortedSet<usize>,
}

//-------------------------------------------
// Functions shared by all selection kinds
//-------------------------------------------

impl<K: SelectionKind> Sel<K> {
    #[inline(always)]
    fn check_index(&self) -> Result<(), SelectionError> {
        K::check_index(&self.index_storage, &self.topology, &self.state)
    }

    pub(crate) fn from_holders_and_index(
        top_holder: Holder<Topology, K>,
        st_holder: Holder<State, K>,
        index: SortedSet<usize>,
    ) -> Result<Self, SelectionError> {
        let sel = Sel {
            topology: top_holder,
            state: st_holder,
            index_storage: index,
        };
        sel.check_index()?;
        Ok(sel)
    }

    // If selection is MutableParallel it will clear used indexes
    // when dropped. This invalidates used indexes in certan situations
    // like splitting and combining selections because indexes are already
    // used by the split fragments or combined selection respectively.
    // To avoid this we can clear index manually so that nothing is cleared
    // upon dropping selection.
    pub(crate) unsafe fn clear_index_before_drop(&mut self) {
        self.index_storage.clear();
    }

    #[inline(always)]
    pub(crate) fn index(&self) -> &SortedSet<usize> {
        self.check_index().unwrap();
        &self.index_storage
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        // We need to check index here manually
        self.check_index().unwrap();
        self.num_atoms()
    }

    // This method doesn't check if the vector is sorted and has duplicates and thus unsafe
    pub(super) unsafe fn subsel_from_sorted_vec_unchecked<KO: SelectionKind>(
        &self,
        index: Vec<usize>,
    ) -> Result<Sel<KO>, SelectionError> {
        if index.len() > 0 {
            Ok(Sel {
                topology: self.topology.clone_into(),
                state: self.state.clone_into(),
                index_storage: SortedSet::from_sorted(index),
            })
        } else {
            Err(SelectionError::FromVec {
                first: index[0],
                last: index[index.len() - 1],
                size: index.len(),
                source: SelectionIndexError::IndexEmpty,
            })
        }
    }

    pub(super) unsafe fn subsel_from_unsorted_vec_unchecked<KO: SelectionKind>(
        &self,
        index: Vec<usize>,
    ) -> Result<Sel<KO>, SelectionError> {
        if index.len() > 0 {
            Ok(Sel {
                topology: self.topology.clone_into(),
                state: self.state.clone_into(),
                index_storage: SortedSet::from_unsorted(index),
            })
        } else {
            Err(SelectionError::FromVec {
                first: index[0],
                last: index[index.len() - 1],
                size: index.len(),
                source: SelectionIndexError::IndexEmpty,
            })
        }
    }

    /// Get a Particle for i-th selection index.
    /// # Safety
    /// This is an unsafe operation that doesn't check if the index is in bounds.
    pub(super) unsafe fn nth_particle_unchecked(&self, i: usize) -> Particle<'_> {
        let ind = *self.index_storage.get_unchecked(i);
        Particle {
            id: ind,
            atom: self.topology.nth_atom_unchecked(ind),
            pos: self.state.nth_pos_unchecked(ind),
        }
    }

    /// Get a mutable Particle for i-th selection index with
    /// # Safety
    /// This is an unsafe operation that doesn't check if the index is in bounds.
    pub(super) unsafe fn nth_particle_unchecked_mut(&self, i: usize) -> ParticleMut<'_> {
        let ind = *self.index_storage.get_unchecked(i);
        ParticleMut {
            id: ind,
            atom: self.topology.nth_atom_unchecked_mut(ind),
            pos: self.state.nth_pos_unchecked_mut(ind),
        }
    }

    pub fn nth_pos(&self, i: usize) -> Option<&Pos> {
        let ind = *self.index().get(i)?;
        Some(unsafe { self.state.nth_pos_unchecked(ind) })
    }

    pub fn nth_pos_mut(&self, i: usize) -> Option<&mut Pos> {
        let ind = *self.index().get(i)?;
        Some(unsafe { self.state.nth_pos_unchecked_mut(ind) })
    }

    //============================
    // Splitting
    //============================

    // Helper splitting function doing actual work
    fn split_gen<RT, F, C, KO>(&self, func: F) -> C
    where
        RT: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(Particle) -> Option<RT>,
        C: FromIterator<Sel<KO>> + Default,
        KO: SelectionKind,
    {
        let mut ids = HashMap::<RT, Vec<usize>>::default();

        for p in self.iter_particle() {
            let i = p.id;
            if let Some(id) = func(p) {
                if let Some(el) = ids.get_mut(&id) {
                    el.push(i);
                } else {
                    ids.insert(id, vec![i]);
                }
            }
        }

        C::from_iter(
            ids.into_values()
                .map(|ind| unsafe { self.subsel_from_unsorted_vec_unchecked(ind).unwrap() }),
            // This should never fail because `ind` can't be empty and can't contain duplicates
        )
    }

    /// Splits selection to pieces that could be disjoint
    /// according to the value of function. Parent selection is consumed.    
    ///
    /// The number of selections correspond to the distinct values returned by `func`.
    /// Selections are stored in a container `C` and has the same kind as parent selection.
    pub fn into_split_disjoint<RT, F, C>(self, func: F) -> C
    where
        RT: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(Particle) -> Option<RT>,
        C: FromIterator<Sel<K>> + Default,
    {
        self.split_gen(func)
    }

    /// Helper method that splits selection into the parts with distinct resids.
    /// Parent selection is consumed.
    pub fn into_split_resid<C>(self) -> C
    where
        C: FromIterator<Sel<K>> + Default,
    {
        self.split_gen(|p| Some(p.atom.resid))
    }

    // pub fn into_fragments(self, func: impl Fn(Self)->Vec<Sel<K>>)->Vec<Sel<K>> {
    //     // Transform self into Serial
    //     func(self)
    // }

    /// Computes the Solvet Accessible Surface Area (SASA).
    pub fn sasa(&self) -> (f32, f32) {
        let (areas, volumes) = molar_powersasa::sasa(
            self.len(),
            0.14,
            |i| unsafe {
                let ind = *self.index().get_unchecked(i);
                self.state.nth_pos_unchecked_mut(ind).coords.as_mut_ptr()
            },
            |i: usize| self.nth_particle(i).unwrap().atom.vdw(),
        );
        (areas.into_iter().sum(), volumes.into_iter().sum())
    }

    /// Returns a copy of the selection index vector.
    pub fn get_index_vec(&self) -> SortedSet<usize> {
        self.index().clone()
    }

    pub fn first_index(&self) -> usize {
        self.index()[0]
    }

    pub fn last_index(&self) -> usize {
        self.index()[self.index().len() - 1]
    }

    pub fn nth_index(&self, i: usize) -> Option<usize> {
        self.index().get(i).cloned()
    }

    /// Get a Particle for the first selection index.
    pub fn first_particle(&self) -> Particle {
        unsafe { self.nth_particle_unchecked(0) }
    }

    /// Get a Particle for the last selection index.
    pub fn last_particle(&self) -> Particle {
        unsafe { self.nth_particle_unchecked(self.index().len() - 1) }
    }

    /// Get a Particle for the n-th selection index.
    /// Index is bound-checked, an error is returned if it is out of bounds.
    pub fn nth_particle(&self, i: usize) -> Option<Particle> {
        if i >= self.len() {
            None
        } else {
            Some(unsafe { self.nth_particle_unchecked(i) })
        }
    }

    pub fn nth_particle_mut(&self, i: usize) -> Option<ParticleMut> {
        if i >= self.len() {
            None
        } else {
            Some(unsafe { self.nth_particle_unchecked_mut(i) })
        }
    }

    /// Return iterator that splits selection into contigous pieces according to the value of function.
    /// Consumes a selection and returns selections of the same kind.
    ///
    /// Whenever `func` returns a value different from the previous one, new selection is created.
    /// Selections are computed lazily when iterating.
    pub fn into_iter_contig_fragments<RT, F>(self, func: F) -> IntoFragmentsIterator<RT, F, K>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(Particle) -> Option<RT>,
    {
        IntoFragmentsIterator::new(self, func)
    }

    /// Return iterator over contigous pieces of selection with distinct contigous resids.
    /// Parent selection is consumed.
    pub fn into_iter_contig_fragments_resindex(
        self,
    ) -> IntoFragmentsIterator<usize, impl Fn(Particle) -> Option<usize>, K> {
        self.into_iter_contig_fragments(|p| Some(p.atom.resindex))
    }

    pub fn into_iter_contig_fragments_chain(
        self,
    ) -> IntoFragmentsIterator<char, impl Fn(Particle) -> Option<char>, K> {
        self.into_iter_contig_fragments(|p| Some(p.atom.chain))
    }
}

impl<K: SelectionKind> Drop for Sel<K> {
    fn drop(&mut self) {
        unsafe {
            K::remove_used(&self.index_storage, &self.topology.get_used());
            K::remove_used(&self.index_storage, &self.state.get_used());
        }
    }
}

//-------------------------------------------
// For all selections except MutableParallel
//-------------------------------------------

impl<K> Sel<K>
where
    K: SelectionKind<UsedIndexesType = ()>,
{
    /// Tests if two selections are from the same source
    pub fn is_same_source(&self, other: &Sel<K>) -> bool {
        self.topology.same_data(&other.topology) && self.state.same_data(&other.state)
    }

    //===================
    // Subselections
    //===================

    fn subselect_internal(&self, index: SortedSet<usize>) -> Result<Self, SelectionError> {
        Self::from_holders_and_index(
            self.topology.clone_with_index(&index)?,
            self.state.clone_with_index(&index)?,
            index,
        )
    }

    /// Subselection from expression
    pub fn subsel_expr(&self, expr: &mut SelectionExpr) -> Result<Sel<K>, SelectionError> {
        let index = index_from_expr_sub(expr, &self.topology, &self.state, &self.index())?;
        Self::from_holders_and_index(
            self.topology.clone_with_index(&index)?,
            self.state.clone_with_index(&index)?,
            index,
        )
    }

    /// Subselection from string
    pub fn subsel_str(&self, sel_str: impl AsRef<str>) -> Result<Sel<K>, SelectionError> {
        let mut expr = SelectionExpr::try_from(sel_str.as_ref())?;
        self.subsel_expr(&mut expr)
    }

    /// Subselection from the range of local selection indexes
    pub fn subsel_local_range(
        &self,
        range: std::ops::Range<usize>,
    ) -> Result<Sel<K>, SelectionError> {
        if range.end >= self.index().len() {
            return Err(SelectionError::LocalRange(
                range.start,
                range.end,
                self.index().len() - 1,
            ));
        }

        // Translate range of local indexes to global indexes
        let index: Vec<usize> = self
            .index()
            .iter()
            .cloned()
            .skip(range.start)
            .take(range.len())
            .collect();

        self.subselect_internal(unsafe { SortedSet::from_sorted(index) })
    }

    /// Subselection from iterator that provides local selection indexes
    pub fn subsel_iter(
        &self,
        iter: impl ExactSizeIterator<Item = usize>,
    ) -> Result<Sel<K>, SelectionError> {
        let ind = self.index();
        let global_vec: Vec<_> = iter.map(|i| ind[i]).collect();
        self.subselect_internal(index_from_vec(global_vec, self.len())?)
    }

    //==============================================
    // Splitting that keeps parent selection alive
    //==============================================

    /// Splits selection to pieces that could be disjoint
    /// according to the value of function. Parent selection is kept intact.
    ///
    /// The number of selections correspond to the distinct values returned by `func`.
    /// Selections are stored in a container `C` and has the same kind as subselections.
    pub fn split_disjoint<RT, F, C>(&self, func: F) -> C
    where
        RT: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(Particle) -> Option<RT>,
        C: FromIterator<Sel<K>> + Default,
    {
        self.split_gen(func)
    }

    /// Helper method that splits selection into the parts with distinct resids.
    /// Parent selection is left alive.
    pub fn split_disjoint_resid<C>(&self) -> C
    where
        C: FromIterator<Sel<K>> + Default,
    {
        self.split_gen(|p| Some(p.atom.resid))
    }

    /// Helper method that splits selection into the parts with distinct resindexes.
    /// Parent selection is left alive.
    pub fn split_disjoint_resindex<C>(&self) -> C
    where
        C: FromIterator<Sel<K>> + Default,
    {
        self.split_gen(|p| Some(p.atom.resindex))
    }

    /// Return serial iterator that splits selection into contigous pieces according to the value of function.
    /// Selection is not consumed. Returned selections are of subselection kind.
    ///
    /// Whenever `func` returns `Some(value)` different from the previous one, new selection is created.
    /// If `func` returns `None` the atom is skipped and do not added to new selection.
    /// Selections are computed lazily when iterating.
    pub fn iter_contig_fragments<RT, F>(&self, func: F) -> SelectionFragmentsIterator<'_, RT, F, K>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(Particle) -> Option<RT>,
    {
        SelectionFragmentsIterator::new(self, func)
    }

    /// Return serial iterator over contigous pieces of selection with distinct contigous resids.
    /// Parent selection is left alive.
    pub fn iter_contig_resindex(
        &self,
    ) -> SelectionFragmentsIterator<'_, usize, fn(Particle) -> Option<usize>, K> {
        self.iter_contig_fragments(|p| Some(p.atom.resindex))
    }

    //=========================================================
    // Splitting into parallel selections (keeps parent alive)
    //=========================================================

    /// Splits selection to pieces that could be processed in parallel.
    /// Pieces may not be contigous and are arranged in random order.
    /// Parent selection is kept intact.
    ///
    /// Each produced selection correspond to the distinct values returned by `split_fn`.
    /// Selections are stored in a special container [ParallelSplit].
    pub fn split_par_disjoint<F, RT>(&self, split_fn: F) -> Result<ParallelSplit, SelectionError>
    where
        F: Fn(Particle) -> Option<RT>,
        RT: Default + std::hash::Hash + std::cmp::Eq,
    {
        Ok(ParallelSplit {
            parts: self.split_gen(split_fn),
            _marker: Default::default(),
        })
    }

    fn split_par_contig_internal<F, RT>(
        &self,
        split_fn: F,
    ) -> Result<Vec<Sel<MutableParallel>>, SelectionError>
    where
        F: Fn(Particle) -> Option<RT>,
        RT: Default + std::hash::Hash + std::cmp::Eq,
    {
        // Iterator oved valid particles for which split_n returns Some
        let mut piter = self.iter_particle().filter_map(|p| {
            let i = p.id;
            split_fn(p).map(|val| (i, val))
        });

        // Get find first valid value or bail with error
        let (mut part_start, mut part_val) =
            piter.next().ok_or_else(|| SelectionError::EmptySplit)?;

        let mut parts = vec![];
        let mut last_valid = part_start;

        for (i, val) in piter {
            if val != part_val {
                // Create new part
                parts.push(unsafe {
                    Sel::from_holders_and_index(
                        self.topology.clone_into(),
                        self.state.clone_into(),
                        SortedSet::from_sorted((part_start..i).collect()),
                    )?
                });
                // Reset cur
                part_val = val;
                part_start = i;
            }
            // Keep track of the last valid value
            last_valid = i;
        }

        // Add last part
        parts.push(unsafe {
            Sel::from_holders_and_index(
                self.topology.clone_into(),
                self.state.clone_into(),
                SortedSet::from_sorted((part_start..=last_valid).collect()),
            )?
        });

        Ok(parts)
    }

    /// Splits selection to pieces that could be processed in parallel.
    /// Pieces may not be contigous and are arranged in random order.
    /// Parent selection is kept intact.
    ///
    /// Each produced selection correspond to the distinct values returned by `split_fn`.
    /// Selections are stored in a special container [ParallelSplit].
    pub fn split_par_contig<F, RT>(&self, split_fn: F) -> Result<ParallelSplit, SelectionError>
    where
        F: Fn(Particle) -> Option<RT>,
        RT: Default + std::hash::Hash + std::cmp::Eq,
    {
        Ok(ParallelSplit {
            parts: self.split_par_contig_internal(split_fn)?,
            _marker: Default::default(),
        })
    }

    pub fn par_iter_contig_fragments<F, RT>(
        &self,
        split_fn: F,
    ) -> Result<rayon::vec::IntoIter<Sel<MutableParallel>>, SelectionError>
    where
        F: Fn(Particle) -> Option<RT>,
        RT: Default + std::hash::Hash + std::cmp::Eq,
    {
        Ok(self.split_par_contig_internal(split_fn)?.into_par_iter())
    }

    //==============================================
    // Getters ans setters
    //==============================================

    pub fn get_topology(&self) -> Holder<Topology, K> {
        self.topology.clone()
    }

    pub fn get_state(&self) -> Holder<State, K> {
        self.state.clone()
    }

    pub fn set_topology(
        &mut self,
        topology: impl Into<Holder<Topology, K>>,
    ) -> Result<Holder<Topology, K>, SelectionError> {
        let topology = topology.into();
        if !self.topology.interchangeable(&topology) {
            return Err(SelectionError::SetState);
        }
        Ok(std::mem::replace(&mut self.topology, topology))
    }

    pub fn set_state(
        &mut self,
        state: impl Into<Holder<State, K>>,
    ) -> Result<Holder<State, K>, SelectionError> {
        let state = state.into();
        if !self.state.interchangeable(&state) {
            return Err(SelectionError::SetState);
        }
        Ok(std::mem::replace(&mut self.state, state))
    }

    pub fn set_same_name(&self, val: &str) {
        for a in self.topology.iter_atoms_mut() {
            a.name = val.into();
        }
    }

    pub fn set_same_resname(&self, val: &str) {
        for a in self.topology.iter_atoms_mut() {
            a.resname = val.into();
        }
    }

    pub fn set_same_resid(&self, val: i32) {
        for a in self.topology.iter_atoms_mut() {
            a.resid = val;
        }
    }

    pub fn set_same_chain(&self, val: char) {
        for a in self.topology.iter_atoms_mut() {
            a.chain = val;
        }
    }

    //==============================================
    // Direct creation of selections without Source
    //==============================================

    /// Creates new selection from an iterator of indexes. Indexes are bound checked, sorted and duplicates are removed.
    /// If any index is out of bounds the error is returned.
    pub fn from_iter(
        topology: &Holder<Topology, K>,
        state: &Holder<State, K>,
        iter: impl Iterator<Item = usize>,
    ) -> Result<Self, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        let vec = index_from_iter(iter, topology.num_atoms())?;
        Self::from_holders_and_index(
            topology.clone_with_index(&vec)?,
            state.clone_with_index(&vec)?,
            vec,
        )
    }

    /// Selects all
    pub fn from_all(
        topology: &Holder<Topology, K>,
        state: &Holder<State, K>,
    ) -> Result<Self, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        let vec = index_from_all(topology.num_atoms());
        Self::from_holders_and_index(
            topology.clone_with_index(&vec)?,
            state.clone_with_index(&vec)?,
            vec,
        )
    }

    /// Creates new selection from a selection expression string. Selection expression is constructed internally but
    /// can't be reused. Consider using [select_expr](Self::select_expr) if you already have selection expression.
    pub fn from_str(
        topology: &Holder<Topology, K>,
        state: &Holder<State, K>,
        selstr: &str,
    ) -> Result<Self, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        let vec = index_from_str(selstr, &topology, &state)?;
        Self::from_holders_and_index(
            topology.clone_with_index(&vec)?,
            state.clone_with_index(&vec)?,
            vec,
        )
    }

    /// Creates new selection from an existing selection expression.
    pub fn from_expr(
        topology: &Holder<Topology, K>,
        state: &Holder<State, K>,
        expr: &mut SelectionExpr,
    ) -> Result<Self, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        let vec = index_from_expr(expr, &topology, &state)?;
        Self::from_holders_and_index(
            topology.clone_with_index(&vec)?,
            state.clone_with_index(&vec)?,
            vec,
        )
    }

    /// Creates new selection from a range of indexes.
    /// If range is out of bounds the error is returned.
    pub fn from_range(
        topology: &Holder<Topology, K>,
        state: &Holder<State, K>,
        range: std::ops::Range<usize>,
    ) -> Result<Self, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        let vec = index_from_range(range, topology.num_atoms())?;
        Self::from_holders_and_index(
            topology.clone_with_index(&vec)?,
            state.clone_with_index(&vec)?,
            vec,
        )
    }

    //==============================================
    // Logic on selections that modify existing ones
    //==============================================
    pub fn append(&mut self, other: &impl IndexProvider) -> Result<(), SelectionError> {
        self.index_storage.extend(other.iter_index());
        self.check_index()?;
        Ok(())
    }

    pub fn exclude(&mut self, other: &impl IndexProvider) {
        let lhs = rustc_hash::FxHashSet::<usize>::from_iter(self.iter_index());
        let rhs = rustc_hash::FxHashSet::<usize>::from_iter(other.iter_index());
        let v: Vec<usize> = lhs.difference(&rhs).cloned().collect();
        self.index_storage = SortedSet::from_unsorted(v);
    }

    pub fn exclude_local(&mut self, other: impl IntoIterator<Item = usize>) {
        let lhs = rustc_hash::FxHashSet::<usize>::from_iter(0..self.len());
        let rhs = rustc_hash::FxHashSet::<usize>::from_iter(other);
        let v: Vec<usize> = lhs
            .difference(&rhs)
            .map(|i| self.index_storage[*i])
            .collect();
        self.index_storage = SortedSet::from_unsorted(v);
    }

    pub fn invert(&mut self) {
        let lhs = rustc_hash::FxHashSet::<usize>::from_iter(0..self.num_atoms());
        let rhs = rustc_hash::FxHashSet::<usize>::from_iter(self.iter_index());
        let v: Vec<usize> = lhs.difference(&rhs).cloned().collect();
        self.index_storage = SortedSet::from_unsorted(v);
    }

    //==============================================
    // Logic on selections that create new ones
    //==============================================

    pub fn union<KO: SelectionKind>(&self, other: &Sel<KO>) -> Self {
        let ind = union_sorted(self.index(), other.index());
        Self {
            topology: self.topology.clone(),
            state: self.state.clone(),
            index_storage: ind,
        }
    }

    pub fn intersection<KO: SelectionKind>(&self, other: &Sel<KO>) -> Result<Self, SelectionError> {
        let ind = intersection_sorted(self.index(), other.index());
        if ind.is_empty() {
            return Err(SelectionError::EmptyIntersection);
        }
        Ok(Self {
            topology: self.topology.clone(),
            state: self.state.clone(),
            index_storage: ind,
        })
    }

    pub fn difference<KO: SelectionKind>(&self, other: &Sel<KO>) -> Result<Self, SelectionError> {
        let ind = difference_sorted(self.index(), other.index());
        if ind.is_empty() {
            return Err(SelectionError::EmptyDifference);
        }
        Ok(Self {
            topology: self.topology.clone(),
            state: self.state.clone(),
            index_storage: ind,
        })
    }

    pub fn complement(&self) -> Result<Self, SelectionError> {
        let ind = difference_sorted(
            unsafe { &SortedSet::from_sorted((0..self.topology.num_atoms()).collect()) },
            self.index(),
        );
        if ind.is_empty() {
            return Err(SelectionError::EmptyComplement);
        }
        Ok(Self {
            topology: self.topology.clone(),
            state: self.state.clone(),
            index_storage: ind,
        })
    }
}

//---------------------------------------------

/// Iterator over the [Particle]s from selection.
pub struct SelectionIterator<'a, K: SelectionKind> {
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
pub struct SelectionIteratorMut<'a, K: SelectionKind> {
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
// IO traots
//---------------------------------------------

impl<K: SelectionKind> WritableToFile for Sel<K> {}

impl<K: SelectionKind> IndexProvider for Sel<K> {
    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize> {
        self.index().iter().cloned()
    }
}

impl<K: SelectionKind> TopologyProvider for Sel<K> {
    fn num_atoms(&self) -> usize {
        self.index().len()
    }
}

impl<K: SelectionKind> StateProvider for Sel<K> {
    fn get_time(&self) -> f32 {
        self.state.get_time()
    }

    fn num_coords(&self) -> usize {
        self.index().len()
    }
}

//==================================================================
// Immutable analysis traits
//==================================================================

impl<K: SelectionKind> BoxProvider for Sel<K> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.state.get_box()
    }
}

impl<K: SelectionKind> MeasurePeriodic for Sel<K> {}

impl<K: SelectionKind> PosProvider for Sel<K> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        unsafe {
            self.index()
                .iter()
                .map(|i| self.state.nth_pos_unchecked(*i))
        }
    }
}

impl<K: SelectionKind> MeasurePos for Sel<K> {}

impl<K: SelectionKind> AtomsProvider for Sel<K> {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        unsafe {
            self.index()
                .iter()
                .map(|i| self.topology.nth_atom_unchecked(*i))
        }
    }
}

impl<K: SelectionKind> MassesProvider for Sel<K> {
    fn iter_masses(&self) -> impl ExactSizeIterator<Item = f32> {
        unsafe {
            self.index()
                .iter()
                .map(|i| self.topology.nth_atom_unchecked(*i).mass)
        }
    }
}

impl<K: SelectionKind> MeasureMasses for Sel<K> {}

impl<K: SelectionKind> LenProvider for Sel<K> {
    fn len(&self) -> usize {
        self.index().len()
    }
}

impl<K: SelectionKind> RandomPos for Sel<K> {
    unsafe fn nth_pos_unchecked(&self, i: usize) -> &Pos {
        let ind = *self.index().get_unchecked(i);
        self.state.nth_pos_unchecked(ind)
    }
}

impl<K: SelectionKind> RandomAtom for Sel<K> {
    fn nth_atom(&self, i: usize) -> Option<&Atom> {
        self.index()
            .get(i)
            .map(|i| unsafe { self.topology.nth_atom_unchecked(*i) })
    }

    unsafe fn nth_atom_unchecked(&self, i: usize) -> &Atom {
        let ind = *self.index().get_unchecked(i);
        self.topology.nth_atom_unchecked(ind)
    }
}

//-------------------------------------------------------
// Mutable analysis traits only for mutable selections
//-------------------------------------------------------

impl<K: MutableSel> PosMutProvider for Sel<K> {
    fn iter_pos_mut(&self) -> impl PosMutIterator<'_> {
        unsafe {
            self.index()
                .iter()
                .map(|i| self.state.nth_pos_unchecked_mut(*i))
        }
    }
}

impl<K: MutableSel> AtomsMutProvider for Sel<K> {
    fn iter_atoms_mut(&self) -> impl AtomMutIterator<'_> {
        unsafe {
            self.index()
                .iter()
                .map(|i| self.topology.nth_atom_unchecked_mut(*i))
        }
    }
}

impl<K: MutableSel> RandomPosMut for Sel<K> {
    fn nth_pos_mut(&self, i: usize) -> Option<&mut Pos> {
        self.index()
            .get(i)
            .map(|i| unsafe { self.state.nth_pos_mut_unchecked(*i) })
    }

    unsafe fn nth_pos_mut_unchecked(&self, i: usize) -> &mut Pos {
        let ind = *self.index().get_unchecked(i);
        self.state.nth_pos_unchecked_mut(ind)
    }
}

impl<K: MutableSel> RandomAtomMut for Sel<K> {
    fn nth_atom_mut(&self, i: usize) -> Option<&mut Atom> {
        self.index()
            .get(i)
            .map(|i| unsafe { self.topology.nth_atom_mut_unchecked(*i) })
    }

    unsafe fn nth_atom_mut_unchecked(&self, i: usize) -> &mut Atom {
        let ind = *self.index().get_unchecked(i);
        self.topology.nth_atom_mut_unchecked(ind)
    }
}

impl<K: MutableSel> ModifyPos for Sel<K> {}
impl<K: MutableSel> ModifyPeriodic for Sel<K> {}
impl<K: MutableSel> ModifyRandomAccess for Sel<K> {}

//----------------------------------------------
