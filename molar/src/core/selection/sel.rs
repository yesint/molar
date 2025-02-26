use crate::prelude::*;
use sorted_vec::SortedSet;
use std::collections::HashMap;

use super::{next_split, utils::*, SplitData};

//═══════════════════
//███  Selection
//═══════════════════

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
    topology: Holder<Topology, K>,
    state: Holder<State, K>,
    index_storage: SortedSet<usize>,
}

//══════════════════════════════════════════════
//███  Functions shared by all selection kinds
//══════════════════════════════════════════════

//━━━━━━━━━━━━━━━━━━━━━━━━━━
//       Public API
//━━━━━━━━━━━━━━━━━━━━━━━━━━

impl<K: SelectionKind> Sel<K> {
    #[inline(always)]
    pub(super) fn index(&self) -> &SortedSet<usize> {
        self.check_index().unwrap();
        &self.index_storage
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        // We need to check index here manually
        self.check_index().unwrap();
        self.num_atoms()
    }

    pub fn nth_pos(&self, i: usize) -> Option<&Pos> {
        let ind = *self.index().get(i)?;
        Some(unsafe { self.state.nth_pos_unchecked(ind) })
    }

    pub fn nth_pos_mut(&self, i: usize) -> Option<&mut Pos> {
        let ind = *self.index().get(i)?;
        Some(unsafe { self.state.nth_pos_unchecked_mut(ind) })
    }

    //=================================
    // Getters, Setters and Accessors
    //=================================

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

    /// Get an [Atom] for the first selection index.
    pub fn first_atom(&self) -> &Atom {
        unsafe { self.nth_atom_unchecked(0) }
    }

    /// Get an [Atom] for the last selection index.
    pub fn last_atom(&self) -> &Atom {
        unsafe { self.nth_atom_unchecked(self.index().len() - 1) }
    }

    /// Get an [Atom] for the n-th selection index.
    /// Index is bound-checked, an error is returned if it is out of bounds.
    pub fn nth_atom(&self, i: usize) -> Option<&Atom> {
        if i >= self.len() {
            None
        } else {
            Some(unsafe { self.nth_atom_unchecked(i) })
        }
    }

    pub fn get_topology(&self) -> Holder<Topology, K> {
        self.topology.clone_with_kind()
    }

    pub fn get_state(&self) -> Holder<State, K> {
        self.state.clone_with_kind()
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

}

//──────────────────────────
//      Private API
//──────────────────────────

impl<K: SelectionKind> Sel<K> {
    #[inline(always)]
    fn check_index(&self) -> Result<(), SelectionError> {
        K::check_index(&self.index_storage, &self.topology, &self.state)
    }

    // Get a Particle for i-th selection index.
    // # Safety
    // This is an unsafe operation that doesn't check if the index is in bounds.
    pub(super) unsafe fn nth_particle_unchecked(&self, i: usize) -> Particle<'_> {
        let ind = *self.index_storage.get_unchecked(i);
        Particle {
            id: ind,
            atom: self.topology.nth_atom_unchecked(ind),
            pos: self.state.nth_pos_unchecked(ind),
        }
    }

    // Get a mutable Particle for i-th selection index with
    // # Safety
    // This is an unsafe operation that doesn't check if the index is in bounds.
    pub(super) unsafe fn nth_particle_unchecked_mut(&self, i: usize) -> ParticleMut<'_> {
        let ind = *self.index_storage.get_unchecked(i);
        ParticleMut {
            id: ind,
            atom: self.topology.nth_atom_unchecked_mut(ind),
            pos: self.state.nth_pos_unchecked_mut(ind),
        }
    }

    // This method doesn't check if the vector is sorted and has duplicates and thus unsafe
    pub(super) unsafe fn subsel_from_sorted_vec_unchecked<KO: SelectionKind>(
        &self,
        index: Vec<usize>,
    ) -> Result<Sel<KO>, SelectionError> {
        if index.len() > 0 {
            Ok(Sel {
                topology: self.topology.clone_with_kind(),
                state: self.state.clone_with_kind(),
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
                topology: self.topology.clone_with_kind(),
                state: self.state.clone_with_kind(),
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

    pub(super) fn from_holders_and_index(
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
}

//══════════════════════════════════════════════
//███  Functions shared by user creatable kinds
//══════════════════════════════════════════════

//──────────────────────────
//      Private API
//──────────────────────────

impl<K: UserCreatableKind> Sel<K> {
    // Helper splitting function doing actual work
    fn split_collect_internal<RT, F, C, KO>(&self, func: F) -> Result<C,SelectionError>
    where
        RT: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(Particle) -> Option<RT>,
        C: FromIterator<Sel<KO>> + Default,
        KO: SelectionKind,
    {
        self.check_index()?;
        let mut values = HashMap::<RT, Vec<usize>>::default();

        for p in self.iter_particle() {
            let i = p.id;
            if let Some(val) = func(p) {
                if let Some(el) = values.get_mut(&val) {
                    el.push(i);
                } else {
                    values.insert(val, vec![i]);
                }
            }
        }

        Ok(C::from_iter(
            values.into_values()
                .map(|ind| unsafe { self.subsel_from_unsorted_vec_unchecked(ind).unwrap() }),
            // This should never fail because `ind` can't be empty and can't contain duplicates
        ))
    }

    fn split_par_contig_internal<F, RT>(
        &self,
        split_fn: F,
    ) -> Result<Vec<Sel<MutableParallel>>, SelectionError>
    where
        F: Fn(Particle) -> Option<RT>,
        RT: Default + std::hash::Hash + std::cmp::Eq,
    {
        let mut data = SplitData {
            func: split_fn,
            cur: 0,
            val: RT::default(),
        };

        let mut parts = vec![];
        while let Some(part) = next_split(&mut data, self) {
            parts.push(part);
        }
        Ok(parts)
    }
}

//━━━━━━━━━━━━━━━━━━━━━━━━━━
//       Public API
//━━━━━━━━━━━━━━━━━━━━━━━━━━

impl<K: UserCreatableKind> Sel<K> {
    //============================
    // Splitting
    //============================

    // Naming convention is
    // split[_into][_<property>][_[par_]iter]

    //----------------------
    // Consuming splitters
    //----------------------

    /// Splits selection to pieces that could be disjoint
    /// according to the value of function. Parent selection is consumed.    
    ///
    /// The number of selections correspond to the distinct values returned by `func`.
    /// Selections are stored in a container `C` and has the same kind as parent selection.
    pub fn split_into<RT, F, C>(self, func: F) -> Result<C,SelectionError>
    where
        RT: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(Particle) -> Option<RT>,
        C: FromIterator<Sel<K>> + Default,
    {
        self.split_collect_internal(func)
    }

    /// Helper method that splits selection into the parts with distinct resids.
    /// Parent selection is consumed.
    pub fn split_resid_into<C>(self) -> Result<C,SelectionError>
    where
        C: FromIterator<Sel<K>> + Default,
    {
        self.split_collect_internal(|p| Some(p.atom.resid))
    }

    /// Return iterator that splits selection into contigous pieces according to the value of function.
    /// Consumes a selection and returns selections of the same kind.
    ///
    /// Whenever `func` returns a value different from the previous one, new selection is created.
    /// Selections are computed lazily when iterating.
    pub fn split_into_iter<RT, F>(self, func: F) -> IntoFragmentsIterator<RT, F, K>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(Particle) -> Option<RT>,
    {
        IntoFragmentsIterator::new(self, func)
    }

    /// Return iterator over contigous pieces of selection with distinct contigous resids.
    /// Parent selection is consumed.
    pub fn split_resindex_into_iter(
        self,
    ) -> IntoFragmentsIterator<usize, impl Fn(Particle) -> Option<usize>, K> {
        self.split_into_iter(|p| Some(p.atom.resindex))
    }

    pub fn split_chain_into_iter(
        self,
    ) -> IntoFragmentsIterator<char, impl Fn(Particle) -> Option<char>, K> {
        self.split_into_iter(|p| Some(p.atom.chain))
    }

    //--------------------------
    // Non-consuming splitters
    //-------------------------

    /// Splits selection to pieces that could be disjoint
    /// according to the value of function. Parent selection is kept intact.
    ///
    /// The number of selections correspond to the distinct values returned by `func`.
    /// Selections are stored in a container `C` and has the same kind as subselections.
    pub fn split<RT, F, C>(&self, func: F) -> Result<C,SelectionError>
    where
        RT: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(Particle) -> Option<RT>,
        C: FromIterator<Sel<K>> + Default,
    {
        self.split_collect_internal(func)
    }

    /// Helper method that splits selection into the parts with distinct resids.
    /// Parent selection is left alive.
    pub fn split_resid<C>(&self) -> Result<C,SelectionError>
    where
        C: FromIterator<Sel<K>> + Default,
    {
        self.split_collect_internal(|p| Some(p.atom.resid))
    }

    /// Helper method that splits selection into the parts with distinct resindexes.
    /// Parent selection is left alive.
    pub fn split_resindex<C>(&self) -> Result<C,SelectionError>
    where
        C: FromIterator<Sel<K>> + Default,
    {
        self.split_collect_internal(|p| Some(p.atom.resindex))
    }

    /// Return serial iterator that splits selection into contigous pieces according to the value of function.
    /// Selection is not consumed.
    ///
    /// Whenever `func` returns `Some(value)` different from the previous one, new selection is created.
    /// If `func` returns `None` the atom is skipped and do not added to new selection.
    /// Selections are computed lazily when iterating.
    pub fn split_iter<'a, RT, F>(&'a self, func: F) -> Result<FragmentsIterator<'a, RT, F, K>,SelectionError>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(Particle) -> Option<RT> + 'a,
    {
        self.check_index()?;
        Ok(FragmentsIterator::new(self, func))
    }

    /// Return serial iterator over contigous pieces of selection with distinct contigous resids.
    /// Parent selection is left alive.
    pub fn split_resindex_iter<'a>(&'a self) -> Result<FragmentsIterator<'a, usize, impl Fn(Particle) -> Option<usize> + 'a , K>,SelectionError> {
        self.split_iter(|p| Some(p.atom.resindex))
    }

    //---------------------------------------------------------
    // Splitting into parallel selections (non-consuming)
    //---------------------------------------------------------

    /// Splits selection to pieces that could be processed in parallel.
    /// Pieces may not be contigous and are arranged in random order.
    /// Parent selection is kept intact.
    ///
    /// Each produced selection correspond to the distinct values returned by `split_fn`.
    /// Selections are stored in a special container [ParallelSplit].
    pub fn split_par_disjoint<F, RT>(&self, split_fn: F) -> Result<ParallelSplit,SelectionError>
    where
        F: Fn(Particle) -> Option<RT>,
        RT: Default + std::hash::Hash + std::cmp::Eq,
    {
        Ok(ParallelSplit {
            parts: self.split_collect_internal(split_fn)?,
            _marker: Default::default(),
        })
    }

    /// Splits selection to pieces that could be processed in parallel.
    /// Pieces are contigous and are arranged in order of appearance.
    /// Parent selection is kept intact.
    ///
    /// New selection starts when `split_fn` returns a value different from the previous one.
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

    pub fn split_molecules_iter(&self) -> Result<MoleculesIterator<K>,SelectionError> {
        Ok(MoleculesIterator {
            sel: self,
            cur: 0,
        })
    }

    //============================
    // Misceleneous
    //============================

    /// Tests if two selections are from the same source
    pub fn is_same_source(&self, other: &Sel<K>) -> bool {
        self.topology.same_data(&other.topology) && self.state.same_data(&other.state)
    }

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

    //===================
    // Subselections
    //===================

    fn subselect_internal(&self, index: SortedSet<usize>) -> Result<Self, SelectionError> {
        Self::from_holders_and_index(self.topology.clone_with_kind(), self.state.clone_with_kind(), index)
    }

    /// Subselection from expression
    pub fn subsel_expr(&self, expr: &mut SelectionExpr) -> Result<Sel<K>, SelectionError> {
        let index = index_from_expr_sub(expr, &self.topology, &self.state, &self.index())?;
        Self::from_holders_and_index(self.topology.clone_with_kind(), self.state.clone_with_kind(), index)
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
        Self::from_holders_and_index(topology.clone_with_kind(), state.clone_with_kind(), vec)
    }

    /// Selects all
    pub fn from_all(
        topology: &Holder<Topology, K>,
        state: &Holder<State, K>,
    ) -> Result<Self, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        let vec = index_from_all(topology.num_atoms());
        Self::from_holders_and_index(topology.clone_with_kind(), state.clone_with_kind(), vec)
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
        Self::from_holders_and_index(topology.clone_with_kind(), state.clone_with_kind(), vec)
    }

    /// Creates new selection from an existing selection expression.
    pub fn from_expr(
        topology: &Holder<Topology, K>,
        state: &Holder<State, K>,
        expr: &mut SelectionExpr,
    ) -> Result<Self, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        let vec = index_from_expr(expr, &topology, &state)?;
        Self::from_holders_and_index(topology.clone_with_kind(), state.clone_with_kind(), vec)
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
        Self::from_holders_and_index(topology.clone_with_kind(), state.clone_with_kind(), vec)
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

    pub fn union<KO: UserCreatableKind>(&self, other: &Sel<KO>) -> Self {
        let ind = union_sorted(self.index(), other.index());
        Self {
            topology: self.topology.clone_with_kind(),
            state: self.state.clone_with_kind(),
            index_storage: ind,
        }
    }

    pub fn intersection<KO: UserCreatableKind>(
        &self,
        other: &Sel<KO>,
    ) -> Result<Self, SelectionError> {
        let ind = intersection_sorted(self.index(), other.index());
        if ind.is_empty() {
            return Err(SelectionError::EmptyIntersection);
        }
        Ok(Self {
            topology: self.topology.clone_with_kind(),
            state: self.state.clone_with_kind(),
            index_storage: ind,
        })
    }

    pub fn difference<KO: UserCreatableKind>(
        &self,
        other: &Sel<KO>,
    ) -> Result<Self, SelectionError> {
        let ind = difference_sorted(self.index(), other.index());
        if ind.is_empty() {
            return Err(SelectionError::EmptyDifference);
        }
        Ok(Self {
            topology: self.topology.clone_with_kind(),
            state: self.state.clone_with_kind(),
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
            topology: self.topology.clone_with_kind(),
            state: self.state.clone_with_kind(),
            index_storage: ind,
        })
    }
}

//══════════════════════════════════════════════
//███  Iterator over Particles
//══════════════════════════════════════════════

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

//══════════════════════════════════════════════
//███  Mutable iterator over Particles
//══════════════════════════════════════════════

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

//══════════════════════════════════════════════
//███  IO traits
//══════════════════════════════════════════════

impl<K: SelectionKind> WritableToFile for Sel<K> {}

impl<K: SelectionKind> IndexProvider for Sel<K> {
    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize> {
        self.index().iter().cloned()
    }
}

impl<K: SelectionKind> TopologyProvider for Sel<K> {}

impl<K: SelectionKind> StateProvider for Sel<K> {
    fn get_time(&self) -> f32 {
        self.state.get_time()
    }
}

//══════════════════════════════════════════════
//███  Immutable analysis traits
//══════════════════════════════════════════════

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

impl<K: SelectionKind> AtomProvider for Sel<K> {
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

impl<K: SelectionKind> RandomPosProvider for Sel<K> {
    fn num_coords(&self) -> usize {
        self.index().len()    
    }

    unsafe fn nth_pos_unchecked(&self, i: usize) -> &Pos {
        let ind = *self.index().get_unchecked(i);
        self.state.nth_pos_unchecked(ind)
    }
}

impl<K: SelectionKind> RandomAtomProvider for Sel<K> {
    fn num_atoms(&self) -> usize {
        self.index().len()
    }

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

impl<K: SelectionKind> MeasureRandomAccess for Sel<K> {}

impl<K: SelectionKind> MoleculesProvider for Sel<K> {
    fn num_molecules(&self) -> usize {
        self.topology.num_molecules()
    }

    fn iter_molecules(&self) -> impl Iterator<Item = &[usize; 2]> {
        self.topology.iter_molecules()
    }

    unsafe fn nth_molecule_unchecked(&self, i: usize) -> &[usize;2] {
        self.topology.nth_molecule_unchecked(i)
    }
}

impl<K: SelectionKind> BondsProvider for Sel<K> {
    fn num_bonds(&self) -> usize {
        self.topology.num_bonds()
    }

    fn iter_bonds(&self) -> impl Iterator<Item = &[usize; 2]> {
        self.topology.iter_bonds()
    }

    unsafe fn nth_bond_unchecked(&self, i: usize) -> &[usize;2] {
        self.topology.nth_bond_unchecked(i)
    }
}

//═══════════════════════════════════════════════════════════
//███  Mutable analysis traits (only for mutable selections)
//═══════════════════════════════════════════════════════════

impl<K: MutableKind> PosMutProvider for Sel<K> {
    fn iter_pos_mut(&self) -> impl PosMutIterator<'_> {
        unsafe {
            self.index()
                .iter()
                .map(|i| self.state.nth_pos_unchecked_mut(*i))
        }
    }
}

impl<K: MutableKind> AtomsMutProvider for Sel<K> {
    fn iter_atoms_mut(&self) -> impl AtomMutIterator<'_> {
        unsafe {
            self.index()
                .iter()
                .map(|i| self.topology.nth_atom_unchecked_mut(*i))
        }
    }
}

impl<K: MutableKind> RandomPosMut for Sel<K> {
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

impl<K: MutableKind> RandomAtomMutProvider for Sel<K> {
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

impl<K: MutableKind> ModifyPos for Sel<K> {}
impl<K: MutableKind> ModifyPeriodic for Sel<K> {}
impl<K: MutableKind> ModifyRandomAccess for Sel<K> {}
