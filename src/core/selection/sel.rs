use crate::prelude::*;
use sorted_vec::SortedSet;
use std::{collections::HashMap, marker::PhantomData};
use triomphe::{Arc, UniqueArc};

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
pub struct Sel<K> {
    topology: triomphe::Arc<Topology>,
    state: triomphe::Arc<State>,
    index_storage: SortedSet<usize>,
    _marker: PhantomData<K>,
}

//-------------------------------------------
impl<K: SelectionKind> Sel<K> {
    #[inline(always)]
    fn index(&self) -> &SortedSet<usize> {
        K::check_index(&self.index_storage, &self.topology, &self.state).unwrap();
        &self.index_storage
    }

    // Only visible in selection module
    pub(super) fn new_internal(
        topology: Arc<Topology>,
        state: Arc<State>,
        index: SortedSet<usize>,
    ) -> Self {
        Self {
            topology,
            state,
            index_storage: index,
            _marker: PhantomData::default(),
        }
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        // We need to check index here manually
        K::check_index(&self.index_storage, &self.topology, &self.state).unwrap();
        self.num_atoms()
    }

    // This method doesn't check if the vector has duplicates and thus unsafe
    pub(super) unsafe fn subsel_from_vec_unchecked<S: SelectionKind>(
        &self,
        index: Vec<usize>,
    ) -> Result<Sel<S>, SelectionError> {
        if index.len() > 0 {
            Ok(Sel::new_internal(
                Arc::clone(&self.topology),
                Arc::clone(&self.state),
                unsafe { SortedSet::from_sorted(index) },
            ))
        } else {
            Err(SelectionError::FromVec {
                first: index[0],
                last: index[index.len() - 1],
                size: index.len(),
                source: SelectionIndexError::IndexEmpty,
            })
        }
    }

    pub fn subsel_from_vec<S: SelectionKind>(
        &self,
        index: Vec<usize>,
    ) -> Result<Sel<S>, SelectionError> {
        Ok(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            index_from_vec(index, self.len())?,
        ))
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

    // Helper splitting function generic over returned selections kind
    fn split_gen<RT, F, C, Kind>(&self, func: F) -> C
    where
        RT: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(Particle) -> Option<RT>,
        C: FromIterator<Sel<Kind>> + Default,
        Kind: SelectionKind,
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
                .map(|ind| unsafe { self.subsel_from_vec_unchecked(ind).unwrap() }),
            // This should never fail because `ind` can't be empty
        )
    }

    /// Splits selection to pieces that could be disjoint
    /// according to the value of function. Parent selection is consumed.
    ///
    /// The number of selections correspond to the distinct values returned by `func`.
    /// Selections are stored in a container `C` and has the same kind as parent selection.
    pub fn split_into<RT, F, C>(self, func: F) -> C
    where
        RT: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(Particle) -> Option<RT>,
        C: FromIterator<Sel<K>> + Default,
    {
        self.split_gen(func)
    }

    /// Helper method that splits selection into the parts with distinct resids.
    /// Parent selection is consumed.
    pub fn split_resid_into<C>(self) -> C
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

    /// Saves selection to file. File type is deduced from extension.
    /*
    pub fn save(&self, fname: &str) -> Result<()> {
        let mut h = FileHandler::create(fname)?;
        h.write(self)
    }
    */

    pub fn first_index(&self) -> usize {
        self.index()[0]
    }

    pub fn last_index(&self) -> usize {
        self.index()[self.index().len() - 1]
    }

    /// Get a Particle for the first selection index.
    pub fn first_particle(&self) -> Particle {
        unsafe { self.nth_particle_unchecked(0) }
    }

    /// Get a Particle for the last selection index.
    pub fn last_particle(&self) -> Particle {
        unsafe { self.nth_particle_unchecked(self.index_storage.len() - 1) }
    }

    /// Get a Particle for the first selection index.
    /// Index is bound-checked, an error is returned if it is out of bounds.
    pub fn nth_particle(&self, i: usize) -> Result<Particle, SelectionError> {
        if i >= self.len() {
            Err(SelectionError::OutOfBounds(i, self.len()))
        } else {
            Ok(unsafe { self.nth_particle_unchecked(i) })
        }
    }

    /// Return iterator that splits selection into contigous pieces according to the value of function.
    /// Consumes a selection and returns selections of the same kind.
    ///
    /// Whenever `func` returns a value different from the previous one, new selection is created.
    /// Selections are computed lazily when iterating.
    pub fn into_split_contig<RT, F>(self, func: F) -> IntoSelectionSplitIterator<RT, F, K>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(Particle) -> Option<RT>,
    {
        IntoSelectionSplitIterator::new(self, func)
    }

    /// Return iterator over contigous pieces of selection with distinct contigous resids.
    /// Parent selection is consumed.
    pub fn into_split_contig_resindex(
        self,
    ) -> IntoSelectionSplitIterator<usize, impl Fn(Particle) -> Option<usize>, K> {
        self.into_split_contig(|p| Some(p.atom.resindex))
    }

    /// Return iterator that splits selection into contigous pieces according to the value of function.
    /// Selection is not consumed. Returned selections are of subselection kind.
    ///
    /// Whenever `func` returns `Some(value)` different from the previous one, new selection is created.
    /// If `func` returns `None` the atom is skipped and do not added to new selection.
    /// Selections are computed lazily when iterating.
    pub fn split_contig<RT, F>(&self, func: F) -> SelectionSplitIterator<'_, RT, F, K>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(Particle) -> Option<RT>,
    {
        SelectionSplitIterator::new(self, func)
    }

    /// Return iterator over contigous pieces of selection with distinct contigous resids.
    /// Parent selection is left alive.
    pub fn split_contig_resindex(
        &self,
    ) -> SelectionSplitIterator<'_, usize, fn(Particle) -> Option<usize>, K> {
        self.split_contig(|p| Some(p.atom.resindex))
    }

    /// Tests if two selections are from the same source
    pub fn same_source(&self, other: &Sel<K>) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology) && Arc::ptr_eq(&self.state, &other.state)
    }
}

impl<K: AllowsSubselect> Sel<K> {
    //===================
    // Subselections
    //===================

    /// Subselection from expression
    pub fn subsel_from_expr(&self, expr: &SelectionExpr) -> Result<Sel<K>, SelectionError> {
        Ok(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            index_from_expr_sub(expr, &self.topology, &self.state, &self.index())?,
        ))
    }

    /// Subselection from string
    pub fn subsel_from_str(&self, sel_str: &str) -> Result<Sel<K>, SelectionError> {
        let expr = SelectionExpr::try_from(sel_str)?;
        self.subsel_from_expr(&expr)
    }

    /// Subselection from the range of local selection indexes
    pub fn subsel_from_local_range(
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

        Ok(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            unsafe { SortedSet::from_sorted(index) },
        ))
    }

    /// Subselection from iterator that provides local selection indexes
    pub fn subsel_from_iter(
        &self,
        iter: impl ExactSizeIterator<Item = usize>,
    ) -> Result<Sel<K>, SelectionError> {
        let ind = self.index();
        let global_vec: Vec<_> = iter.map(|i| ind[i]).collect();
        Ok(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            index_from_vec(global_vec, self.len())?,
        ))
    }

    //==============================================
    // Splitting that keeps parent selection alive
    //==============================================

    /// Splits selection to pieces that could be disjoint
    /// according to the value of function. Parent selection is kept intact.
    ///
    /// The number of selections correspond to the distinct values returned by `func`.
    /// Selections are stored in a container `C` and has the same kind as subselections.
    pub fn split<RT, F, C>(&self, func: F) -> C
    where
        RT: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(Particle) -> Option<RT>,
        C: FromIterator<Sel<K>> + Default,
    {
        self.split_gen(func)
    }

    /// Helper method that splits selection into the parts with distinct resids.
    /// Parent selection is left alive.
    pub fn split_resid<C>(&self) -> C
    where
        C: FromIterator<Sel<K>> + Default,
    {
        self.split_gen(|p| Some(p.atom.resid))
    }

    /// Helper method that splits selection into the parts with distinct resindexes.
    /// Parent selection is left alive.
    pub fn split_resindex<C>(&self) -> C
    where
        C: FromIterator<Sel<K>> + Default,
    {
        self.split_gen(|p| Some(p.atom.resindex))
    }
}

impl<K: SerialSel> Sel<K> {
    pub fn get_shared_topology(&self) -> Holder<Topology, K> {
        Holder::from_arc(Arc::clone(&self.topology))
    }

    pub fn get_shared_state(&self) -> Holder<State, K> {
        Holder::from_arc(Arc::clone(&self.state))
    }

    pub fn set_shared_topology(
        &mut self,
        topology: Holder<Topology, K>,
    ) -> Result<Holder<Topology, K>, SelectionError> {
        if !self.topology.interchangeable(&topology) {
            return Err(SelectionError::SetTopology);
        }
        Ok(Holder::from_arc(std::mem::replace(
            &mut self.topology,
            topology.into_arc(),
        )))
    }

    pub fn set_shared_state(
        &mut self,
        state: Holder<State, K>,
    ) -> Result<Holder<State, K>, SelectionError> {
        if !self.state.interchangeable(&state) {
            return Err(SelectionError::SetState);
        }
        Ok(Holder::from_arc(std::mem::replace(
            &mut self.state,
            state.into_arc(),
        )))
    }

    pub fn set_state(&mut self, state: UniqueArc<State>) -> Result<Arc<State>, SelectionError> {
        if !self.state.interchangeable(&state) {
            return Err(SelectionError::SetState);
        }
        Ok(std::mem::replace(&mut self.state, state.shareable()))
    }

    //---------------------------------------------------------------
    // For serial selections direct creation is available
    //---------------------------------------------------------------

    /// Creates new selection from an iterator of indexes. Indexes are bound checked, sorted and duplicates are removed.
    /// If any index is out of bounds the error is returned.
    pub fn from_iter(
        topology: &Arc<Topology>,
        state: &Arc<State>,
        iter: impl Iterator<Item = usize>,
    ) -> Result<Self, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        let vec = index_from_iter(iter, topology.num_atoms())?;
        Ok(Sel::new_internal(
            Arc::clone(&topology),
            Arc::clone(&state),
            vec,
        ))
    }

    /// Selects all
    pub fn all(topology: &Arc<Topology>, state: &Arc<State>) -> Result<Self, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        let vec = index_from_all(topology.num_atoms());
        Ok(Sel::new_internal(
            Arc::clone(&topology),
            Arc::clone(&state),
            vec,
        ))
    }

    /// Creates new selection from a selection expression string. Selection expression is constructed internally but
    /// can't be reused. Consider using [select_expr](Self::select_expr) if you already have selection expression.
    pub fn from_str(
        topology: &Arc<Topology>,
        state: &Arc<State>,
        selstr: &str,
    ) -> Result<Self, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        let vec = index_from_str(selstr, &topology, &state)?;
        Ok(Sel::new_internal(
            Arc::clone(&topology),
            Arc::clone(&state),
            vec,
        ))
    }

    /// Creates new selection from an existing selection expression.
    pub fn from_expr(
        topology: &Arc<Topology>,
        state: &Arc<State>,
        expr: &SelectionExpr,
    ) -> Result<Self, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        let vec = index_from_expr(expr, &topology, &state)?;
        Ok(Sel::new_internal(
            Arc::clone(&topology),
            Arc::clone(&state),
            vec,
        ))
    }

    /// Creates new selection from a range of indexes.
    /// If rangeis out of bounds the error is returned.
    pub fn from_range(
        topology: &Arc<Topology>,
        state: &Arc<State>,
        range: &std::ops::Range<usize>,
    ) -> Result<Self, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        let vec = index_from_range(range, topology.num_atoms())?;
        Ok(Sel::new_internal(
            Arc::clone(&topology),
            Arc::clone(&state),
            vec,
        ))
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
// Implement analysis traits

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
