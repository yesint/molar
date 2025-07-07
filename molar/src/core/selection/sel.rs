use crate::prelude::*;
use itertools::Itertools;
use molar_powersasa::SasaResults;
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
/// There are four kinds of selections set by the generic marker parameter:
/// ### Serial mutable ([MutableSerial])
/// * May overlap.
/// * Mutable.
/// * Could only be accessed from the same thread where they were created (they are neither [Send] nor [Sync]).
/// ### Builder selections ([BuilderSerial])
/// * May overlap.
/// * Mutable.
/// * Could only be accessed from the same thread where they were created (they are neither [Send] nor [Sync]).
/// * Should be used together with `[Source::<BuilderSerial>]` to add or delete atoms from the system 
/// (other selection kinds don't allow this). 
/// <section class="warning">
/// This selection kind involves additional checks to ensure that selection
/// indexes are valid after potential deletion of atoms, so marginally slower than other selection kinds.
/// </section>
/// ### Parallel immutable ([ImmutableParallel])
/// * May overlap.
/// * Immutable.
/// * Could be processed in parallel from multiple threads and sent between threads.
/// ### Parallel mutable ([MutableParallel])
/// * _Can't_ overlap
/// * Mutable
/// * Could be processed in parallel from multiple threads and sent between threads.
/// * Can't be created directly. Obtained by splitting MutableSerial selections only.
///
/// Immutable selections implement the traits that allow to query properties, while
/// mutable also implement traits that allow to modify atoms and coordinates.
///
/// # Subselections
/// It is possible to select within existing selection using `subsel` methods.

#[derive(Debug)]
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
        self.validate_index().unwrap();
        &self.index_storage
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        // We need to validate index here manually
        self.validate_index().unwrap();
        self.num_atoms()
    }

    //=================================
    // Getters, Setters and Accessors
    //=================================

    /// Returns reference to index vector.
    pub fn get_index_vec(&self) -> &SortedSet<usize> {
        &self.index()
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
    
    /// Get a Particle for the last selection index.
    pub fn last_particle(&self) -> Particle {
        unsafe { self.nth_particle_unchecked(self.index().len() - 1) }
    }

    pub fn get_topology(&self) -> Holder<Topology, K> {
        unsafe {self.topology.new_ref_with_kind()}
    }

    pub fn get_state(&self) -> Holder<State, K> {
        unsafe {self.state.new_ref_with_kind()}
    }

    pub fn set_topology(
        &mut self,
        topology: impl Into<Holder<Topology, K>>,
    ) -> Result<Holder<Topology, K>, SelectionError> {
        let topology = topology.into();
        if !self.topology.interchangeable(&topology) {
            return Err(SelectionError::IncompatibleState);
        }
        Ok(std::mem::replace(&mut self.topology, topology))
    }

    /// Sets new [State] in [Sel].
    /// This is "shallow" operation, which doesn't affect other
    /// [Sel]s that point to the same [State].
    pub fn set_state(
        &mut self,
        state: impl Into<Holder<State, K>>,
    ) -> Result<Holder<State, K>, SelectionError> {
        let state: Holder<State,K> = state.into();
        if !self.state.interchangeable(&state) {
            return Err(SelectionError::IncompatibleState);
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

    pub fn set_same_mass(&self, val: f32) {
        for a in self.topology.iter_atoms_mut() {
            a.mass = val;
        }
    }
}

//──────────────────────────
//      Private API
//──────────────────────────

impl<K: SelectionKind> Sel<K> {
    pub(super) unsafe fn into_other_kind<KO: UserCreatableKind>(self) -> Sel<KO> {
        Sel {
            topology: self.topology.new_ref_with_kind(),
            state: self.state.new_ref_with_kind(),
            index_storage: self.index_storage,
        }
    }

    #[inline(always)]
    fn validate_index(&self) -> Result<(), SelectionError> {
        K::validate_index(&self.index_storage, &self.topology, &self.state)
    }

    // This method doesn't check if the vector is sorted and has duplicates and thus unsafe
    pub(super) unsafe fn subsel_from_sorted_vec_unchecked<KO: SelectionKind>(
        &self,
        index: Vec<usize>,
    ) -> Result<Sel<KO>, SelectionError> {
        if index.len() > 0 {
            Ok(Sel {
                topology: self.topology.new_ref_with_kind(),
                state: self.state.new_ref_with_kind(),
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
                topology: self.topology.new_ref_with_kind(),
                state: self.state.new_ref_with_kind(),
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
        sel.validate_index()?;
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
    fn split_collect_internal<RT, F, C, KO>(&self, func: F) -> Result<C, SelectionError>
    where
        RT: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(Particle) -> Option<RT>,
        C: FromIterator<Sel<KO>> + Default,
        KO: SelectionKind,
    {
        self.validate_index()?;
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
            values
                .into_values()
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
    pub fn split_into<RT, F, C>(self, func: F) -> Result<C, SelectionError>
    where
        RT: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(Particle) -> Option<RT>,
        C: FromIterator<Sel<K>> + Default,
    {
        self.split_collect_internal(func)
    }

    /// Helper method that splits selection into the parts with distinct resids.
    /// Parent selection is consumed.
    pub fn split_resid_into<C>(self) -> Result<C, SelectionError>
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

    pub fn split_mol_into<C>(self) -> Result<C, SelectionError>
    where
        C: FromIterator<Sel<K>> + Default,
    {
        // Error if no molecules
        if self.topology.num_molecules() == 0 {
            return Err(SelectionError::NoMolecules);
        }

        // Iterate over molecules and find those inside selection
        let first = self.first_index();
        let last = self.last_index();

        Ok(self
            .topology
            .iter_molecules()
            .cloned()
            .filter_map(|[b, e]| {
                if b < first && e >= first && e <= last {
                    // molecule starts before Sel
                    Some(0..=e - first)
                } else if b >= first && e <= last {
                    // molecule inside Sel
                    Some(b - first..=e - first)
                } else if b >= first && b <= last && e > last {
                    // molecule ends after Sel
                    Some(b - first..=last - first)
                } else {
                    None
                }
            })
            .map(|r| self.subsel(r))
            .collect::<Result<C, SelectionError>>()?)
    }

    //--------------------------
    // Non-consuming splitters
    //-------------------------

    /// Splits selection to pieces that could be disjoint
    /// according to the value of function. Parent selection is kept intact.
    ///
    /// The number of selections correspond to the distinct values returned by `func`.
    /// Selections are stored in a container `C` and has the same kind as subselections.
    pub fn split<RT, F, C>(&self, func: F) -> Result<C, SelectionError>
    where
        RT: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(Particle) -> Option<RT>,
        C: FromIterator<Sel<K>> + Default,
    {
        self.split_collect_internal(func)
    }

    /// Helper method that splits selection into the parts with distinct resids.
    /// Parent selection is left alive.
    pub fn split_resid<C>(&self) -> Result<C, SelectionError>
    where
        C: FromIterator<Sel<K>> + Default,
    {
        self.split_collect_internal(|p| Some(p.atom.resid))
    }

    /// Helper method that splits selection into the parts with distinct resindexes.
    /// Parent selection is left alive.
    pub fn split_resindex<C>(&self) -> Result<C, SelectionError>
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
    pub fn split_iter<'a, RT, F>(
        &'a self,
        func: F,
    ) -> Result<FragmentsIterator<'a, RT, F, K>, SelectionError>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(Particle) -> Option<RT> + 'a,
    {
        self.validate_index()?;
        Ok(FragmentsIterator::new(self, func))
    }

    /// Return serial iterator over contigous pieces of selection with distinct contigous resids.
    /// Parent selection is left alive.
    pub fn split_resindex_iter<'a>(
        &'a self,
    ) -> Result<
        FragmentsIterator<'a, usize, impl Fn(Particle) -> Option<usize> + 'a, K>,
        SelectionError,
    > {
        self.split_iter(|p| Some(p.atom.resindex))
    }

    pub fn split_molecules_iter(&self) -> Result<MoleculesIterator<K>, SelectionError> {
        Ok(MoleculesIterator { sel: self, cur: 0 })
    }

    //============================
    // Misceleneous
    //============================

    /// Tests if two selections are from the same source
    pub fn is_same_source(&self, other: &Sel<K>) -> bool {
        self.topology.same_data(&other.topology) && self.state.same_data(&other.state)
    }

    /// Computes the Solvet Accessible Surface Area (SASA).
    pub fn sasa(&self) -> SasaResults {
        molar_powersasa::compute_sasa(
            self.len(),
            0.14,
            |i| unsafe {
                let ind = *self.index().get_unchecked(i);
                self.state.nth_pos_mut_unchecked(ind).coords.as_mut_ptr()
            },
            |i: usize| self.nth_particle(i).unwrap().atom.vdw(),
        )
    }

    //===================
    // Subselections
    //===================

    /// Subselection from expression
    pub fn subsel(&self, def: impl SelectionDef) -> Result<Sel<K>, SelectionError> {
        Self::from_holders_and_index(
            unsafe {self.topology.new_ref_with_kind()},
            unsafe {self.state.new_ref_with_kind()},
            def.into_sel_index(&self.topology, &self.state, Some(self.index().as_slice()))?,
        )
    }

    //==============================================
    // Direct creation of selections without Source
    //==============================================
    pub fn new(
        topology: impl Into<Holder<Topology, K>>,
        state: impl Into<Holder<State, K>>,
        def: impl SelectionDef,
    ) -> Result<Self, SelectionError> {
        let topology = topology.into();
        let state = state.into();
        check_topology_state_sizes(&topology, &state)?;
        Self::from_holders_and_index(
            unsafe {topology.new_ref_with_kind()},
            unsafe {state.new_ref_with_kind()},
            def.into_sel_index(&topology, &state, None)?,
        )
    }

    /// Selects all
    pub fn new_all(
        topology: impl Into<Holder<Topology, K>>,
        state: impl Into<Holder<State, K>>,
    ) -> Result<Self, SelectionError> {
        let topology = topology.into();
        let state = state.into();
        check_topology_state_sizes(&topology, &state)?;
        Self::from_holders_and_index(
            unsafe {topology.new_ref_with_kind()},
            unsafe {state.new_ref_with_kind()},
            unsafe { SortedSet::from_sorted((0..topology.num_atoms()).collect()) },
        )
    }

    //==============================================
    // Logic on selections that modify existing ones
    //==============================================
    pub fn add(&mut self, other: impl SelectionDef) -> Result<(), SelectionError> {
        self.index_storage.extend(
            other
                .into_sel_index(&self.topology, &self.state, None)?
                .iter()
                .cloned(),
        );
        self.validate_index()?;
        Ok(())
    }

    pub fn remove_global(&mut self, other: impl SelectionDef) -> Result<(), SelectionError> {
        let lhs = rustc_hash::FxHashSet::<usize>::from_iter(self.iter_index());
        let rhs = rustc_hash::FxHashSet::<usize>::from_iter(
            other
                .into_sel_index(&self.topology, &self.state, None)?
                .iter()
                .cloned(),
        );
        let v: Vec<usize> = lhs.difference(&rhs).cloned().collect();
        self.index_storage = SortedSet::from_unsorted(v);
        Ok(())
    }

    pub fn remove_local(&mut self, other: impl IntoIterator<Item = usize>) {
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
            topology: unsafe {self.topology.new_ref_with_kind()},
            state: unsafe {self.state.new_ref_with_kind()},
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
            topology: unsafe {self.topology.new_ref_with_kind()},
            state: unsafe {self.state.new_ref_with_kind()},
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
            topology: unsafe {self.topology.new_ref_with_kind()},
            state: unsafe {self.state.new_ref_with_kind()},
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
            topology: unsafe {self.topology.new_ref_with_kind()},
            state: unsafe {self.state.new_ref_with_kind()},
            index_storage: ind,
        })
    }

    pub fn to_gromacs_ndx_str(&self, name: impl AsRef<str>) -> String {
        let name = name.as_ref();
        let mut s = format!("[ {} ]\n", name);
        for chunk in &self.iter_index().chunks(15) {
            let line: String = chunk.map(|i| (i + 1).to_string()).join(" ");
            s.push_str(&line);
            s.push('\n');
        }
        s
    }

    pub fn from_gromacs_ndx_str(
        topology: &Holder<Topology, K>,
        state: &Holder<State, K>,
        ndx_str: impl AsRef<str>,
        group_name: impl AsRef<str>,
    ) -> Result<Self, SelectionError> {
        let group_name = group_name.as_ref();
        let ndx_str = ndx_str.as_ref();

        // Find the group header
        let group_header = format!("[ {} ]", group_name);
        let mut found_group = false;
        let mut numbers = Vec::new();

        for line in ndx_str.lines() {
            let line = line.trim();

            if line == group_header {
                found_group = true;
                continue;
            }

            // If we hit another group header, stop reading
            if found_group && line.starts_with('[') {
                break;
            }

            // Parse numbers if we're in the right group
            if found_group && !line.is_empty() {
                numbers.extend(
                    line.split_whitespace()
                        .map(|s| s.parse::<usize>())
                        .map_ok(|i| i - 1) // Convert to zero-based
                        .collect::<Result<Vec<_>, _>>()
                        .map_err(|e| NdxError::Parse(group_name.into(), e))?,
                );
            }
        }

        if !found_group {
            return Err(NdxError::NoGroup(group_name.into()))?;
        }

        if numbers.is_empty() {
            return Err(NdxError::EmptyGroup(group_name.into()))?;
        }

        // Create and return new selection
        Self::from_holders_and_index(
            unsafe {topology.new_ref_with_kind()},
            unsafe {state.new_ref_with_kind()},
            numbers.into(),
        )
    }
}

impl<K: SerialKind> Sel<K> {
    //---------------------------------------------------------
    // Splitting into parallel selections (non-consuming)
    //---------------------------------------------------------

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
}

impl Sel<ImmutableParallel> {
    //---------------------------------------------------------
    // Splitting into parallel selections (non-consuming)
    //---------------------------------------------------------

    /// Splits selection to pieces that could be processed in parallel.
    /// Pieces may not be contigous and are arranged in random order.
    /// Parent selection is kept intact.
    ///
    /// Each produced selection correspond to the distinct values returned by `split_fn`.
    /// Selections are stored in a special container [ParallelSplit].
    pub unsafe fn split_par_disjoint<F, RT>(&self, split_fn: F) -> Result<ParallelSplit, SelectionError>
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
    pub unsafe fn split_par_contig<F, RT>(&self, split_fn: F) -> Result<ParallelSplit, SelectionError>
    where
        F: Fn(Particle) -> Option<RT>,
        RT: Default + std::hash::Hash + std::cmp::Eq,
    {
        Ok(ParallelSplit {
            parts: self.split_par_contig_internal(split_fn)?,
            _marker: Default::default(),
        })
    }
}

//══════════════════════════════════════════════
//███  Iterator over Particles
//══════════════════════════════════════════════
/// Iterator over [particles](Particle) in selection
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

impl<K: SelectionKind> ParticleIterProvider for Sel<K> {
    fn iter_particle(&self) -> impl ExactSizeIterator<Item = Particle<'_>> {
        SelectionIterator { sel: self, cur: 0 }
    }
}

//══════════════════════════════════════════════
//███  Mutable iterator over Particles
//══════════════════════════════════════════════
/// Iterator over [mutable particles](ParticleMut) in selection
pub struct SelectionIteratorMut<'a, K> {
    sel: &'a Sel<K>,
    cur: usize,
}

impl<'a, K: MutableKind> Iterator for SelectionIteratorMut<'a, K> {
    type Item = ParticleMut<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.cur < self.sel.len() {
            let ret = unsafe { self.sel.nth_particle_mut_unchecked(self.cur) };
            self.cur += 1;
            Some(ret)
        } else {
            None
        }
    }
}

impl<K: MutableKind> ExactSizeIterator for SelectionIteratorMut<'_, K> {
    fn len(&self) -> usize {
        self.sel.len()
    }
}

impl<K: MutableKind> ParticleIterMutProvider for Sel<K> {
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

impl<K: SelectionKind> TopologyIoProvider for Sel<K> {}

impl<K: SelectionKind> StateIoProvider for Sel<K> {}

impl<K: SelectionKind> TimeProvider for Sel<K> {
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

impl<K: SelectionKind> PosIterProvider for Sel<K> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        unsafe {
            self.index()
                .iter()
                .map(|i| self.state.nth_pos_unchecked(*i))
        }
    }
}

impl<K: SelectionKind> MeasurePos for Sel<K> {}

impl<K: SelectionKind> AtomIterProvider for Sel<K> {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        unsafe {
            self.index()
                .iter()
                .map(|i| self.topology.nth_atom_unchecked(*i))
        }
    }
}

impl<K: SelectionKind> MassIterProvider for Sel<K> {
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
    fn num_pos(&self) -> usize {
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

    unsafe fn nth_atom_unchecked(&self, i: usize) -> &Atom {
        let ind = *self.index().get_unchecked(i);
        self.topology.nth_atom_unchecked(ind)
    }
}

impl<K: SelectionKind> RandomAtomMutProvider for Sel<K> {
    unsafe fn nth_atom_mut_unchecked(&self, i: usize) -> &mut Atom {
        let ind = *self.index().get_unchecked(i);
        self.topology.nth_atom_mut_unchecked(ind)
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

    unsafe fn nth_molecule_unchecked(&self, i: usize) -> &[usize; 2] {
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

    unsafe fn nth_bond_unchecked(&self, i: usize) -> &[usize; 2] {
        self.topology.nth_bond_unchecked(i)
    }
}

impl<K: SelectionKind> RandomParticleProvider for Sel<K> {
    unsafe fn nth_particle_unchecked(&self, i: usize) -> Particle<'_> {
        let ind = *self.index_storage.get_unchecked(i);
        Particle {
            id: ind,
            atom: self.topology.nth_atom_unchecked(ind),
            pos: self.state.nth_pos_unchecked(ind),
        }
    }
}

//═══════════════════════════════════════════════════════════
//███  Mutable analysis traits (only for mutable selections)
//═══════════════════════════════════════════════════════════

impl<K: MutableKind> PosIterMutProvider for Sel<K> {
    fn iter_pos_mut(&self) -> impl PosMutIterator<'_> {
        unsafe {
            self.index()
                .iter()
                .map(|i| self.state.nth_pos_mut_unchecked(*i))
        }
    }
}

impl<K: MutableKind> AtomsIterMutProvider for Sel<K> {
    fn iter_atoms_mut(&self) -> impl AtomMutIterator<'_> {
        unsafe {
            self.index()
                .iter()
                .map(|i| self.topology.nth_atom_mut_unchecked(*i))
        }
    }
}

impl<K: MutableKind> RandomPosMutProvider for Sel<K> {
    fn nth_pos_mut(&self, i: usize) -> Option<&mut Pos> {
        self.index()
            .get(i)
            .map(|i| unsafe { self.state.nth_pos_mut_unchecked(*i) })
    }

    unsafe fn nth_pos_mut_unchecked(&self, i: usize) -> &mut Pos {
        let ind = *self.index().get_unchecked(i);
        self.state.nth_pos_mut_unchecked(ind)
    }
}

impl<K: MutableKind> RandomParticleMutProvider for Sel<K> {
    unsafe fn nth_particle_mut_unchecked(&self, i: usize) -> ParticleMut {
        let ind = *self.index_storage.get_unchecked(i);
        ParticleMut {
            id: ind,
            atom: self.topology.nth_atom_mut_unchecked(ind),
            pos: self.state.nth_pos_mut_unchecked(ind),
        }
    }
}


impl<K: MutableKind> ModifyPos for Sel<K> {}
impl<K: MutableKind> ModifyPeriodic for Sel<K> {}
impl<K: MutableKind> ModifyRandomAccess for Sel<K> {}
