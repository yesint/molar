use super::{
    measure::{MeasureMasses, MeasurePeriodic, MeasurePos},
    modify::{ModifyPos, ModifyRandomAccess},
    providers::*,
    Atom, AtomIterator, PeriodicBox, Pos, PosIterator, PosMutIterator, State, StateArc, StateUArc,
    Topology, TopologyArc, TopologyUArc,
};
use crate::io::{FileHandler, IndexProvider, StateProvider, TopologyProvider};
use anyhow::{anyhow, bail, Result};
use itertools::Itertools;
use std::{collections::HashMap, marker::PhantomData, ops::Range};

pub use super::selection_parser::SelectionExpr;
//-----------------------------------------------------------------
pub trait SelectionKind {
    type SubselType;
    const NEED_CHECK_OVERLAP: bool;
}
pub trait MayOverlap: SelectionKind {}
pub trait MutableSel: SelectionKind {}

// Read-write unique marker (multithreaded)
pub struct NonOverlappingMut {}
impl SelectionKind for NonOverlappingMut {
    // Subseletions may overlap but won't be Send
    type SubselType = OverlappingMut;
    const NEED_CHECK_OVERLAP: bool = true;
}
impl MutableSel for NonOverlappingMut {}

// Read-write marker (single-threaded)
pub struct OverlappingMut(*const ());
impl SelectionKind for OverlappingMut {
    type SubselType = OverlappingMut;
    const NEED_CHECK_OVERLAP: bool = false;
}
impl MayOverlap for OverlappingMut {}
impl MutableSel for OverlappingMut {}

// Read-only marker (multithreaded )
pub struct Overlapping {}
impl SelectionKind for Overlapping {
    type SubselType = Overlapping;
    const NEED_CHECK_OVERLAP: bool = false;
}
impl MayOverlap for Overlapping {}

//=====================================================

pub struct Source<T> {
    topology: TopologyArc,
    state: StateArc,
    used: rustc_hash::FxHashSet<usize>,
    _marker: PhantomData<T>,
}

impl Source<()> {
    pub fn new_non_overlapping_mut(
        topology: TopologyUArc,
        state: StateUArc,
    ) -> Result<Source<NonOverlappingMut>> {
        check_sizes(&topology, &state)?;
        Ok(Source {
            topology: topology.shareable(),
            state: state.shareable(),
            used: Default::default(),
            _marker: Default::default(),
        })
    }

    pub fn new_overlapping_mut(
        topology: TopologyUArc,
        state: StateUArc,
    ) -> Result<Source<OverlappingMut>> {
        check_sizes(&topology, &state)?;
        Ok(Source {
            topology: topology.shareable(),
            state: state.shareable(),
            used: Default::default(),
            _marker: Default::default(),
        })
    }

    pub fn new_overlapping(
        topology: TopologyUArc,
        state: StateUArc,
    ) -> Result<Source<Overlapping>> {
        check_sizes(&topology, &state)?;
        Ok(Source {
            topology: topology.shareable(),
            state: state.shareable(),
            used: Default::default(),
            _marker: Default::default(),
        })
    }
}

impl<T> Source<T> {
    pub fn release(self) -> anyhow::Result<(TopologyUArc, StateUArc)> {
        Ok((
            triomphe::Arc::try_unique(self.topology)
                .or_else(|_| bail!("Multiple references are active!"))?,
            triomphe::Arc::try_unique(self.state)
                .or_else(|_| bail!("Multiple references are active!"))?,
        ))
    }
}

// Conversions to other builder types
impl Source<NonOverlappingMut> {
    pub fn to_overlapping_mut(self) -> anyhow::Result<Source<OverlappingMut>> {
        let (top, st) = self.release()?;
        Ok(Source::new_overlapping_mut(top, st)?)
    }

    pub fn to_overlapping(self) -> anyhow::Result<Source<Overlapping>> {
        let (top, st) = self.release()?;
        Ok(Source::new_overlapping(top, st)?)
    }
}

impl Source<OverlappingMut> {
    pub fn to_non_overlapping_mut(self) -> anyhow::Result<Source<NonOverlappingMut>> {
        let (top, st) = self.release()?;
        Ok(Source::new_non_overlapping_mut(top, st)?)
    }

    pub fn to_overlapping(self) -> anyhow::Result<Source<Overlapping>> {
        let (top, st) = self.release()?;
        Ok(Source::new_overlapping(top, st)?)
    }
}

impl Source<Overlapping> {
    pub fn to_non_overlapping_mut(self) -> anyhow::Result<Source<NonOverlappingMut>> {
        let (top, st) = self.release()?;
        Ok(Source::new_non_overlapping_mut(top, st)?)
    }

    pub fn to_overlapping(self) -> anyhow::Result<Source<OverlappingMut>> {
        let (top, st) = self.release()?;
        Ok(Source::new_overlapping_mut(top, st)?)
    }
}

//---------------------------------
// Creating selections
//---------------------------------

impl<T: SelectionKind> Source<T> {
    pub fn get_topology(&self) -> TopologyArc {
        triomphe::Arc::clone(&self.topology)
    }

    pub fn get_state(&self) -> StateArc {
        triomphe::Arc::clone(&self.state)
    }

    fn check_overlap_if_needed(&mut self, index: &Vec<usize>) -> anyhow::Result<()> {
        if T::NEED_CHECK_OVERLAP {
            for i in index.iter() {
                if !self.used.insert(*i) {
                    bail!("Index {i} is already used!");
                }
            }
        }
        Ok(())
    }

    pub fn set_state(&mut self, state: StateUArc) -> Result<StateUArc> {
        // Check if the states are compatible
        if !self.state.interchangeable(&state) {
            bail!("States are incompatible!")
        }

        let ret = state.shareable();
        unsafe {
            std::ptr::swap(
                self.state.as_ptr() as *mut State,
                ret.as_ptr() as *mut State,
            )
        }
        triomphe::Arc::try_unique(ret).or_else(|_| bail!("Multiple references are active!"))
    }

    pub fn set_topology(&mut self, topology: TopologyUArc) -> Result<TopologyUArc> {
        // Check if the states are compatible
        if !self.topology.interchangeable(&topology) {
            bail!("Topologies are incompatible!")
        }

        let ret = topology.shareable();
        unsafe {
            std::ptr::swap(
                self.topology.as_ptr() as *mut State,
                ret.as_ptr() as *mut State,
            )
        }
        triomphe::Arc::try_unique(ret).or_else(|_| bail!("Multiple references are active!"))
    }

    pub fn select_from_iter(
        &mut self,
        iter: impl Iterator<Item = usize>,
    ) -> anyhow::Result<Sel<T>> {
        let vec = index_from_iter(iter, self.topology.num_atoms())?;
        self.check_overlap_if_needed(&vec)?;
        Ok(Sel {
            topology: triomphe::Arc::clone(&self.topology),
            state: triomphe::Arc::clone(&self.state),
            index: vec,
            _marker: PhantomData::default(),
        })
    }

    pub fn select_all(&mut self) -> anyhow::Result<Sel<T>> {
        let vec = index_from_all(self.topology.num_atoms());
        self.check_overlap_if_needed(&vec)?;
        Ok(Sel {
            topology: triomphe::Arc::clone(&self.topology),
            state: triomphe::Arc::clone(&self.state),
            index: vec,
            _marker: PhantomData::default(),
        })
    }

    pub fn select_str(&mut self, selstr: &str) -> anyhow::Result<Sel<T>> {
        let vec = index_from_str(selstr, &self.topology, &self.state)?;
        self.check_overlap_if_needed(&vec)?;
        Ok(Sel {
            topology: triomphe::Arc::clone(&self.topology),
            state: triomphe::Arc::clone(&self.state),
            index: vec,
            _marker: PhantomData::default(),
        })
    }

    pub fn select_expr(&mut self, expr: &SelectionExpr) -> anyhow::Result<Sel<T>> {
        let vec = index_from_expr(expr, &self.topology, &self.state)?;
        self.check_overlap_if_needed(&vec)?;
        Ok(Sel {
            topology: triomphe::Arc::clone(&self.topology),
            state: triomphe::Arc::clone(&self.state),
            index: vec,
            _marker: PhantomData::default(),
        })
    }

    pub fn select_range(&mut self, range: &std::ops::Range<usize>) -> anyhow::Result<Sel<T>> {
        let vec = index_from_range(range, self.topology.num_atoms())?;
        self.check_overlap_if_needed(&vec)?;
        Ok(Sel {
            topology: triomphe::Arc::clone(&self.topology),
            state: triomphe::Arc::clone(&self.state),
            index: vec,
            _marker: PhantomData::default(),
        })
    }
}

// IO traits
impl<T> TopologyProvider for Source<T> {
    fn num_atoms(&self) -> usize {
        self.topology.num_atoms()
    }
}

impl<T> AtomsProvider for Source<T> {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        self.topology.iter_atoms()
    }
}

impl<T> StateProvider for Source<T> {
    fn get_time(&self) -> f32 {
        self.state.get_time()
    }

    fn num_coords(&self) -> usize {
        self.state.num_coords()
    }
}

impl<T> BoxProvider for Source<T> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.state.get_box()
    }
}

impl<T> PosProvider for Source<T> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.state.iter_pos()
    }
}

//-----------------------------------------------------------------
fn check_sizes(topology: &Topology, state: &State) -> Result<()> {
    let n1 = topology.num_atoms();
    let n2 = state.num_coords();
    match n1 == n2 {
        true => Ok(()),
        false => bail!(
            "Structure and State are incompatible (sizes {} and {})",
            n1,
            n2
        ),
    }
}

fn index_from_all(n: usize) -> Vec<usize> {
    (0..n).collect()
}

fn index_from_expr(expr: &SelectionExpr, topology: &Topology, state: &State) -> Result<Vec<usize>> {
    let index = expr.apply_whole(&topology, &state)?;
    if index.len() > 0 {
        Ok(index)
    } else {
        bail!("Selection is empty")
    }
}

fn index_from_str(selstr: &str, topology: &Topology, state: &State) -> Result<Vec<usize>> {
    let index = SelectionExpr::try_from(selstr)?.apply_whole(&topology, &state)?;
    if index.len() > 0 {
        Ok(index)
    } else {
        bail!("Selection is empty")
    }
}

fn index_from_range(range: &Range<usize>, n: usize) -> Result<Vec<usize>> {
    if range.start > n - 1 || range.end > n - 1 {
        bail!(
            "Index range {}:{} is invalid, 0:{} is allowed.",
            range.start,
            range.end,
            n
        );
    }
    if range.len() > 0 {
        Ok(range.clone().collect())
    } else {
        bail!("Selection is empty")
    }
}

fn index_from_iter(it: impl Iterator<Item = usize>, n: usize) -> Result<Vec<usize>> {
    let index: Vec<usize> = it.sorted().dedup().collect();
    if index.is_empty() {
        bail!("Selection is empty")
    }
    if index[0] > n - 1 || index[index.len() - 1] > n - 1 {
        bail!(
            "Index range {}:{} is invalid, 0:{} is allowed.",
            index[0],
            index[index.len() - 1],
            n
        );
    }
    Ok(index)
}

//---------------------------------------
// Selection
//---------------------------------------

pub struct Sel<T> {
    topology: TopologyArc,
    state: StateArc,
    index: Vec<usize>,
    _marker: PhantomData<T>,
}

impl<T: SelectionKind> Sel<T> {
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.num_atoms()
    }

    //===================
    // Subselections
    //===================

    /// Subselection from expression
    pub fn subsel_from_expr(&self, expr: &SelectionExpr) -> Result<Sel<T::SubselType>> {
        let index = expr.apply_subset(&self.topology, &self.state, self.index.iter().cloned())?;
        if index.len() > 0 {
            Ok(Sel {
                topology: triomphe::Arc::clone(&self.topology),
                state: triomphe::Arc::clone(&self.state),
                index,
                _marker: Default::default(),
            })
        } else {
            bail!("Selection is empty")
        }
    }

    /// Subselection from string
    pub fn subsel_from_str(&self, sel_str: &str) -> Result<Sel<T::SubselType>> {
        let expr = SelectionExpr::try_from(sel_str)?;
        self.subsel_from_expr(&expr)
    }

    /// Subselection from the range of local selection indexes
    pub fn subsel_from_local_range(
        &self,
        range: std::ops::Range<usize>,
    ) -> Result<Sel<T::SubselType>> {
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

        Ok(Sel {
            topology: triomphe::Arc::clone(&self.topology),
            state: triomphe::Arc::clone(&self.state),
            index,
            _marker: Default::default(),
        })
    }

    /// Subselection from iterator that provides local selection indexes
    pub fn subsel_from_iter(
        &self,
        iter: impl ExactSizeIterator<Item = usize>,
    ) -> Result<Sel<T::SubselType>> {
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
    unsafe fn subsel_from_vec_unchecked<S>(&self, index: Vec<usize>) -> Result<Sel<S>> {
        if index.len() > 0 {
            Ok(Sel {
                topology: triomphe::Arc::clone(&self.topology),
                state: triomphe::Arc::clone(&self.state),
                index,
                _marker: Default::default(),
            })
        } else {
            bail!("Selection is empty")
        }
    }

    pub unsafe fn nth_unchecked(&self, i: usize) -> (usize, &Atom, &Pos) {
        let ind = *self.index.get_unchecked(i);
        (
            ind,
            self.topology.nth_atom_unchecked(ind),
            self.state.nth_pos_unchecked(ind),
        )
    }

    pub unsafe fn nth_unchecked_mut(&self, i: usize) -> (usize, &mut Atom, &mut Pos) {
        let ind = *self.index.get_unchecked(i);
        (
            ind,
            self.topology.nth_atom_unchecked_mut(ind),
            self.state.nth_pos_unchecked_mut(ind),
        )
    }

    // Helper splitting function generic over returned selections kind
    fn split_gen<RT, F, C, Kind>(&self, func: F) -> C
    where
        RT: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(usize, &Atom, &Pos) -> RT,
        C: FromIterator<Sel<Kind>> + Default,
    {
        let mut ids = HashMap::<RT, Vec<usize>>::default();

        for i in 0..self.index.len() {
            let (ind, at, pos) = unsafe { self.nth_unchecked(i) };
            let id = func(ind, at, pos);
            if let Some(el) = ids.get_mut(&id) {
                el.push(ind);
            } else {
                ids.insert(id, vec![ind]);
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
    /// The number of selections correspond to the distinct values returned by `func`.
    /// Selections are stored in a container C and has the same kind as subselections.
    pub fn split<RT, F, C>(&self, func: F) -> C
    where
        RT: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(usize, &Atom, &Pos) -> RT,
        C: FromIterator<Sel<T::SubselType>> + Default,
    {
        self.split_gen(func)
    }

    /// Splits selection to pieces that could be disjoint
    /// according to the value of function. Parent selection is consumed.
    /// The number of selections correspond to the distinct values returned by `func`.
    /// Selections are stored in a container C and has the same kind as parent selection.
    pub fn split_into<RT, F, C>(self, func: F) -> C
    where
        RT: Default + std::hash::Hash + std::cmp::Eq,
        F: Fn(usize, &Atom, &Pos) -> RT,
        C: FromIterator<Sel<T>> + Default,
    {
        self.split_gen(func)
    }

    pub fn split_resid<C>(&self) -> C
    where
        C: FromIterator<Sel<T::SubselType>> + Default,
    {
        self.split_gen(|_, at, _| at.resid)
    }

    pub fn split_resid_into<C>(self) -> C
    where
        C: FromIterator<Sel<T>> + Default,
    {
        self.split_gen(|_, at, _| at.resid)
    }

    // Sasa
    pub fn sasa(&self) -> (f32, f32) {
        let (areas, volumes) = molar_powersasa::sasa(
            self.len(),
            0.14,
            |i| unsafe {
                let ind = *self.index.get_unchecked(i);
                self.state.nth_pos_unchecked_mut(ind).coords.as_mut_ptr()
            },
            |i: usize| self.nth(i).unwrap().1.vdw(),
        );
        (areas.into_iter().sum(), volumes.into_iter().sum())
    }

    pub fn get_index_vec(&self) -> Vec<usize> {
        self.index.clone()
    }

    pub fn save(&self, fname: &str) -> Result<()> {
        let mut h = FileHandler::create(fname)?;
        h.write(self)
    }

    pub fn first(&self) -> (usize, &Atom, &Pos) {
        unsafe { self.nth_unchecked(0) }
    }

    pub fn nth(&self, i: usize) -> Result<(usize, &Atom, &Pos)> {
        if i > self.len() {
            bail!("Index {} is beyond the allowed range [0:{}]", i, self.len())
        }
        Ok(unsafe { self.nth_unchecked(i) })
    }

    pub fn iter(&self) -> SelectionIterator<T> {
        SelectionIterator { sel: self, cur: 0 }
    }

    /// Return iterator that splits selection into contigous pieces according to the value of function.
    /// Whenever `func` returns a value different from the previous one, new selection is created.
    /// Selections are computed lazily when iterating.
    /// This consumer a selection and returns
    /// selections of the same kind.
    pub fn into_split_contig<RT, F>(self, func: F) -> IntoSelectionSplitIterator<RT, F, T>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(usize, &Atom, &Pos) -> RT,
    {
        IntoSelectionSplitIterator::new(self, func)
    }

    // Pre-defined splitters
    pub fn into_split_contig_resid(
        self,
    ) -> IntoSelectionSplitIterator<i32, fn(usize, &Atom, &Pos) -> i32, T> {
        self.into_split_contig(|_, at, _| at.resid)
    }

    /// Return iterator that splits selection into contigous pieces according to the value of function.
    /// Whenever `func` returns a value different from the previous one, new selection is created.
    /// Selections are computed lazily when iterating.
    /// Selection is not consumed. Returned selections are of subselection kind.
    pub fn split_contig<RT, F>(&self, func: F) -> SelectionSplitIterator<'_, RT, F, T>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(usize, &Atom, &Pos) -> RT,
    {
        SelectionSplitIterator::new(self, func)
    }

    // Pre-defined splitters
    pub fn split_contig_resid(
        &self,
    ) -> SelectionSplitIterator<'_, i32, fn(usize, &Atom, &Pos) -> i32, T> {
        self.split_contig(|_, at, _| at.resid)
    }
}

//---------------------------------------------
pub struct SelectionIterator<'a, T> {
    sel: &'a Sel<T>,
    cur: usize,
}

impl<'a, T: SelectionKind> Iterator for SelectionIterator<'a, T> {
    type Item = (usize, &'a Atom, &'a Pos);
    fn next(&mut self) -> Option<Self::Item> {
        if self.cur < self.sel.len() {
            let ret = unsafe { self.sel.nth_unchecked(self.cur) };
            self.cur += 1;
            Some(ret)
        } else {
            None
        }
    }
}

//---------------------------------------------
// Implement traits for IO

impl<T> IndexProvider for Sel<T> {
    fn iter_index(&self) -> impl Iterator<Item = usize> {
        self.index.iter().cloned()
    }
}

impl<T> TopologyProvider for Sel<T> {
    fn num_atoms(&self) -> usize {
        self.index.len()
    }
}

impl<T> StateProvider for Sel<T> {
    fn get_time(&self) -> f32 {
        self.state.get_time()
    }

    fn num_coords(&self) -> usize {
        self.index.len()
    }
}

//==================================================================
// Implement analysis traits

impl<T> BoxProvider for Sel<T> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.state.get_box()
    }
}

impl<T> MeasurePeriodic for Sel<T> {}

impl<T> PosProvider for Sel<T> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        unsafe { self.index.iter().map(|i| self.state.nth_pos_unchecked(*i)) }
    }
}

impl<T> MeasurePos for Sel<T> {}

impl<T> AtomsProvider for Sel<T> {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        unsafe {
            self.index
                .iter()
                .map(|i| self.topology.nth_atom_unchecked(*i))
        }
    }
}

impl<T> MassesProvider for Sel<T> {
    fn iter_masses(&self) -> impl ExactSizeIterator<Item = f32> {
        unsafe {
            self.index
                .iter()
                .map(|i| self.topology.nth_atom_unchecked(*i).mass)
        }
    }
}

impl<T> MeasureMasses for Sel<T> {}

//-------------------------------------------------------
// Mutable analysis traits only for mutable selections

impl<T: MutableSel> PosMutProvider for Sel<T> {
    fn iter_pos_mut(&self) -> impl PosMutIterator<'_> {
        unsafe {
            self.index
                .iter()
                .map(|i| self.state.nth_pos_unchecked_mut(*i))
        }
    }
}

impl<T: MutableSel> RandomPosMutProvider for Sel<T> {
    #[inline(always)]
    unsafe fn nth_pos_unchecked_mut(&self, i: usize) -> &mut Pos {
        let ind = *self.index.get_unchecked(i);
        self.state.nth_pos_unchecked_mut(ind)
    }
}

impl<T: MutableSel> ModifyPos for Sel<T> {}
impl<T: MutableSel> ModifyRandomAccess for Sel<T> {}

//-------------------------------------------------------
// Splitting iterator
//-------------------------------------------------------
pub struct SplitData<RT, F> {
    func: F,
    counter: usize,
    id: RT,
}

pub struct SelectionSplitIterator<'a, RT, F, S> {
    sel: &'a Sel<S>,
    data: SplitData<RT, F>,
}

impl<S> SelectionSplitIterator<'_, (), (), S> {
    pub fn new<RT, F>(sel: &Sel<S>, func: F) -> SelectionSplitIterator<'_, RT, F, S>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(usize, &Atom, &Pos) -> RT,
    {
        SelectionSplitIterator {
            sel,
            data: SplitData {
                func,
                counter: 0,
                id: RT::default(),
            },
        }
    }
}

impl<RT, F, S> Iterator for SelectionSplitIterator<'_, RT, F, S>
where
    RT: Default + std::cmp::PartialEq,
    F: Fn(usize, &Atom, &Pos) -> RT,
    S: SelectionKind,
{
    // Non-consuming splitter returns subselections
    type Item = Sel<S::SubselType>;

    fn next(&mut self) -> Option<Self::Item> {
        next_split(&mut self.data, self.sel)
    }
}

//--------------------------------------------
pub struct IntoSelectionSplitIterator<RT, F, S> {
    sel: Sel<S>,
    data: SplitData<RT, F>,
}

impl<S> IntoSelectionSplitIterator<(), (), S> {
    pub fn new<RT, F>(sel: Sel<S>, func: F) -> IntoSelectionSplitIterator<RT, F, S>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(usize, &Atom, &Pos) -> RT,
    {
        IntoSelectionSplitIterator {
            sel,
            data: SplitData {
                func,
                counter: 0,
                id: RT::default(),
            },
        }
    }
}

impl<RT, F, S> Iterator for IntoSelectionSplitIterator<RT, F, S>
where
    RT: Default + std::cmp::PartialEq,
    F: Fn(usize, &Atom, &Pos) -> RT,
    S: SelectionKind,
{
    // Consuming splitter always return the same selection kind as parent
    type Item = Sel<S>;

    fn next(&mut self) -> Option<Self::Item> {
        next_split(&mut self.data, &self.sel)
    }
}

// Actual function that does the splitting
fn next_split<RT, F, S, SR>(data: &mut SplitData<RT, F>, sel: &Sel<S>) -> Option<Sel<SR>>
where
    RT: Default + std::cmp::PartialEq,
    F: Fn(usize, &Atom, &Pos) -> RT,
    S: SelectionKind,
{
    let mut index = Vec::<usize>::new();
    while data.counter < sel.len() {
        let (i, at, pos) = unsafe { sel.nth_unchecked(data.counter) };
        let id = (data.func)(i, at, pos);

        if id == data.id {
            // Current selection continues. Add current index
            index.push(i);
        } else if index.is_empty() {
            // The very first id is not default, this is Ok, add index
            // and update self.id
            data.id = id;
            index.push(i);
        } else {
            // The end of current selection
            data.id = id; // Update self.id for the next selection
            return unsafe { Some(sel.subsel_from_vec_unchecked(index).unwrap()) };
        }
        // Next element
        data.counter += 1;
    }

    // Return any remaining index as last selection
    if !index.is_empty() {
        return unsafe { Some(sel.subsel_from_vec_unchecked(index).unwrap()) };
    }

    // If we are here stop iterating
    None
}
//############################################################
//#  Tests
//############################################################

#[cfg(test)]
mod tests {

    use rayon::iter::{ParallelBridge, ParallelIterator};

    use super::{OverlappingMut, Sel, Source};
    use crate::{
        core::{
            providers::PosProvider, selection::AtomsProvider, MeasureMasses, MeasurePeriodic,
            MeasurePos, ModifyPos, ModifyRandomAccess, Pos, State, Topology, Vector3f, PBC_FULL,
        },
        io::*,
    };

    fn read_test_pdb() -> (triomphe::UniqueArc<Topology>, triomphe::UniqueArc<State>) {
        let mut h = FileHandler::open("tests/no_ATP.pdb").unwrap();
        let (top, st) = h.read_raw().unwrap();
        (triomphe::UniqueArc::new(top), triomphe::UniqueArc::new(st))
    }

    #[test]
    fn builder_overlap() -> anyhow::Result<()> {
        let (top, st) = read_test_pdb();
        let mut b = Source::new_overlapping_mut(top, st)?;
        // Create two overlapping selections
        let _sel1 = b.select_from_iter(0..10)?;
        let _sel2 = b.select_from_iter(5..15)?;
        Ok(())
    }

    #[test]
    fn builder_par_no_overlap() {
        let (top, st) = read_test_pdb();
        let mut b = Source::new_non_overlapping_mut(top, st).unwrap();
        // Create two non-overlapping selections.
        let _sel1 = b.select_from_iter(0..10).unwrap();
        let _sel2 = b.select_from_iter(11..15).unwrap();
    }

    #[test]
    #[should_panic]
    fn builder_par_overlap() {
        let (top, st) = read_test_pdb();
        let mut b = Source::new_non_overlapping_mut(top, st).unwrap();
        // Create two overlapping selections. This must fail!
        let _sel1 = b.select_from_iter(0..10).unwrap();
        let _sel2 = b.select_from_iter(5..15).unwrap();
    }

    #[test]
    fn builder_par_threads() -> anyhow::Result<()> {
        let (top, st) = read_test_pdb();
        let mut b = Source::new_non_overlapping_mut(top, st)?;
        // Create two valid non-overlapping selections.
        let sel1 = b.select_from_iter(0..10).unwrap();
        let sel2 = b.select_from_iter(11..15).unwrap();
        // Pass them to threads
        let t1 = std::thread::spawn(move || -> anyhow::Result<Pos> {
            sel1.translate(Vector3f::new(1.0, 2.0, 3.0));
            Ok(sel1.center_of_mass()?)
        });

        let t2 = std::thread::spawn(move || -> anyhow::Result<Pos> {
            sel2.translate(Vector3f::new(3.0, 2.0, 1.0));
            Ok(sel2.center_of_mass()?)
        });

        let cm1 = t1.join().unwrap()?;
        let cm2 = t2.join().unwrap()?;
        println!("cm1 = {cm1}; cm2={cm2}");
        Ok(())
    }

    #[test]
    #[should_panic]
    fn convert_builders_fail() {
        let (top, st) = read_test_pdb();
        // Create parallel builder
        let mut b = Source::new_non_overlapping_mut(top, st).unwrap();
        // Create two valid non-overlapping selections.
        let _sel1 = b.select_from_iter(0..10).unwrap();
        let _sel2 = b.select_from_iter(11..15).unwrap();

        // Create serial builder. This will fail since selections are still alive
        let _b = b.to_overlapping().unwrap();
    }

    #[test]
    fn convert_builders() {
        let (top, st) = read_test_pdb();
        // Create parallel builder
        let mut b = Source::new_non_overlapping_mut(top, st).unwrap();
        // Create two valid non-overlapping selections.
        {
            let _sel1 = b.select_from_iter(0..10).unwrap();
            let _sel2 = b.select_from_iter(11..15).unwrap();
        }
        // Selections are now dropped

        // Create serial builder.
        let mut b = b.to_overlapping().unwrap();
        let _sel1 = b.select_from_iter(0..10).unwrap();
        let _sel2 = b.select_from_iter(11..15).unwrap();
    }

    fn make_sel_all() -> anyhow::Result<Sel<OverlappingMut>> {
        let (top, st) = read_test_pdb();
        let mut b = Source::new_overlapping_mut(top, st)?;
        let sel = b.select_all()?;
        Ok(sel)
    }

    fn make_sel_prot() -> anyhow::Result<Sel<OverlappingMut>> {
        let (top, st) = read_test_pdb();
        let mut b = Source::new_overlapping_mut(top, st)?;
        let sel = b.select_str("not resname TIP3 POT CLA")?;
        Ok(sel)
    }

    #[test]
    fn test_measure() -> anyhow::Result<()> {
        let sel = make_sel_all()?;
        println!("before {}", sel.iter_pos().next().unwrap());

        let (minv, maxv) = sel.min_max();
        println!("{minv}:{maxv}");

        //sel.translate(&Vector3f::new(10.0,10.0,10.0));
        println!("after {}", sel.iter_pos().next().unwrap());
        println!("{:?}", sel.min_max());
        Ok(())
    }

    #[test]
    fn test_measure_pbc() -> anyhow::Result<()> {
        let sel = make_sel_all()?;

        let cm = sel.center_of_mass()?;
        println!("{cm}");
        Ok(())
    }

    #[test]
    fn test_translate() -> anyhow::Result<()> {
        let sel = make_sel_all()?;

        println!("before {}", sel.iter_pos().next().unwrap());
        sel.translate(Vector3f::new(10.0, 10.0, 10.0));
        println!("after {}", sel.iter_pos().next().unwrap());
        Ok(())
    }

    #[test]
    fn test_write_to_file() -> anyhow::Result<()> {
        let sel = make_sel_all()?;

        let mut h = FileHandler::create("f.pdb")?;
        h.write_topology(&sel)?;
        h.write_state(&sel)?;
        Ok(())
    }

    #[test]
    fn test_unwrap_connectivity_1() -> anyhow::Result<()> {
        let sel = make_sel_prot()?;
        sel.unwrap_connectivity_dim(0.2, &PBC_FULL)?;

        let mut h = FileHandler::create("unwrapped.pdb")?;
        h.write_topology(&sel)?;
        h.write_state(&sel)?;
        Ok(())
    }

    #[test]
    fn eigen_test() -> anyhow::Result<()> {
        let sel1 = make_sel_prot()?;
        let sel2 = make_sel_prot()?;

        sel2.rotate(&Vector3f::x_axis(), 80.0_f32.to_radians());

        sel2.save("tests/sel2.pdb")?;
        sel1.save("tests/sel1_before.pdb")?;
        println!("Initial RMSD:{}", Sel::rmsd_mw(&sel1, &sel2)?);

        let m = Sel::fit_transform(&sel1, &sel2)?;
        println!("{m}");

        sel1.apply_transform(&m);

        sel1.save("tests/sel1_after.pdb")?;
        println!("Final RMSD:{}", Sel::rmsd_mw(&sel1, &sel2)?);
        Ok(())
    }

    #[test]
    fn split_test() -> anyhow::Result<()> {
        let sel1 = make_sel_prot()?;
        for res in sel1.split_contig(|_, at, _| at.resid) {
            println!("Res: {}", res.iter_atoms().next().unwrap().resid)
        }
        Ok(())
    }

    #[test]
    fn sasa_test() -> anyhow::Result<()> {
        let sel1 = make_sel_all()?;
        let (a, v) = sel1.sasa();
        println!("Sasa: {a}, Volume: {v}");
        Ok(())
    }

    #[test]
    fn tets_gyration() -> anyhow::Result<()> {
        let sel1 = make_sel_prot()?;
        let g = sel1.gyration_pbc()?;
        println!("Gyration radius: {g}");
        Ok(())
    }

    #[test]
    fn test_inertia() -> anyhow::Result<()> {
        let sel1 = make_sel_all()?;
        //sel1.rotate(
        //    &UnitVector3::new_normalize(Vector3f::new(1.0,2.0,3.0)),
        //    0.45
        //);
        let (moments, axes) = sel1.inertia_pbc()?;
        println!("Inertia moments: {moments}");
        println!("Inertia axes: {axes:?}");
        Ok(())
        // {0.7308828830718994 0.3332606256008148 -0.5956068634986877}
        // {0.6804488301277161 -0.28815552592277527 0.6737624406814575}
        // {-0.05291106179356575 0.8977214694023132 0.4373747706413269}
    }

    #[test]
    fn test_principal_transform() -> anyhow::Result<()> {
        let sel1 = make_sel_prot()?;
        let tr = sel1.principal_transform()?;
        println!("Transform: {tr}");

        let (_, axes) = sel1.inertia()?;
        println!("Axes before: {axes}");

        sel1.apply_transform(&tr);

        let (_, axes) = sel1.inertia()?;
        println!("Axes after: {axes}");

        sel1.save("oriented.pdb")?;

        Ok(())
    }

    #[test]
    fn test_swap_state() -> anyhow::Result<()> {
        let top = FileHandler::open("tests/no_ATP.pdb")?.read_topology()?;
        let mut traj = FileHandler::open("tests/no_ATP.xtc")?;
        let st1 = traj.read_state()?.unwrap();
        let st2 = traj.read_state()?.unwrap();

        let mut source = Source::new_overlapping(top, st1)?;
        let sel = source.select_all()?;
        println!("Before swap: {}", source.get_state().num_coords());
        println!("First point: {:?}", sel.iter_pos().next().unwrap());

        source.set_state(st2)?;

        println!("After swap: {}", source.get_state().num_coords());
        println!("First point: {:?}", sel.iter_pos().next().unwrap());
        Ok(())
    }

    #[test]
    fn test_swap_traj() -> anyhow::Result<()> {
        let top = FileHandler::open("tests/no_ATP.pdb")?.read_topology()?;
        let mut traj = FileHandler::open("tests/no_ATP.xtc")?.into_iter();

        let mut source = Source::new_overlapping(top, traj.next().unwrap())?;
        let sel = source.select_all()?;
        println!("First time {}: {:?}", sel.get_time(), sel.center_of_mass()?);

        for st in traj {
            source.set_state(st)?;
            println!("Time {}: {:?}", sel.get_time(), sel.center_of_mass()?);
        }

        Ok(())
    }

    #[test]
    fn test_sel_par_proc() -> anyhow::Result<()> {
        let (top, st) = FileHandler::open("tests/colored.pdb")?.read()?;

        let mut source = Source::new_non_overlapping_mut(top, st)?;

        let (sender, receiver) = std::sync::mpsc::channel();

        source.select_all()?
        .into_split_contig_resid()
        .enumerate()
        .par_bridge()
        .try_for_each_with(sender,|s,(i,sel)| -> anyhow::Result<()> {
            println!("Doing {i}");
            //sel.unwrap_connectivity(0.2)?;
            s.send(sel.center_of_mass()?)?;
            Ok(())
        })?;

        let res: Vec<_> = receiver.iter().collect();
        println!("cms: {:?}",res);
        Ok(())
    }
}
