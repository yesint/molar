use std::marker::PhantomData;
use sorted_vec::SortedSet;
use crate::{core::{State, Topology}, io::{StateProvider, TopologyProvider}};

use super::SelectionError;

/// Trait for kinds of selections
pub trait SelectionKind {
    #[inline(always)]
    #[allow(unused_variables)]
    fn check_index(index: &SortedSet<usize>, system: &Topology, state: &State) -> Result<(), SelectionError> {
        Ok(())
    }

    #[inline(always)]
    #[allow(unused_variables)]
    fn check_overlap(index: &SortedSet<usize>, used: &mut rustc_hash::FxHashSet<usize>) -> Result<(), SelectionError> {
        Ok(())
    }
}

/// Trait marking non-overlapping selections
pub trait MutableSel: SelectionKind {}
/// Trait marking parallel selections
pub trait ParallelSel: SelectionKind + Send + Sync {}
/// Trait marking serial selections
pub trait SerialSel: SelectionKind {}
/// Trait marking selections that allow subselections
pub trait AllowsSubselect: SelectionKind {}


/// Marker type for possibly overlapping mutable selection (single-threaded)
pub struct MutableSerial(PhantomData<*const ()>);
impl SelectionKind for MutableSerial {}
impl MutableSel for MutableSerial {}
impl SerialSel for MutableSerial {}
impl AllowsSubselect for MutableSerial {}

/// Marker type for possibly overlapping builder selection (single-threaded)
pub struct BuilderSerial(PhantomData<*const ()>);
impl SelectionKind for BuilderSerial {
    #[inline(always)]
    fn check_index(index: &SortedSet<usize>, topology: &Topology, state: &State) -> Result<(), SelectionError> {
        let first = index[0];
        let last = index[index.len()-1];
        let n1 = topology.num_atoms();
        let n2 = state.num_coords();
        if first >= n1 || last >= n1 ||  first >= n2 || last >= n2 {
            return Err(SelectionError::IndexCheck(first,last,n1));
        } else {
            Ok(())
        }
    }    
}
impl MutableSel for BuilderSerial {}
impl SerialSel for BuilderSerial {}
impl AllowsSubselect for BuilderSerial {}

/// Marker type for non-overlapping mutable selection (multi-threaded)
pub struct MutableParallel {}

impl SelectionKind for MutableParallel {
    #[inline(always)]
    fn check_overlap(index: &SortedSet<usize>, used: &mut rustc_hash::FxHashSet<usize>) -> Result<(), SelectionError> {
        for i in index.iter() {
            if !used.insert(*i) {
                return Err(SelectionError::OverlapCheck(*i));
            }
        }
        Ok(())
    }
}
impl MutableSel for MutableParallel {}
impl ParallelSel for MutableParallel {}
// Doesn't allow subselecting!

/// Marker type for possibly overlapping immutable selection (multi-threaded)
pub struct ImmutableParallel {}
impl SelectionKind for ImmutableParallel {}
impl ParallelSel for ImmutableParallel {}
impl AllowsSubselect for ImmutableParallel {}
