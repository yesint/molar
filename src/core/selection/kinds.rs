use std::marker::PhantomData;
use sorted_vec::SortedSet;
use crate::{core::{State, Topology}, io::{IndexProvider, StateProvider, TopologyProvider}};

use super::SelectionError;

/// Trait for kinds of selections
pub trait SelectionKind {
    type UsedIndexType: Clone;

    #[inline(always)]
    #[allow(unused_variables)]
    fn check_index(index: &SortedSet<usize>, system: &super::System) -> Result<(), SelectionError> {
        Ok(())
    }

    #[inline(always)]
    #[allow(unused_variables)]
    fn try_add_used(index: &impl IndexProvider, used: &Self::UsedIndexType) -> Result<(), SelectionError> {
        Ok(())
    }

    #[inline(always)]
    #[allow(unused_variables)]
    fn remove_used(index: &impl IndexProvider, used: &Self::UsedIndexType) {
        ()
    }
}

/// Trait marking non-overlapping selections
pub trait MutableSel: SelectionKind {}
/// Trait marking parallel selections
pub trait ParallelSel: SelectionKind + Send + Sync {}
/// Trait marking serial selections
pub trait SerialSel: SelectionKind {}


/// Marker type for possibly overlapping mutable selection (single-threaded)
pub struct MutableSerial(PhantomData<*const ()>);
impl SelectionKind for MutableSerial {
    type UsedIndexType = ();
}
impl MutableSel for MutableSerial {}
impl SerialSel for MutableSerial {}

/// Marker type for possibly overlapping builder selection (single-threaded)
pub struct BuilderSerial(PhantomData<*const ()>);
impl SelectionKind for BuilderSerial {
    type UsedIndexType = ();

    #[inline(always)]
    fn check_index(index: &SortedSet<usize>, system: &super::System) -> Result<(), SelectionError> {
        let first = index[0];
        let last = index[index.len()-1];
        let n = system.state.num_coords();
        if first >= n || last >= n {
            Err(SelectionError::IndexCheck(first,last,n))
        } else {
            Ok(())
        }
    }    
}
impl MutableSel for BuilderSerial {}
impl SerialSel for BuilderSerial {}

/// Marker type for non-overlapping mutable selection (multi-threaded)
pub struct MutableParallel {}

impl SelectionKind for MutableParallel {
    type UsedIndexType = super::UsedHashMap;

    #[inline(always)]
    fn try_add_used(index: &impl IndexProvider, used: &Self::UsedIndexType) -> Result<(), SelectionError> {
        let mut g = used.lock().unwrap();
        for i in index.iter_index() {
            if g.contains(&i) {
                return Err(SelectionError::OverlapCheck(i));
            }
        }
        // If all indexes are clear and not used add them
        for i in index.iter_index() {
            g.insert(i);
        }
        Ok(())
    }

    fn remove_used(index: &impl IndexProvider, used: &Self::UsedIndexType) {
        // Obtain a lock for used
        let mut g = used.lock().unwrap();
        for i in index.iter_index() {
            g.remove(&i);
        }
    }
}
impl MutableSel for MutableParallel {}
impl ParallelSel for MutableParallel {}

/// Marker type for possibly overlapping immutable selection (multi-threaded)
pub struct ImmutableParallel {}
impl SelectionKind for ImmutableParallel {
    type UsedIndexType = ();
}
impl ParallelSel for ImmutableParallel {}
