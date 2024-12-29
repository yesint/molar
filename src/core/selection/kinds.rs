use std::{marker::PhantomData, sync::Mutex};
use sorted_vec::SortedSet;
use crate::{core::{State, Topology}, io::{IndexProvider, StateProvider, TopologyProvider}};

use super::{holder::HolderOverlapCheckError, SelectionError};

/// Trait for kinds of selections
pub trait SelectionKind {
    type UsedIndexesType: Clone + Default;

    #[inline(always)]
    #[allow(unused_variables)]
    fn check_index(index: &SortedSet<usize>, topology: &Topology, state: &State) -> Result<(), SelectionError> {
        Ok(())
    }

    #[inline(always)]
    #[allow(unused_variables)]
    fn try_add_to_used(index: &impl IndexProvider, used: &Self::UsedIndexesType) -> Result<(), HolderOverlapCheckError> {
        Ok(())
    }

    #[inline(always)]
    #[allow(unused_variables)]
    fn remove_used(index: &impl IndexProvider, used: &Self::UsedIndexesType) {
        ()
    }
}

/// Trait marking mutable selections
pub trait MutableSel: SelectionKind {}

//------------------------------------------------------

/// Marker type for possibly overlapping mutable selection (single-threaded)
pub struct MutableSerial(PhantomData<*const ()>);

impl SelectionKind for MutableSerial {
    type UsedIndexesType = ();
}

impl MutableSel for MutableSerial {}

/// Marker type for possibly overlapping builder selection (single-threaded)
pub struct BuilderSerial(PhantomData<*const ()>);

impl SelectionKind for BuilderSerial {
    type UsedIndexesType = ();

    #[inline(always)]
    fn check_index(index: &SortedSet<usize>, topology: &Topology, state: &State) -> Result<(), SelectionError> {
        let first = unsafe { *index.get_unchecked(0) };
        let last = unsafe { *index.get_unchecked(index.len()-1) };
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

/// Marker type for non-overlapping mutable selection (multi-threaded)
pub struct MutableParallel {}

impl SelectionKind for MutableParallel {
    type UsedIndexesType = triomphe::Arc<Mutex<rustc_hash::FxHashSet<usize>>>;

    #[inline(always)]
    fn try_add_to_used(index: &impl IndexProvider, used: &Self::UsedIndexesType) -> Result<(), HolderOverlapCheckError> {
        let mut g = used.lock().unwrap();
        for i in index.iter_index() {
            if g.contains(&i) {
                return Err(HolderOverlapCheckError(i));
            }
        }
        // If all indexes are clear and not used add them
        for i in index.iter_index() {
            g.insert(i);
        }
        Ok(())
    }

    fn remove_used(index: &impl IndexProvider, used: &Self::UsedIndexesType) {
        // Obtain a lock for used
        let mut g = used.lock().unwrap();
        for i in index.iter_index() {
            g.remove(&i);
        }
    }
}

impl MutableSel for MutableParallel {}
// Doesn't allow subselecting!

/// Marker type for possibly overlapping immutable selection (multi-threaded)
pub struct ImmutableParallel {}

impl SelectionKind for ImmutableParallel {
    type UsedIndexesType = ();
}