use std::marker::PhantomData;
use sorted_vec::SortedSet;
use crate::prelude::*;
use super::SelectionError;

/// Trait for kinds of selections
pub trait SelectionKind {
    #[inline(always)]
    #[allow(unused_variables)]
    fn check_index(index: &SortedSet<usize>, topology: &Topology, state: &State) -> Result<(), SelectionError> {
        Ok(())
    }
}

/// Trait marking mutable selections
pub trait MutableKind: SelectionKind {}
/// Trait marking selections that can be created by user
pub trait UserCreatableKind: SelectionKind {}

//------------------------------------------------------

/// Marker type for possibly overlapping mutable selection (single-threaded)
pub struct MutableSerial(PhantomData<*const ()>);

impl SelectionKind for MutableSerial {}
impl MutableKind for MutableSerial {}
impl UserCreatableKind for MutableSerial {}

/// Marker type for possibly overlapping builder selection (single-threaded)
pub struct BuilderSerial(PhantomData<*const ()>);

impl SelectionKind for BuilderSerial {
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

impl MutableKind for BuilderSerial {}
impl UserCreatableKind for BuilderSerial {}

/// Marker type for non-overlapping mutable selection (multi-threaded)
pub struct MutableParallel {}

impl SelectionKind for MutableParallel {}
impl MutableKind for MutableParallel {}
// Not user creatable!

/// Marker type for possibly overlapping immutable selection (multi-threaded)
pub struct ImmutableParallel {}

impl SelectionKind for ImmutableParallel {}
impl UserCreatableKind for ImmutableParallel {}