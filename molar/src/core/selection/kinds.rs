use std::marker::PhantomData;
use sorted_vec::SortedSet;
use crate::prelude::*;
use super::SelectionError;

/// Trait for kinds of selections
pub trait SelectionKind {
    #[inline(always)]
    #[allow(unused_variables)]
    fn validate_index(index: &SortedSet<usize>, topology: &Topology, state: &State) -> Result<(), SelectionError> {
        Ok(())
    }
}

/// Trait marking mutable selections
pub trait MutableKind: SelectionKind {}
/// Trait marking selections that can be created by user
pub trait UserCreatableKind: SelectionKind {}
pub trait SerialKind: UserCreatableKind {}

//---------------------------------------------------------------------------
/// Marker type for possibly overlapping mutable selection (single-threaded)
//---------------------------------------------------------------------------
pub struct MutableSerial(PhantomData<*const ()>);

impl SelectionKind for MutableSerial {}
impl MutableKind for MutableSerial {}
impl UserCreatableKind for MutableSerial {}
impl SerialKind for MutableSerial {}
//---------------------------------------------------------------------------
/// Marker type for possibly overlapping builder selection (single-threaded)
//---------------------------------------------------------------------------
pub struct BuilderSerial(PhantomData<*const ()>);

impl SelectionKind for BuilderSerial {
    #[inline(always)]
    fn validate_index(index: &SortedSet<usize>, topology: &Topology, state: &State) -> Result<(), SelectionError> {
        // Index is guaranteed to be non-empty
        let first = unsafe { *index.get_unchecked(0) };
        let last = unsafe { *index.get_unchecked(index.len()-1) };
        let ntop = topology.num_atoms();
        let nst = state.num_pos();
        if first >= ntop || last >= ntop ||  first >= nst || last >= nst {
            return Err(SelectionError::IndexValidation(first,last,ntop));
        } else {
            Ok(())
        }
    }    
}

impl MutableKind for BuilderSerial {}
impl UserCreatableKind for BuilderSerial {}
impl SerialKind for BuilderSerial {}

//---------------------------------------------------------------------------
/// Marker type for non-overlapping mutable selection (multi-threaded)
//---------------------------------------------------------------------------
pub struct MutableParallel {}

impl SelectionKind for MutableParallel {}
impl MutableKind for MutableParallel {}
// Not user creatable!

//---------------------------------------------------------------------------
/// Marker type for possibly overlapping immutable selection (multi-threaded)
//---------------------------------------------------------------------------
pub struct ImmutableParallel {}

impl SelectionKind for ImmutableParallel {}
impl UserCreatableKind for ImmutableParallel {}