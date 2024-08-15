use std::marker::PhantomData;
use sorted_vec::SortedSet;
use crate::{core::{State, Topology}, io::{StateProvider, TopologyProvider}};

use super::{SelectionError, TopologyStateSizes};

/// Trait for kinds of selections
pub trait SelectionKind {
    type SubselKind: SelectionKind;    

    #[inline(always)]
    #[allow(unused_variables)]
    fn check_index<K>(index: &SortedSet<usize>, top: &Topology<K>, st: &State<K>) -> Result<(), SelectionError> {
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


/// Marker type for possibly overlapping mutable selection (single-threaded)
pub struct MutableSerial(PhantomData<*const ()>);
impl SelectionKind for MutableSerial {
    type SubselKind = MutableSerial;
}
impl MutableSel for MutableSerial {}
impl SerialSel for MutableSerial {}

/// Marker type for possibly overlapping builder selection (single-threaded)
pub struct BuilderSerial(PhantomData<*const ()>);
impl SelectionKind for BuilderSerial {
    type SubselKind = BuilderSerial;    
    
    #[inline(always)]
    fn check_index<K>(index: &SortedSet<usize>, top: &Topology<K>, st: &State<K>) -> Result<(), SelectionError> {
        let first = index[0];
        let last = index[index.len()-1];
        //if st.num_coords() != top.num_atoms() {
        //    return Err(SelectionError::DifferentSizes(TopologyStateSizes(st.num_coords(),top.num_atoms())))
        //}
        let n1 = top.num_atoms();
        let n2 = st.num_coords();
        if first >= n1 || last >= n1 {
            return Err(SelectionError::IndexCheck(first,last,n1))
        } 
        if first >= n2 || last >= n2 {
            return Err(SelectionError::IndexCheck(first,last,n2))
        }
        Ok(())
    }    
}
impl MutableSel for BuilderSerial {}
impl SerialSel for BuilderSerial {}

/// Marker type for non-overlapping mutable selection (multi-threaded)
pub struct MutableParallel {}

impl SelectionKind for MutableParallel {
    // Subseletions may overlap but won't be Send
    type SubselKind = MutableSerial;    

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

/// Marker type for possibly overlapping immutable selection (multi-threaded)
pub struct ImmutableParallel {}
impl SelectionKind for ImmutableParallel {
    type SubselKind = ImmutableParallel;
}
impl ParallelSel for ImmutableParallel {}
