use std::marker::PhantomData;
use sorted_vec::SortedSet;
use crate::io::StateProvider;

/// Trait for kinds of selections
pub trait SelectionKind {
    type SubselKind: SelectionKind;
    const NEED_CHECK_OVERLAP: bool;

    #[inline(always)]
    #[allow(unused_variables)]
    fn check_index(index: &SortedSet<usize>, system: &super::System) {}
}

/// Trait marking selections that can overlap
pub trait MayOverlap: SelectionKind {}

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
    const NEED_CHECK_OVERLAP: bool = false;
}
impl MayOverlap for MutableSerial {}
impl MutableSel for MutableSerial {}
impl SerialSel for MutableSerial {}

/// Marker type for possibly overlapping builder selection (single-threaded)
pub struct BuilderSerial(PhantomData<*const ()>);
impl SelectionKind for BuilderSerial {
    type SubselKind = BuilderSerial;
    const NEED_CHECK_OVERLAP: bool = false;
    #[inline(always)]
    fn check_index(index: &SortedSet<usize>, system: &super::System) {
        let first = index[0];
        let last = index[index.len()-1];
        let n = system.state.num_coords();
        if first >= n || last >= n {
            panic!(
                "Builder selection indexes [{}:{}] are out of allowed range [0:{}]",
                first,last,n
            );
        }
    }
}
impl MayOverlap for BuilderSerial {}
impl MutableSel for BuilderSerial {}
impl SerialSel for BuilderSerial {}

/// Marker type for non-overlapping mutable selection (multi-threaded)
pub struct MutableParallel {}

impl SelectionKind for MutableParallel {
    // Subseletions may overlap but won't be Send
    type SubselKind = MutableSerial;
    const NEED_CHECK_OVERLAP: bool = true;
}
impl MutableSel for MutableParallel {}
impl ParallelSel for MutableParallel {}

/// Marker type for possibly overlapping immutable selection (multi-threaded)
pub struct ImmutableParallel {}
impl SelectionKind for ImmutableParallel {
    type SubselKind = ImmutableParallel;
    const NEED_CHECK_OVERLAP: bool = false;
}
impl MayOverlap for ImmutableParallel {}
impl ParallelSel for ImmutableParallel {}
