use std::marker::PhantomData;

pub(crate) trait CheckedIndex {
    fn check_index(&self) -> &Vec<usize>;
}

/// Trait for kinds of selections
pub trait SelectionKind {
    type SubselKind: SelectionKind;
    const NEED_CHECK_OVERLAP: bool;
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
