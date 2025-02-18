use std::{marker::PhantomData, ops::Deref};
use super::{BuilderSerial, ImmutableParallel, MutableParallel, MutableSerial, SelectionError, SelectionKind};

/// Smart pointer wrapper for sharing [Topology] and [State] between serial selections.
/// Acts like Rc parameterized by selection kind, so the user can't accidentally mix incompatible
/// selections pointing to the same data.
/// Can't be sent to other threads.
/// Normally this type should not be used directly by the user.
pub struct Holder<T, K> {
    pub(super) arc: triomphe::Arc<T>,
    _kind: PhantomData<K>,
}

// Conversion from T to Holder<T,K>
macro_rules! impl_from_t_for_holder {
    ( $k:ty ) => {
        impl<T> From<T> for Holder<T, $k> {
            fn from(value: T) -> Self {
                Self {
                    arc: triomphe::Arc::new(value),
                    _kind: Default::default(),
                }
            }
        }
    }
}

macro_rules! impl_clone_for_holder {
    ( $k:ty ) => {
        impl<T> Clone for Holder<T, $k> {
            fn clone(&self) -> Self {
                Self {
                    arc: self.arc.clone(),
                    _kind: Default::default(),
                }
            }
        }
    };
}

impl_from_t_for_holder!(MutableSerial);
impl_from_t_for_holder!(BuilderSerial);
impl_from_t_for_holder!(MutableParallel);
impl_from_t_for_holder!(ImmutableParallel);

impl_clone_for_holder!(MutableSerial);
impl_clone_for_holder!(BuilderSerial);
impl_clone_for_holder!(ImmutableParallel);


impl<T,K: SelectionKind> Holder<T, K> {
    /// Create new [Holder]
    pub fn new(value: T) -> Self {
        Self {
            arc: triomphe::Arc::new(value),
            _kind: Default::default(),
        }
    }

    /// Check if two holders point to the same data
    pub fn same_data(&self,other: &Self) -> bool {
        triomphe::Arc::ptr_eq(&self.arc, &other.arc)
    }

    // Clone has limited visibility and can clone to other kinds
    pub(super) fn clone_with_kind<KO: SelectionKind>(&self) -> Holder<T,KO> {
        Holder {
            arc: self.arc.clone(),
            _kind: Default::default(),
        }
    }

    /// Unsafely swaps allocations of two holders without any checks
    pub unsafe fn swap_unchecked(
        &mut self,
        other: &mut Self,
    ) {
        let p1 = self.arc.as_ptr() as *mut T;
        let p2 = other.arc.as_ptr() as *mut T;
        // We physically spap memory at these locations
        // this is cheap because coordinates are allocated on heap
        // and only pointers to allocations are swapped
        std::ptr::swap(p1, p2);
    }

    /// Replace arc in this Holder without any checks
    pub unsafe fn replace_arc<KO: SelectionKind>(&mut self, other: Holder<T,KO>) {
        self.arc = other.arc;
    }

}

/// Holders are dereferenced as usual smart pointers
impl<T, K: SelectionKind> Deref for Holder<T, K> {
    type Target = T;
    fn deref(&self) -> &Self::Target {
        &self.arc
    }
}

impl<T: Clone, K: SelectionKind> Holder<T, K> {
    /// Holder can release a wrapped type if it is uniquesly owned
    pub fn release(self) -> Result<T, SelectionError> {
        Ok( triomphe::Arc::try_unwrap(self.arc).map_err(|_| SelectionError::Release)? )
    }
}