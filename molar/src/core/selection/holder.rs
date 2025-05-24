use std::{marker::PhantomData, ops::Deref};
use super::{SelectionError, SelectionKind, UserCreatableKind};

/// Smart pointer wrapper for sharing [Topology](crate::core::Topology) or [State](crate::core::State) between selections.
/// 
/// Acts like [Rc](std::rc::Rc) parameterized by [SelectionKind], so the user can't accidentally mix incompatible
/// selections pointing to the same data.
/// Can't be sent to other threads.
/// Normally this type should not be created directly by the user.
#[derive(Debug)]
pub struct Holder<T, K> {
    pub(super) arc: triomphe::Arc<T>,
    _kind: PhantomData<K>,
}

impl<T,K: SelectionKind> From<T> for Holder<T,K> {
    /// Create [Holder] from data
    fn from(value: T) -> Self {
        Self {
            arc: triomphe::Arc::new(value),
            _kind: Default::default(),
        }
    }
}

impl<T,K: UserCreatableKind> Holder<T, K> {
    /// Shallow-clones the [Holder] creating new reference to its data.
    pub fn new_ref(&self) -> Self {
        Self {
            arc: triomphe::Arc::clone(&self.arc),
            _kind: Default::default(),
        }
    }
}

impl<T,K: SelectionKind> Holder<T, K> {
    pub fn ref_count(&self) -> usize {
        triomphe::Arc::count(&self.arc)
    }

    /// Check if two holders point to the same data
    pub fn same_data(&self,other: &Self) -> bool {
        triomphe::Arc::ptr_eq(&self.arc, &other.arc)
    }

    // Clone to other kinds unsafely
    pub unsafe fn new_ref_with_kind<KO: SelectionKind>(&self) -> Holder<T,KO> {
        Holder {
            arc: triomphe::Arc::clone(&self.arc),
            _kind: Default::default(),
        }
    }

    /// Release a wrapped type if it is uniquesly owned
    pub fn release(self) -> Result<T, SelectionError> {
        Ok( triomphe::Arc::try_unwrap(self.arc).map_err(|_| SelectionError::Release)? )
    }
}

/// Holders are dereferenced as usual smart pointers
impl<T, K: SelectionKind> Deref for Holder<T, K> {
    type Target = T;
    fn deref(&self) -> &Self::Target {
        &self.arc
    }
}