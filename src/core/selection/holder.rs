use std::{marker::PhantomData, ops::Deref};
use crate::io::IndexProvider;
use thiserror::Error;
use super::{BuilderSerial, ImmutableParallel, MutableParallel, MutableSerial, SelectionError, SelectionKind};

/// Smart pointer wrapper for sharing [Topology] and [State] between serial selections.
/// Acts like Rc parameterized by selection kind, so the user can't accidentally mix incompatible
/// selections pointing to the same data.
/// Can't be sent to other threads.
/// Normally this type should not be used directly by the user.
pub struct Holder<T, K: SelectionKind> {
    arc: triomphe::Arc<T>,
    used: K::UsedIndexesType,
    _kind: PhantomData<K>,
}

// This should not fail because there is no legal way to make
// non-unique untyped holder.
// Just in case we check this and panic if this is violated.
macro_rules! impl_from_t_for_holder {
    ( $t:ty ) => {
        impl<T> From<T> for Holder<T, $t> {
            fn from(value: T) -> Self {
                Self {
                    arc: triomphe::Arc::new(value),
                    used: Default::default(),
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

// All holders except MutableParallel allowSelection cloning
impl<T,K> Clone for Holder<T,K> 
where 
    K: SelectionKind<UsedIndexesType = ()>
{
    fn clone(&self) -> Self {
        Self {
            arc: self.arc.clone(),
            used: (),
            _kind: Default::default(),
        }
    }
}

// All holders define clone_with_index function
// which just redirects to the function of marker type K
impl<T,K: SelectionKind> Holder<T, K> {
    pub(crate) fn clone_with_index(
        &self,
        ind: &impl IndexProvider,
    ) -> Result<Self, HolderOverlapCheckError> {
        K::try_add_used(ind, &self.used)?;
        Ok(Self {
            arc: self.arc.clone(),
            used: self.used.clone(),
            _kind: Default::default(),
        })
    }

    // Unsafe access to used indexes
    pub(crate) unsafe fn get_used(&self) -> &K::UsedIndexesType {
        &self.used
    }

    pub fn same_data(&self,other: &Self) -> bool {
        triomphe::Arc::ptr_eq(&self.arc, &other.arc)
    }
}

// All holders are dereferenced as usual smart pointers
impl<T, K: SelectionKind> Deref for Holder<T, K> {
    type Target = T;
    fn deref(&self) -> &Self::Target {
        &self.arc
    }
}

// All holders can release a type if
// they are uniquesly owned
impl<T: Clone, K: SelectionKind> Holder<T, K> {
    pub fn release(self) -> Result<T, SelectionError> {
        if self.arc.is_unique() {
            Ok(triomphe::Arc::unwrap_or_clone(self.arc))
        } else {
            Err(SelectionError::Release)
        }
    }
}

//----------------------------------------------------------

#[derive(Error,Debug)]
#[error("index {0} is already used")]
pub struct HolderOverlapCheckError(pub usize);
