
use std::sync::{Arc,RwLock, RwLockReadGuard, RwLockWriteGuard};

#[derive(Debug,Clone)]
pub struct SharedHandle<T>(
    Arc<RwLock<T>>
);

impl<T> SharedHandle<T> {
    pub fn read(&self) -> RwLockReadGuard<'_,T> {
        self.0.read().expect("SharedHandle poisoned")
    }

    pub fn write(&self) -> RwLockWriteGuard<'_,T> {
        self.0.write().expect("SharedHandle poisoned")
    }
}

impl<T> From<T> for SharedHandle<T> {
    fn from(value: T) -> Self {
        Self(Arc::new(RwLock::new(value)))
    }
}