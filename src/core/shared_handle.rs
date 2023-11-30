
use std::sync::{Arc,RwLock, RwLockReadGuard, RwLockWriteGuard};
use anyhow::bail;

#[derive(Debug,Clone)]
pub struct SharedHandle<T>(
    Arc<RwLock<T>>
);

impl<T> SharedHandle<T> {
    pub fn read(&self) -> anyhow::Result<RwLockReadGuard<'_,T>> {
        match self.0.read() {
            Ok(h) => Ok(h),
            Err(e) => bail!(e.to_string()),
        }
    }

    pub fn write(&self) -> anyhow::Result<RwLockWriteGuard<'_,T>> {
        match self.0.write() {
            Ok(h) => Ok(h),
            Err(e) => bail!(e.to_string()),
        }
    }
}

impl<T> From<T> for SharedHandle<T> {
    fn from(value: T) -> Self {
        Self(Arc::new(RwLock::new(value)))
    }
}