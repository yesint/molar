
use std::{sync::{Arc,RwLock, RwLockReadGuard, RwLockWriteGuard}, rc::Rc, cell::{RefCell, Ref, RefMut}, ops::Deref};

use super::State;

//---------------------------------------------------------
#[derive(Debug,Clone)]
pub struct Handle<T>(
    Rc<RefCell<T>>,
);

impl<T> Handle<T> {
    pub fn read(&self) -> Ref<'_,T> {
        self.0.borrow()
    }

    pub fn write(&self) -> RefMut<'_,T> {
        self.0.borrow_mut()
    }
}

impl<T> From<T> for Handle<T> {
    fn from(value: T) -> Self {
        Self(Rc::new(RefCell::new(value)))
    }
}

//---------------------------------------------------------
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