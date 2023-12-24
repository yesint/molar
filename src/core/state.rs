use std::{rc::Rc, cell::RefCell, sync::{RwLock, Arc}};

use crate::io::IndexAndStateProvider;

use super::{PeriodicBox, Pos};
//use super::handle::{SharedHandle, Handle};

#[derive(Debug, Default,Clone)]
pub struct State {
    pub coords: Vec<Pos>,
    pub time: f32,
    pub box_: Option<PeriodicBox>,
}

impl State {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn to_rc(self) -> Rc<RefCell<Self>> {
        Rc::new(RefCell::new(self))
    }

    pub fn to_arc(self) -> Arc<RwLock<Self>> {
        Arc::new(RwLock::new(self))
    }
}

impl IndexAndStateProvider for State {
    fn get_index_and_state(&self) -> (impl super::IndexIterator, &State) {
        (0..self.coords.len(), &self)
    }
}
