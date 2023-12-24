use std::{rc::Rc, cell::RefCell, sync::{RwLock, Arc}};

use crate::io::{IndexProvider, StateProvider};

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

impl IndexProvider for State {
    fn get_index(&self) -> impl super::IndexIterator {
        0..self.coords.len()
    }
}

impl StateProvider for State {
    fn get_state(&self) -> &State {
        self
    }
}
