use std::{rc::Rc, cell::RefCell, sync::{RwLock, Arc}};

use crate::io::IndexAndStateProvider;
use anyhow::{Result, anyhow};
use super::{PeriodicBox, Pos, BoxProvider, Measure, MeasurePos};
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

impl BoxProvider for State {
    fn get_box(&self) -> Result<&PeriodicBox> {
        let r = self.box_
            .as_ref()
            .ok_or(anyhow!("No periodic box"))?;
        Ok(&r)
    }
}

impl MeasurePos for State {
    fn iter_pos(&self) -> impl super::PosIterator<'_> {
        self.coords.iter()
    }
}