use std::{rc::Rc, cell::RefCell, sync::{RwLock, Arc}, ops::Deref};
use crate::io::{IoIndexProvider, IoStateProvider};
use anyhow::{Result, anyhow};
use super::{PeriodicBox, Pos, MeasureBox, MeasurePos, IdPosIterator, IndexIterator};
//use super::handle::{SharedHandle, Handle};

#[derive(Debug, Default,Clone)]
pub struct State {
    pub coords: Vec<Pos>,
    pub time: f32,
    pub box_: Option<PeriodicBox>,
}

pub type StateRc = Rc<RefCell<State>>;

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

    pub fn iter_id_pos_indexed<'a>(&'a self, index: impl IndexIterator) -> impl IdPosIterator<'a> {
        index.map(|i| (i,&self.coords[i]))
    }
}

impl IoIndexProvider for State {
    fn get_index(&self) -> impl super::IndexIterator {
        0..self.coords.len()
    }
}

impl IoStateProvider for State {
    fn get_state(&self) -> impl Deref<Target = State> {
        self
    }
}

impl IoIndexProvider for StateRc {
    fn get_index(&self) -> impl super::IndexIterator {
        0..self.borrow().coords.len()
    }
}

impl IoStateProvider for StateRc {
    fn get_state(&self) -> impl Deref<Target = State> {
        self.borrow()
    }
}

impl MeasureBox for State {
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