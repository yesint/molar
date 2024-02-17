use std::{cell::{Ref, RefCell, RefMut}, rc::Rc, sync::{Arc, RwLock}};
use crate::io::{IoIndexProvider, IoStateProvider};
use anyhow::{Result, anyhow};
use super::{BoxProvider, IdPosIterator, IndexIterator, Measure, PeriodicBox, Pos, PosProvider};
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

impl Measure for StateRc {
    type Provider<'a> = Ref<'a,State>;
    fn get_provider<'a>(&'a self) -> Self::Provider<'a> {
        self.borrow()
    }
}

impl IoIndexProvider for Ref<'_,State> {
    fn get_index(&self) -> impl super::IndexIterator {
        0..self.coords.len()
    }
}

impl IoStateProvider for Ref<'_,State> {
    #[allow(refining_impl_trait)]
    fn get_state(&self) -> &State {
        self
    }
}

impl PosProvider for Ref<'_,State> {
    fn iter_pos(&self) -> impl super::PosIterator<'_> {
        self.coords.iter()
    }
}

impl BoxProvider for Ref<'_,State> {
    fn get_box(&self) -> Result<&PeriodicBox> {
        let r = self.box_
            .as_ref()
            .ok_or(anyhow!("No periodic box"))?;
        Ok(&r)
    }
}

impl BoxProvider for RefMut<'_,State> {
    fn get_box(&self) -> Result<&PeriodicBox> {
        let r = self.box_
            .as_ref()
            .ok_or(anyhow!("No periodic box"))?;
        Ok(&r)
    }
}
