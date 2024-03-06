use std::{cell::{Ref, RefCell, RefMut}, rc::Rc};
use crate::io::StateProvider;
use super::{measure::GuardedQuery, providers::{BoxProvider, PosProvider}, PeriodicBox, Pos};
//use super::handle::{SharedHandle, Handle};

#[derive(Debug, Default,Clone)]
pub struct State {
    pub coords: Vec<Pos>,
    pub time: f32,
    pub pbox: Option<PeriodicBox>,
}

pub type StateRc = Rc<RefCell<State>>;

impl State {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn to_rc(self) -> Rc<RefCell<Self>> {
        Rc::new(RefCell::new(self))
    }
}

impl GuardedQuery for StateRc {
    type Guard<'a> = Ref<'a,State>;
    
    #[inline(always)]
    fn guard<'a>(&'a self) -> Self::Guard<'a> {
        self.borrow()
    }
}

impl StateProvider for Ref<'_,State> {
    fn get_time(&self) -> f32 {
        self.time
    }

    fn num_coords(&self) -> usize {
        self.coords.len()
    }
}

impl PosProvider for Ref<'_,State> {
    fn iter_pos(&self) -> impl super::PosIterator<'_> {
        self.coords.iter()
    }
}

impl BoxProvider for Ref<'_,State> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.pbox.as_ref()
    }
}

impl BoxProvider for RefMut<'_,State> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.pbox.as_ref()
    }
}
