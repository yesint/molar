use super::{PeriodicBox, Pos};
use super::shared_handle::SharedHandle;

#[derive(Debug, Default)]
pub struct State {
    pub coords: Vec<Pos>,
    pub time: f32,
    pub box_: Option<PeriodicBox>,
}

impl State {
    pub fn new() -> Self {
        Default::default()
    }
}

pub type StateHandle = SharedHandle<State>;

