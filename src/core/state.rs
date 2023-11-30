use super::{PeriodicBox, Pos};
use std::sync::{Arc,RwLock};
use anyhow::bail;

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

pub type StateHandle = Arc<RwLock<State>>;
