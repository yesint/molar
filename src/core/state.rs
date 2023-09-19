use super::PeriodicBox;

#[derive(Debug, Default)]
pub struct State {
    pub coords: Vec<[f32; 3]>,
    pub time: f32,
    pub box_: PeriodicBox, 
}

impl State {
    pub fn new() -> Self {
        Default::default()
    }
}
