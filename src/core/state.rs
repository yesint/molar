use super::PeriodicBox;
use nalgebra::Point3;

#[derive(Debug, Default)]
pub struct State {
    pub coords: Vec<Point3<f32>>,
    pub time: f32,
    pub box_: PeriodicBox,
}

impl State {
    pub fn new() -> Self {
        Default::default()
    }
}

