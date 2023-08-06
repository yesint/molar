use nalgebra::Matrix3;

#[derive(Debug, Default)]
pub struct State {
    pub coords: Vec<[f32; 3]>,
    pub time: f32,
    pub box_: Matrix3<f32>, 
}

impl State {
    pub fn new() -> Self {
        Default::default()
    }
}
