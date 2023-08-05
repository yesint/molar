
#[derive(Debug, Default)]
pub struct State {
    pub coords: Vec<[f32; 3]>,
}

impl State {
    pub fn new() -> Self {
        Default::default()
    }
}
