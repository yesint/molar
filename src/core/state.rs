
#[derive(Debug, Default)]
pub struct State {
    pub coords: Vec<[f32; 3]>,
    pub step: i32,
    pub time: f32,
    //pub box_: 
}

impl State {
    pub fn new() -> Self {
        Default::default()
    }
}
