use ascii::{AsciiChar, AsciiString};


#[allow(dead_code)]
#[derive(Debug, Default, Clone)]
pub struct Atom {
    pub name: AsciiString,
    pub resname: AsciiString,
    pub resid: isize,
    //resindex: usize,
    pub mass: f32,
    pub charge: f32,
    pub chain: AsciiChar,
}

impl Atom {
    pub fn new() -> Self {
        Default::default()
    }
}
