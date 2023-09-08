use super::Atom;

#[allow(dead_code)]
#[derive(Debug, Default, Clone)]
pub struct Structure {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<[usize; 2]>,
    pub molecules: Vec<[usize; 2]>,
}

impl Structure {
    pub fn new() -> Self {
        Default::default()
    }
}
