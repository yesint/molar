use super::Atom;

#[allow(dead_code)]
#[derive(Debug, Default, Clone)]
pub struct Structure {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<[usize; 2]>,
    pub angles: Vec<[usize; 3]>,
    pub dihedrals: Vec<[usize; 4]>,
}

impl Structure {
    pub fn new() -> Self {
        Default::default()
    }
}
