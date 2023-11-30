use super::Atom;
use super::shared_handle::SharedHandle;

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

    pub fn assign_resindex(&mut self) {
        let mut resindex = 0usize;
        let mut cur_resid = self.atoms[0].resid;
        for at in self.atoms.iter_mut() {
            if at.resid != cur_resid {
                cur_resid = at.resid;
                resindex += 1;
            }
            at.resindex = resindex;
        }
    }
}

pub type StructureHandle = SharedHandle<Structure>;