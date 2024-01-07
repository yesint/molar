use std::{rc::Rc, cell::RefCell, sync::{Arc, RwLock}, ops::Deref};

use crate::io::IoIndexAndTopologyProvider;

use super::{Atom, IndexIterator};
//use super::handle::{SharedHandle, Handle};

#[allow(dead_code)]
#[derive(Debug, Default, Clone)]
pub struct Topology {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<[usize; 2]>,
    pub molecules: Vec<[usize; 2]>,
}

impl Topology {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn to_rc(self) -> Rc<RefCell<Self>> {
        Rc::new(RefCell::new(self))
    }

    pub fn to_arc(self) -> Arc<RwLock<Self>> {
        Arc::new(RwLock::new(self))
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

impl IoIndexAndTopologyProvider for Topology {
    fn get_index_and_topology(&self) -> (impl IndexIterator, impl Deref<Target=Topology>) {
        (0..self.atoms.len(), self)
    }
}

impl IoIndexAndTopologyProvider for Rc<RefCell<Topology>> {
    fn get_index_and_topology(&self) -> (impl IndexIterator, impl Deref<Target=Topology>) {
        (0..self.borrow().atoms.len(), self.borrow())
    }
}
