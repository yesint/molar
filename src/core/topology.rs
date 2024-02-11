use std::{rc::Rc, cell::RefCell, sync::{Arc, RwLock}, ops::Deref};

use crate::io::{IoIndexProvider, IoTopologyProvider};

use super::Atom;
//use super::handle::{SharedHandle, Handle};

#[allow(dead_code)]
#[derive(Debug, Default, Clone)]
pub struct Topology {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<[usize; 2]>,
    pub molecules: Vec<[usize; 2]>,
}

pub type TopologyRc = Rc<RefCell<Topology>>;

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

impl IoIndexProvider for Topology {
    fn get_index(&self) -> impl super::IndexIterator {
        0..self.atoms.len()
    }
}

impl IoTopologyProvider for Topology {
    #[allow(refining_impl_trait)]
    fn get_topology(&self) -> &Topology {
        self
    }
}

impl IoIndexProvider for TopologyRc {
    fn get_index(&self) -> impl super::IndexIterator {
        0..self.borrow().atoms.len()
    }
}

impl IoTopologyProvider for TopologyRc {
    fn get_topology(&self) -> impl Deref<Target = Topology> {
        self.borrow()
    }
}