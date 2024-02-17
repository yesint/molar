use std::{cell::{Ref, RefCell}, rc::Rc};

use crate::io::{IoIndexProvider, IoTopologyProvider};

use super::{Atom, GuardedQuery};
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

impl GuardedQuery for TopologyRc {
    type Guard<'a> = Ref<'a,Topology>;
    fn guard<'a>(&'a self) -> Self::Guard<'a> {
        self.borrow()
    }
}

impl IoIndexProvider for Ref<'_,Topology> {
    fn get_index(&self) -> impl super::IndexIterator {
        0..self.atoms.len()
    }
}

impl IoTopologyProvider for Ref<'_,Topology> {
    #[allow(refining_impl_trait)]
    fn get_topology(&self) -> &Topology {
        self
    }
}