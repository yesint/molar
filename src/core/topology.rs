use std::{cell::{Ref, RefCell, RefMut}, rc::Rc};

use crate::io::{IoIndexProvider, IoTopologyProvider};

use super::{providers::{AtomsMutProvider, AtomsProvider, MassesProvider}, Atom, modify::GuardedModify, measure::GuardedQuery};
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

impl GuardedModify for TopologyRc {
    type GuardMut<'a> = RefMut<'a,Topology>;
    fn guard_mut<'a>(&'a self) -> Self::GuardMut<'a> {
        self.borrow_mut()
    }
}

impl IoIndexProvider for Ref<'_,Topology> {
    fn get_index(&self) -> impl super::IndexIterator {
        0..self.atoms.len()
    }
}

impl IoTopologyProvider for Ref<'_,Topology> {
    fn get_topology(&self) -> &Topology {
        self
    }
}

impl AtomsProvider for Ref<'_,Topology> {
    fn iter_atoms(&self) -> impl super::AtomIterator<'_> {
        self.atoms.iter()
    }
}

impl AtomsMutProvider for RefMut<'_,Topology> {
    fn iter_atoms_mut(&mut self) -> impl super::AtomMutIterator<'_> {
        self.atoms.iter_mut()
    }
}

impl MassesProvider for Ref<'_,Topology> {
    fn iter_masses(&self) -> impl ExactSizeIterator<Item = f32> {
        self.atoms.iter().map(|at| at.mass)
    }
}