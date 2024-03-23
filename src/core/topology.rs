use std::ops::Deref;
use sync_unsafe_cell::SyncUnsafeCell;

use crate::io::TopologyProvider;

use super::{providers::{AtomsMutProvider, AtomsProvider, MassesProvider}, Atom};
//use super::handle::{SharedHandle, Handle};

#[allow(dead_code)]
#[derive(Debug, Default, Clone)]
pub struct TopologyStorage {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<[usize; 2]>,
    pub molecules: Vec<[usize; 2]>,
}

pub struct Topology(SyncUnsafeCell<TopologyStorage>);

impl Clone for Topology {
    fn clone(&self) -> Self {
        Self(SyncUnsafeCell::new(self.get().clone()))
    }
}

impl From<TopologyStorage> for Topology {
    fn from(value: TopologyStorage) -> Self {
        Self(SyncUnsafeCell::new(value))
    }
}

impl Topology {
    // Private convenience accessors
    #[inline(always)]
    fn get(&self) -> &TopologyStorage {
        unsafe {&*self.0.get()}
    }

    #[inline(always)]
    fn get_mut(&self) -> &mut TopologyStorage {
        unsafe {&mut *self.0.get()}
    }

    //--------------------------------------------

    pub fn to_rc(self) -> triomphe::UniqueArc<Self> {
        triomphe::UniqueArc::new(self)
    }

    #[inline(always)]
    pub unsafe fn nth_atom_unchecked(&self, i: usize) -> &Atom {
        self.get().atoms.get_unchecked(i)
    }

    #[inline(always)]
    pub unsafe fn nth_atom_unchecked_mut(&self, i: usize) -> &mut Atom {
        self.get_mut().atoms.get_unchecked_mut(i)
    }

    #[inline(always)]
    pub fn nth_atom(&self, i: usize) -> Option<&Atom> {
        self.get().atoms.get(i)
    }

    #[inline(always)]
    pub fn nth_atom_mut(&self, i: usize) -> Option<&mut Atom> {
        self.get_mut().atoms.get_mut(i)
    }

    pub fn assign_resindex(&mut self) {
        let mut resindex = 0usize;
        let mut cur_resid = unsafe{self.nth_atom_unchecked(0)}.resid;
        for at in self.iter_atoms_mut() {
            if at.resid != cur_resid {
                cur_resid = at.resid;
                resindex += 1;
            }
            at.resindex = resindex;
        }
    }
}

// Impls for Topology itself
impl TopologyProvider for Topology {
    fn num_atoms(&self) -> usize {
        self.get().atoms.len()
    }
}

impl AtomsProvider for Topology {
    fn iter_atoms(&self) -> impl super::AtomIterator<'_> {
        self.get().atoms.iter()
    }
}

impl AtomsMutProvider for Topology {
    fn iter_atoms_mut(&self) -> impl super::AtomMutIterator<'_> {
        self.get_mut().atoms.iter_mut()
    }
}

impl MassesProvider for Topology {
    fn iter_masses(&self) -> impl ExactSizeIterator<Item = f32> {
        self.get().atoms.iter().map(|at| at.mass)
    }
}

//--------------------------
// Impls for smart pointers
//--------------------------
impl<T: Deref<Target=Topology>> TopologyProvider for T {
    fn num_atoms(&self) -> usize {
        self.get().atoms.len()
    }
}

impl<T: Deref<Target=Topology>> AtomsProvider for T {
    fn iter_atoms(&self) -> impl super::AtomIterator<'_> {
        self.get().atoms.iter()
    }
}

impl<T: Deref<Target=Topology>> AtomsMutProvider for T {
    fn iter_atoms_mut(&self) -> impl super::AtomMutIterator<'_> {
        self.get_mut().atoms.iter_mut()
    }
}

impl<T: Deref<Target=Topology>> MassesProvider for T {
    fn iter_masses(&self) -> impl ExactSizeIterator<Item = f32> {
        self.get().atoms.iter().map(|at| at.mass)
    }
}