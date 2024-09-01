use sync_unsafe_cell::SyncUnsafeCell;
use thiserror::Error;
use triomphe::{Arc, UniqueArc};
use crate::io::TopologyProvider;
use super::{providers::{AtomsMutProvider, AtomsProvider, MassesProvider}, Atom};

#[doc(hidden)]
#[derive(Debug, Default, Clone)]
pub(crate) struct TopologyStorage {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<[usize; 2]>,
    pub molecules: Vec<[usize; 2]>,
}

#[derive(Error,Debug)]
pub enum BuilderError {
    #[error("indexes to remove {0}:{1} are out of allowed range 0:{2}")]
    RemoveIndexes(usize,usize,usize)
}

impl TopologyStorage {
    pub(crate) fn add_atoms<'a>(&'a mut self, atoms: impl super::AtomIterator<'a>) {
        self.atoms.extend(atoms.cloned());
    }

    pub(crate) fn add_bonds(&mut self, added: impl Iterator<Item = [usize; 2]>) {
        todo!("Check if bond already exists");
        self.bonds.extend(added);
    }

    pub(crate) fn add_molecules(&mut self, added: impl Iterator<Item = [usize; 2]>) {
        self.molecules.extend(added);
    }

    pub(crate) fn remove_atoms(&mut self, removed: impl Iterator<Item = usize>) -> Result<(),BuilderError> {
        let mut ind = removed.collect::<Vec<_>>();
        if ind.len()==0 {
            return Ok(());
        }
        ind.sort_unstable();
        ind.dedup();
        if ind[0] >= self.atoms.len() || ind[ind.len()-1] >= self.atoms.len() {
            return Err(BuilderError::RemoveIndexes(ind[0],ind[ind.len()-1],self.atoms.len()));
        }
        for i in ind.iter().rev().cloned() {
            self.atoms.remove(i);
            // Remove affected bonds
            self.bonds.retain(|b| b[0] != i && b[1] != i);
            // Modify molecules and remove those, which become invalid
            self.molecules.retain_mut(|m| {
                if m[0]==i {m[0] += 1}
                if m[1]==i {
                    if m[1]>0 {
                        m[0] -= 1
                    } else {
                        // Remove whole molecule
                        return false
                    }
                }
                m[1]>=m[0]
            });
        }
        Ok(())
    }
}

/// Topology of the molecular system: atoms, bonds, molecules, etc.
/// 
/// [Topology] is typically read from structure of trajectory file and is not intended
/// to be manipulated directly by the user. Insead [State](super::State) and [Topology] 
/// are used to create atom selections, which give an access to the properties of
/// individual atoms and allow to query various properties.
#[derive(Default)]
pub struct Topology(SyncUnsafeCell<TopologyStorage>);

impl Clone for Topology {
    fn clone(&self) -> Self {
        Self(SyncUnsafeCell::new(self.get_storage().clone()))
    }
}

impl From<TopologyStorage> for UniqueArc<Topology> {
    fn from(value: TopologyStorage) -> Self {
        UniqueArc::new(Topology(SyncUnsafeCell::new(value)))
    }
}

impl Topology {
    // Private convenience accessors
    #[inline(always)]
    pub(super) fn get_storage(&self) -> &TopologyStorage {
        unsafe {&*self.0.get()}
    }

    #[inline(always)]
    pub(super) fn get_storage_mut(&self) -> &mut TopologyStorage {
        unsafe {&mut *self.0.get()}
    }

    #[inline(always)]
    pub unsafe fn nth_atom_unchecked(&self, i: usize) -> &Atom {
        self.get_storage().atoms.get_unchecked(i)
    }

    #[inline(always)]
    pub unsafe fn nth_atom_unchecked_mut(&self, i: usize) -> &mut Atom {
        self.get_storage_mut().atoms.get_unchecked_mut(i)
    }

    #[inline(always)]
    pub fn nth_atom(&self, i: usize) -> Option<&Atom> {
        self.get_storage().atoms.get(i)
    }

    #[inline(always)]
    pub fn nth_atom_mut(&self, i: usize) -> Option<&mut Atom> {
        self.get_storage_mut().atoms.get_mut(i)
    }

    pub fn assign_resindex(&self) {
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

    pub fn interchangeable(&self, other: &Topology) -> bool {
            self.get_storage().atoms.len() == other.get_storage().atoms.len()
        &&  self.get_storage().bonds.len() == other.get_storage().bonds.len()
        &&  self.get_storage().molecules.len() == other.get_storage().molecules.len()
    }
}

//---------------------------
macro_rules! impl_topology_traits {
    ( $t:ty ) => {
        impl TopologyProvider for $t {
            fn num_atoms(&self) -> usize {
                self.get_storage().atoms.len()
            }
        }
        
        impl AtomsProvider for $t {
            fn iter_atoms(&self) -> impl super::AtomIterator<'_> {
                self.get_storage().atoms.iter()
            }
        }
        
        impl AtomsMutProvider for $t {
            fn iter_atoms_mut(&self) -> impl super::AtomMutIterator<'_> {
                self.get_storage_mut().atoms.iter_mut()
            }
        }
        
        impl MassesProvider for $t {
            fn iter_masses(&self) -> impl ExactSizeIterator<Item = f32> {
                self.get_storage().atoms.iter().map(|at| at.mass)
            }
        }        
    }
}

// Impls for Topology itself
impl_topology_traits!(Topology);
// Impls for smart pointers
impl_topology_traits!(Arc<Topology>);
impl_topology_traits!(UniqueArc<Topology>);
