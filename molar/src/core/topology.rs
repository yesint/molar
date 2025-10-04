use std::sync::{RwLock, RwLockReadGuard, RwLockWriteGuard};

use crate::prelude::*;
use thiserror::Error;

#[doc(hidden)]
#[derive(Debug, Default, Clone)]
pub(crate) struct TopologyStorage {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<[usize; 2]>,
    pub molecules: Vec<[usize; 2]>,
}

/// Errors related to builder sources
#[derive(Error, Debug)]
pub enum BuilderError {
    #[error("indexes to remove {0}:{1} are out of allowed range 0:{2}")]
    RemoveIndexes(usize, usize, usize),
}

impl TopologyStorage {
    pub(crate) fn add_atoms<'a>(&'a mut self, atoms: impl Iterator<Item = Atom>) {
        self.atoms.extend(atoms);
    }

    pub(crate) fn remove_atoms(
        &mut self,
        removed: impl Iterator<Item = usize>,
    ) -> Result<(), BuilderError> {
        let mut ind = removed.collect::<Vec<_>>();
        if ind.len() == 0 {
            return Ok(());
        }
        ind.sort_unstable();
        ind.dedup();
        if ind[0] >= self.atoms.len() || ind[ind.len() - 1] >= self.atoms.len() {
            return Err(BuilderError::RemoveIndexes(
                ind[0],
                ind[ind.len() - 1],
                self.atoms.len(),
            ));
        }

        let mut it = ind.iter().cloned();
        let mut to_remove = it.next().unwrap_or(usize::MAX);
        let mut i = 0;
        self.atoms.retain(|_| {
            let ok = i != to_remove;
            i += 1;
            if !ok {
                to_remove = it.next().unwrap_or(usize::MAX);
            }
            ok
        });

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
pub struct Topology(RwLock<TopologyStorage>);

pub struct TopologyReadGuard<'a>(pub(crate) RwLockReadGuard<'a, TopologyStorage>);
pub struct TopologyWriteGuard<'a>(pub(crate) RwLockWriteGuard<'a, TopologyStorage>);

impl Topology {
    pub fn read(&self) -> TopologyReadGuard<'_> {
        TopologyReadGuard(self.0.read().unwrap())
    }

    pub fn write(&self) -> TopologyWriteGuard<'_> {
        TopologyWriteGuard(self.0.write().unwrap())
    }

    pub fn interchangeable(&self, other: &Topology) -> bool {
        let g1 = self.read();
        let g2 = other.read();
        g1.0.atoms.len() == g2.0.atoms.len()
            && g1.0.bonds.len() == g2.0.bonds.len()
            && g1.0.molecules.len() == g2.0.molecules.len()
    }
}

impl TopologyWriteGuard<'_> {
    pub fn assign_resindex(&mut self) {
        let mut resindex = 0usize;
        let mut cur_resid = unsafe { self.get_atom_unchecked(0) }.resid;
        for at in self.iter_atoms_mut() {
            if at.resid != cur_resid {
                cur_resid = at.resid;
                resindex += 1;
            }
            at.resindex = resindex;
        }
    }
}

//---------------------------
macro_rules! impl_topology_traits {
    ( $t:ty ) => {
        impl AtomIterProvider for $t {
            fn iter_atoms(&self) -> impl super::AtomIterator<'_> {
                self.0.atoms.iter()
            }
        }

        impl LenProvider for $t {
            fn len(&self) -> usize {
                self.0.atoms.len()
            }
        }

        impl RandomAtomProvider for $t {
            unsafe fn get_atom_unchecked(&self, i: usize) -> &Atom {
                self.0.atoms.get_unchecked(i)
            }
        }

        impl RandomBondProvider for $t {
            fn num_bonds(&self) -> usize {
                self.0.bonds.len()
            }

            unsafe fn get_bond_unchecked(&self, i: usize) -> &[usize; 2] {
                self.0.bonds.get_unchecked(i)
            }
        }

        impl BondIterProvider for $t {
            fn iter_bonds(&self) -> impl Iterator<Item = &[usize; 2]> {
                self.0.bonds.iter()
            }
        }

        impl RandomMoleculeProvider for $t {
            fn num_molecules(&self) -> usize {
                self.0.molecules.len()
            }

            unsafe fn get_molecule_unchecked(&self, i: usize) -> &[usize; 2] {
                self.0.molecules.get_unchecked(i)
            }
        }

        impl MoleculeIterProvider for $t {
            fn iter_molecules(&self) -> impl Iterator<Item = &[usize; 2]> {
                self.0.molecules.iter()
            }
        }

        impl TopologyWrite for $t {}
    };
}

macro_rules! impl_topology_mut_traits {
    ( $t:ty ) => {
        impl AtomIterMutProvider for $t {
            fn iter_atoms_mut(&mut self) -> impl super::AtomMutIterator<'_> {
                self.0.atoms.iter_mut()
            }
        }

        impl RandomAtomMutProvider for $t {
            fn get_atom_mut(&mut self, i: usize) -> Option<&mut Atom> {
                self.0.atoms.get_mut(i)
            }

            unsafe fn get_atom_mut_unchecked(&mut self, i: usize) -> &mut Atom {
                self.0.atoms.get_unchecked_mut(i)
            }
        }
    }
}

// Impls for Topology itself
impl_topology_traits!(TopologyReadGuard<'_>);
impl_topology_traits!(TopologyWriteGuard<'_>);
impl_topology_mut_traits!(TopologyWriteGuard<'_>);
// Impls for smart pointers
//impl_topology_traits!(triomphe::Arc<Topology>);
