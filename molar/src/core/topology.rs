use crate::prelude::*;
use sync_unsafe_cell::SyncUnsafeCell;
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

    // pub(crate) fn add_bonds(&mut self, added: impl Iterator<Item = [usize; 2]>) {
    //     todo!("Check if bond already exists");
    //     self.bonds.extend(added);
    // }

    // pub(crate) fn add_molecules(&mut self, added: impl Iterator<Item = [usize; 2]>) {
    //     self.molecules.extend(added);
    // }

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
pub struct Topology(SyncUnsafeCell<TopologyStorage>);

impl std::fmt::Debug for Topology {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = self.get_storage();
        f.debug_struct("Topology")
            .field("atoms", &s.atoms)
            .field("bonds", &s.bonds)
            .field("molecules", &s.molecules)
            .finish()
    }
}

impl Clone for Topology {
    fn clone(&self) -> Self {
        Self(SyncUnsafeCell::new(self.get_storage().clone()))
    }
}

impl From<TopologyStorage> for Topology {
    fn from(value: TopologyStorage) -> Self {
        Topology(SyncUnsafeCell::new(value))
    }
}

impl Topology {
    // Private convenience accessors
    #[inline(always)]
    pub(super) fn get_storage(&self) -> &TopologyStorage {
        unsafe { &*self.0.get() }
    }

    #[inline(always)]
    pub(super) fn get_storage_mut(&self) -> &mut TopologyStorage {
        unsafe { &mut *self.0.get() }
    }

    // #[inline(always)]
    // pub unsafe fn get_atom_unchecked(&self, i: usize) -> &Atom {
    //     self.get_storage().atoms.get_unchecked(i)
    // }

    // #[inline(always)]
    // pub unsafe fn get_atom_mut_unchecked(&self, i: usize) -> &mut Atom {
    //     self.get_storage_mut().atoms.get_unchecked_mut(i)
    // }

    // #[inline(always)]
    // pub fn get_atom(&self, i: usize) -> Option<&Atom> {
    //     self.get_storage().atoms.get(i)
    // }

    // #[inline(always)]
    // pub fn get_atom_mut(&self, i: usize) -> Option<&mut Atom> {
    //     self.get_storage_mut().atoms.get_mut(i)
    // }

    pub fn assign_resindex(&self) {
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

    pub fn interchangeable(&self, other: &Topology) -> bool {
        self.get_storage().atoms.len() == other.get_storage().atoms.len()
            && self.get_storage().bonds.len() == other.get_storage().bonds.len()
            && self.get_storage().molecules.len() == other.get_storage().molecules.len()
    }
}

//---------------------------
macro_rules! impl_topology_traits {
    ( $t:ty ) => {
        impl TopologyIoProvider for $t {}

        impl AtomIterProvider for $t {
            fn iter_atoms(&self) -> impl super::AtomIterator<'_> {
                self.get_storage().atoms.iter()
            }
        }

        impl AtomIterMutProvider for $t {
            fn iter_atoms_mut(&self) -> impl super::AtomMutIterator<'_> {
                self.get_storage_mut().atoms.iter_mut()
            }
        }

        impl MassIterProvider for $t {
            fn iter_masses(&self) -> impl Iterator<Item = f32> {
                self.get_storage().atoms.iter().map(|at| at.mass)
            }
        }

        impl LenProvider for $t {
            fn len(&self) -> usize {
                self.get_storage_mut().atoms.len()
            }
        }

        impl RandomAtomProvider for $t {
            unsafe fn get_atom_unchecked(&self, i: usize) -> &Atom {
                self.get_storage().atoms.get_unchecked(i)
            }
        }

        impl RandomAtomMutProvider for $t {
            fn get_atom_mut(&self, i: usize) -> Option<&mut Atom> {
                self.get_storage_mut().atoms.get_mut(i)
            }

            unsafe fn get_atom_mut_unchecked(&self, i: usize) -> &mut Atom {
                self.get_storage_mut().atoms.get_unchecked_mut(i)
            }
        }

        impl BondsProvider for $t {
            fn num_bonds(&self) -> usize {
                self.get_storage().bonds.len()
            }

            unsafe fn get_bond_unchecked(&self, i: usize) -> &[usize;2] {
                self.get_storage().bonds.get_unchecked(i)
            }

            fn iter_bonds(&self) -> impl Iterator<Item = &[usize; 2]> {
                self.get_storage().bonds.iter()
            }
        }
        
        impl MoleculesProvider for $t {
            fn num_molecules(&self) -> usize {
                self.get_storage().molecules.len()
            }
            
            unsafe fn get_molecule_unchecked(&self, i: usize) -> &[usize;2] {
                self.get_storage().molecules.get_unchecked(i)
            }

            fn iter_molecules(&self) -> impl Iterator<Item = &[usize; 2]> {
                self.get_storage().molecules.iter()
            }
        }
    };
}

// Impls for Topology itself
impl_topology_traits!(Topology);
// Impls for smart pointers
impl_topology_traits!(triomphe::Arc<Topology>);