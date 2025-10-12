use crate::prelude::*;
use thiserror::Error;

#[doc(hidden)]
#[derive(Debug, Default, Clone)]
pub(crate) struct Topology {
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

impl Topology {
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

impl Topology {
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

    pub fn interchangeable(&self, other: &Topology) -> bool {
        self.atoms.len() == other.atoms.len()
            && self.bonds.len() == other.bonds.len()
            && self.molecules.len() == other.molecules.len()
    }
}

//---------------------------
impl TopologyWrite for Topology {}

impl AtomIterProvider for Topology {
    fn iter_atoms(&self) -> impl super::AtomIterator<'_> {
        self.atoms.iter()
    }
}

impl AtomIterMutProvider for Topology {
    fn iter_atoms_mut(&mut self) -> impl super::AtomMutIterator<'_> {
        self.atoms.iter_mut()
    }
}

impl LenProvider for Topology {
    fn len(&self) -> usize {
        self.atoms.len()
    }
}

impl RandomAtomProvider for Topology {
    unsafe fn get_atom_unchecked(&self, i: usize) -> &Atom {
        self.atoms.get_unchecked(i)
    }
}

impl RandomAtomMutProvider for Topology {
    fn get_atom_mut(&mut self, i: usize) -> Option<&mut Atom> {
        self.atoms.get_mut(i)
    }

    unsafe fn get_atom_mut_unchecked(&mut self, i: usize) -> &mut Atom {
        self.atoms.get_unchecked_mut(i)
    }
}

impl RandomBondProvider for Topology {
    fn num_bonds(&self) -> usize {
        self.bonds.len()
    }

    unsafe fn get_bond_unchecked(&self, i: usize) -> &[usize; 2] {
        self.bonds.get_unchecked(i)
    }
}

impl BondIterProvider for Topology {
    fn iter_bonds(&self) -> impl Iterator<Item = &[usize; 2]> {
        self.bonds.iter()
    }
}

impl RandomMoleculeProvider for Topology {
    fn num_molecules(&self) -> usize {
        self.molecules.len()
    }

    unsafe fn get_molecule_unchecked(&self, i: usize) -> &[usize; 2] {
        self.molecules.get_unchecked(i)
    }
}

impl MoleculeIterProvider for Topology {
    fn iter_molecules(&self) -> impl Iterator<Item = &[usize; 2]> {
        self.molecules.iter()
    }
}
