use crate::prelude::*;
use thiserror::Error;

pub trait TopLike: Default + BondProvider + MolProvider + SaveTopology {
    // Construction methods (to use with format readers)
    fn add_atom_dyn(&mut self, atom: &dyn AtomLike);
    fn add_bond(&mut self, bond: &[usize;2]);
    fn add_molecule(&mut self, bond: &[usize;2]);
    // Retrieval methods (to use with format writers) are provided by base traits
}

/// Topology of the molecular system: atoms, bonds, molecules, etc.
///
/// [Topology] is typically read from structure of trajectory file and is not intended
/// to be manipulated directly by the user. Insead [State](super::State) and [Topology]
/// are used to create atom selections, which give an access to the properties of
/// individual atoms and allow to query various properties.

#[derive(Debug, Default, Clone)]
pub struct Topology {
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
    pub fn add_atoms<'a>(&'a mut self, atoms: impl Iterator<Item = Atom>) {
        self.atoms.extend(atoms);
    }

    pub fn remove_atoms(
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
impl SaveTopology for Topology {
    fn iter_atoms_dyn(&self) -> Box<dyn Iterator<Item = &Atom> + '_> {
        Box::new(self.iter_atoms())
    }
    fn iter_bonds_dyn<'a>(&'a self) -> Box<dyn Iterator<Item = &'a [usize; 2]> + 'a> {
        Box::new(BondProvider::iter_bonds(self))
    }
    fn num_bonds(&self) -> usize {
        BondProvider::num_bonds(self)
    }
}

impl LenProvider for Topology {
    fn len(&self) -> usize {
        self.atoms.len()
    }
}

/// Identity index provider for Topology (index i → atom i)
impl IndexProvider for Topology {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        i
    }

    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize> {
        0..self.atoms.len()
    }
}

impl AtomProvider for Topology {
    unsafe fn atoms_ptr(&self) -> *const Atom {
        self.atoms.as_ptr()
    }
}

impl AtomMutProvider for Topology {
    unsafe fn atoms_ptr_mut(&mut self) -> *mut Atom {
        self.atoms.as_mut_ptr()
    }
}

impl BondProvider for Topology {
    fn num_bonds(&self) -> usize {
        self.bonds.len()
    }

    unsafe fn get_bond_unchecked(&self, i: usize) -> &[usize; 2] {
        self.bonds.get_unchecked(i)
    }

    fn iter_bonds(&self) -> impl Iterator<Item = &[usize; 2]> {
        self.bonds.iter()
    }
}

impl MolProvider for Topology {
    fn num_molecules(&self) -> usize {
        self.molecules.len()
    }

    unsafe fn get_molecule_unchecked(&self, i: usize) -> &[usize; 2] {
        self.molecules.get_unchecked(i)
    }

    fn iter_molecules(&self) -> impl Iterator<Item = &[usize; 2]> {
        self.molecules.iter()
    }
}
