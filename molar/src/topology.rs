use crate::prelude::*;
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Chemical bond order. File formats that don't record it (PDB, GRO, XYZ, …) yield
/// [`BondOrder::Unspecified`]; formats that do (e.g. SDF, in the future) set the
/// concrete order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum BondOrder {
    /// The source didn't record an order (the common case for guessed/PDB bonds).
    #[default]
    Unspecified,
    Single,
    Double,
    Triple,
    Aromatic,
}

/// A bond between two atoms (by global index) with an optional chemical order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Bond {
    pub i1: usize,
    pub i2: usize,
    pub order: BondOrder,
}

impl Bond {
    /// A bond between `i1` and `i2` with an unspecified order.
    pub fn new(i1: usize, i2: usize) -> Self {
        Self { i1, i2, order: BondOrder::Unspecified }
    }

    /// A bond with an explicit order.
    pub fn with_order(i1: usize, i2: usize, order: BondOrder) -> Self {
        Self { i1, i2, order }
    }

    /// The two atom indices as a pair (for code that just needs the endpoints).
    pub fn pair(&self) -> [usize; 2] {
        [self.i1, self.i2]
    }

    /// Whether `idx` is one of this bond's endpoints.
    pub fn contains(&self, idx: usize) -> bool {
        self.i1 == idx || self.i2 == idx
    }
}

pub trait TopLike: Default + BondProvider + MolProvider + SaveTopology {
    // Construction methods (to use with format readers)
    fn add_atom_dyn(&mut self, atom: &dyn AtomLike);
    fn add_bond(&mut self, bond: Bond);
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
    pub atoms: AtomStorage,
    pub bonds: Vec<Bond>,
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

        self.atoms.retain_by_index(&ind);
        Ok(())
    }
}

impl Topology {
    pub fn assign_resindex(&mut self) {
        let mut resindex = 0usize;
        let mut cur_resid = unsafe { self.get_atom_unchecked(0) }.get_resid();
        for mut at in self.iter_atoms_mut() {
            if at.get_resid() != cur_resid {
                cur_resid = at.get_resid();
                resindex += 1;
            }
            at.set_resindex(resindex);
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
    fn iter_atoms_dyn(&self) -> Box<dyn Iterator<Item = AtomRef<'_>> + '_> {
        Box::new(self.iter_atoms())
    }
    fn iter_bonds_dyn<'a>(&'a self) -> Box<dyn Iterator<Item = &'a Bond> + 'a> {
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
    fn atom_storage(&self) -> &AtomStorage {
        &self.atoms
    }
}

impl AtomMutProvider for Topology {
    fn atom_storage_mut(&mut self) -> &mut AtomStorage {
        &mut self.atoms
    }
}

impl BondProvider for Topology {
    fn num_bonds(&self) -> usize {
        self.bonds.len()
    }

    unsafe fn get_bond_unchecked(&self, i: usize) -> &Bond { unsafe {
        self.bonds.get_unchecked(i)
    }}

    fn iter_bonds(&self) -> impl Iterator<Item = &Bond> {
        self.bonds.iter()
    }
}

impl MolProvider for Topology {
    fn num_molecules(&self) -> usize {
        self.molecules.len()
    }

    unsafe fn get_molecule_unchecked(&self, i: usize) -> &[usize; 2] { unsafe {
        self.molecules.get_unchecked(i)
    }}

    fn iter_molecules(&self) -> impl Iterator<Item = &[usize; 2]> {
        self.molecules.iter()
    }
}
