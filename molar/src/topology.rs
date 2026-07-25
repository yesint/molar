use crate::prelude::*;
use thiserror::Error;

/// Topology of the molecular system: atoms, bonds, molecules, etc.
///
/// [Topology] is typically read from structure of trajectory file and is not intended
/// to be manipulated directly by the user. Insead [State](super::State) and [Topology]
/// are used to create atom selections, which give an access to the properties of
/// individual atoms and allow to query various properties.

#[derive(Debug, Default, Clone)]
pub struct Topology {
    pub atoms: AtomStorage,
    pub bonds: BondStorage,
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
        // The bonded adjacency is sized by the atom count, so a longer atom array leaves a
        // short (stale) index behind. Nothing here adds bonds, so only the count changed.
        self.bonds.invalidate_adjacency();
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

        // Drop bonds incident on a removed atom and renumber the survivors *before* the atom
        // columns shrink, since both take pre-removal indices.
        self.bonds.remove_by_index(&ind, self.atoms.len());
        self.atoms.remove_by_index(&ind);
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
    fn iter_bonds_dyn<'a>(&'a self) -> Box<dyn Iterator<Item = BondRef<'a>> + 'a> {
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

    unsafe fn get_bond_unchecked(&self, i: usize) -> BondRef<'_> { unsafe {
        self.bonds.get_unchecked(i)
    }}

    fn iter_bonds(&self) -> impl Iterator<Item = BondRef<'_>> {
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

#[cfg(test)]
mod tests {
    use crate::prelude::*;

    fn chain(n: usize) -> Topology {
        let mut t = Topology::default();
        for i in 0..n {
            t.atoms.push_row(&Atom::new().with_name("C").with_resid(i as i32));
        }
        for i in 0..n - 1 {
            t.bonds.push(&Bond::new(i, i + 1));
        }
        t
    }

    /// Before the columnar-bond migration `remove_atoms` shrank only the atom columns, leaving
    /// bonds pointing at dead or out-of-range indices.
    #[test]
    fn remove_atoms_drops_and_renumbers_bonds() {
        // 0-1-2-3-4; removing atom 2 must kill bonds 1-2 and 2-3.
        let mut t = chain(5);
        t.remove_atoms([2].into_iter()).unwrap();

        assert_eq!(t.atoms.len(), 4);
        assert_eq!(t.bonds.len(), 2, "both bonds incident on atom 2 are gone");
        for b in t.bonds.iter() {
            let [i, j] = b.pair();
            assert!(i < t.atoms.len() && j < t.atoms.len(), "bond {i}-{j} is out of range");
        }
        // The survivors are old 0-1 and old 3-4, renumbered onto the 4-atom space.
        let pairs: Vec<[usize; 2]> = t.bonds.iter().map(|b| b.pair()).collect();
        assert_eq!(pairs, vec![[0, 1], [2, 3]]);
        // ...and still connect the same atoms, identified by the resid they carried.
        let resid = |i: usize| t.atoms.get(i).unwrap().get_resid();
        assert_eq!((resid(2), resid(3)), (3, 4), "old atoms 3 and 4 moved down by one");
    }

    #[test]
    fn remove_atoms_invalidates_the_bond_adjacency() {
        let mut t = chain(5);
        t.bonds.ensure_adjacency(5);
        t.remove_atoms([2].into_iter()).unwrap();
        assert!(t.bonds.get_adjacency().is_none());
    }

    /// `add_atoms` grows the atom count without touching bonds, which would leave the
    /// adjacency's `offsets` short and panic on `neighbors(new_atom)`.
    #[test]
    fn add_atoms_invalidates_the_bond_adjacency() {
        let mut t = chain(3);
        t.bonds.ensure_adjacency(3);
        t.add_atoms(std::iter::once(Atom::new().with_name("X")));
        assert!(
            t.bonds.get_adjacency().is_none(),
            "a stale adjacency would be shorter than the atom count"
        );
        // Rebuilding covers the new atom, which is simply isolated.
        assert!(t.bonds.ensure_adjacency(t.atoms.len()).neighbors(3).is_empty());
    }
}
