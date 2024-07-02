use std::ops::Deref;
use anyhow::bail;
use sync_unsafe_cell::SyncUnsafeCell;

use crate::io::TopologyProvider;

use super::{providers::{AtomsMutProvider, AtomsProvider, MassesProvider}, Atom, TopologyUArc};

#[doc(hidden)]
#[derive(Debug, Default, Clone)]
pub(crate) struct TopologyStorage {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<[usize; 2]>,
    pub molecules: Vec<[usize; 2]>,
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
        Self(SyncUnsafeCell::new(self.get().clone()))
    }
}

impl From<TopologyStorage> for TopologyUArc {
    fn from(value: TopologyStorage) -> Self {
        Self::new(Topology(SyncUnsafeCell::new(value)))
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

    pub fn interchangeable(&self, other: &Topology) -> bool {
            self.get().atoms.len() == other.get().atoms.len()
        &&  self.get().bonds.len() == other.get().bonds.len()
        &&  self.get().molecules.len() == other.get().molecules.len()
    }

    pub fn add_atoms(&mut self, added: impl Iterator<Item = Atom>) {
        self.get_mut().atoms.extend(added);
    }

    pub fn add_bonds(&mut self, added: impl Iterator<Item = [usize; 2]>) {
        self.get_mut().bonds.extend(added);
    }

    pub fn add_molecules(&mut self, added: impl Iterator<Item = [usize; 2]>) {
        self.get_mut().molecules.extend(added);
    }

    pub fn remove(&mut self, removed: impl Iterator<Item = usize>) -> anyhow::Result<()> {
        let mut ind = removed.collect::<Vec<_>>();
        if ind.len()==0 {
            return Ok(());
        }
        ind.sort_unstable();
        ind.dedup();
        if ind[0] >= self.num_atoms() || ind[ind.len()-1] >= self.num_atoms() {
            bail!(
                "Indexes to remove [{}:{}] are out of allowed range [0:{}]",
                ind[0],ind[ind.len()-1],self.num_atoms()
            );
        }
        for i in ind.iter().rev().cloned() {
            self.get_mut().atoms.remove(i);
            // Remove affected bonds
            self.get_mut().bonds.retain(|b| b[0] != i && b[1] != i);
            // Modify molecules and remove those, which become invalid
            self.get_mut().molecules.retain_mut(|m| {
                if m[0]==i {m[0] += 1}
                if m[1]==i {
                    if m[1]>0 {
                        m[0] -= 1
                    } else {
                        return false;
                    }
                }
                m[1]>=m[0]
            });
        }
        Ok(())
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