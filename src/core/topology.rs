use std::marker::PhantomData;

use sync_unsafe_cell::SyncUnsafeCell;
use thiserror::Error;
use triomphe::Arc;

use crate::io::TopologyProvider;

use super::{providers::{AtomsMutProvider, AtomsProvider, MassesProvider}, Atom, BuilderSerial, ImmutableParallel, MutableParallel, MutableSerial, SelectionError, SelectionKind};

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
pub struct Topology<K=()>{
    pub(crate) arc: Arc<SyncUnsafeCell<TopologyStorage>>,
    _marker: PhantomData<K>,
}

impl<K> Clone for Topology<K> {
    fn clone(&self) -> Self {
        Self{
            arc: Arc::new(SyncUnsafeCell::new(self.get_storage().clone())),
            _marker: Default::default(),
        }
    }
}

fn convert<K1,K2>(value: Topology<K1>) -> Result<Topology<K2>, SelectionError> {
    if Arc::count(&value.arc) == 1 {
        Ok(Topology::<K2> {
            arc: value.arc,
            _marker: Default::default(),
        })
    } else {
        Err(SelectionError::Release)
    }
}

impl TryFrom<Topology<()>> for Topology<MutableSerial> {
    type Error = SelectionError;
    fn try_from(value: Topology<()>) -> Result<Self, Self::Error> {
        convert(value)
    }
}

impl TryFrom<Topology<()>> for Topology<ImmutableParallel> {
    type Error = SelectionError;
    fn try_from(value: Topology<()>) -> Result<Self, Self::Error> {
        convert(value)
    }
}

impl TryFrom<Topology<()>> for Topology<MutableParallel> {
    type Error = SelectionError;
    fn try_from(value: Topology<()>) -> Result<Self, Self::Error> {
        convert(value)
    }
}

impl TryFrom<Topology<()>> for Topology<BuilderSerial> {
    type Error = SelectionError;
    fn try_from(value: Topology<()>) -> Result<Self, Self::Error> {
        convert(value)
    }
}

impl From<TopologyStorage> for Topology<()> {
    fn from(value: TopologyStorage) -> Self {
        Topology{
            arc: Arc::new(SyncUnsafeCell::new(value)),
            _marker: Default::default(),
        }
    }
}

impl<K> Topology<K> {
    pub fn new_arc<K2: SelectionKind>(&self) -> Topology<K2> {
        Topology::<K2>{
            arc: Arc::clone(&self.arc),
            _marker: Default::default(),
        }
    }

    // Private convenience accessors
    #[inline(always)]
    pub(super) fn get_storage(&self) -> &TopologyStorage {
        unsafe {&*self.arc.get()}
    }

    #[inline(always)]
    pub(super) fn get_storage_mut(&self) -> &mut TopologyStorage {
        unsafe {&mut *self.arc.get()}
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

    pub fn interchangeable(&self, other: &Topology<K>) -> bool {
            self.get_storage().atoms.len() == other.get_storage().atoms.len()
        &&  self.get_storage().bonds.len() == other.get_storage().bonds.len()
        &&  self.get_storage().molecules.len() == other.get_storage().molecules.len()
    }
}

// Impls for Topology itself
impl<K> TopologyProvider for Topology<K> {
    fn num_atoms(&self) -> usize {
        self.get_storage().atoms.len()
    }
}

impl<K> AtomsProvider for Topology<K> {
    fn iter_atoms(&self) -> impl super::AtomIterator<'_> {
        self.get_storage().atoms.iter()
    }
}

impl<K> AtomsMutProvider for Topology<K> {
    fn iter_atoms_mut(&self) -> impl super::AtomMutIterator<'_> {
        self.get_storage_mut().atoms.iter_mut()
    }
}

impl<K> MassesProvider for Topology<K> {
    fn iter_masses(&self) -> impl ExactSizeIterator<Item = f32> {
        self.get_storage().atoms.iter().map(|at| at.mass)
    }
}
