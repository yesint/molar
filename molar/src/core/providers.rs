use crate::prelude::*;
use sorted_vec::SortedSet;

//--------------------------------------------------------------
// Basic providers
//--------------------------------------------------------------

/// Trait for providing both topology I/O operations
pub trait TopologyIoProvider: RandomAtomProvider + MoleculesProvider + BondsProvider {}

/// Trait for providing state I/O operations
pub trait StateIoProvider: RandomPosProvider + BoxProvider {
    fn get_time(&self) -> f32;
}

/// Trait for providing file writing for topology and state data
pub trait WritableToFile: TopologyIoProvider + StateIoProvider
where
    Self: Sized,
{
    fn save(&self, fname: &str) -> Result<(), FileIoError> {
        let mut h = FileHandler::create(fname)?;
        h.write(self)
    }
}

//--------------------------------------------------------------
// Index
//--------------------------------------------------------------

/// Trait for providing iteration over selected indices
pub trait IndexProvider {
    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize>;
}

impl IndexProvider for SortedSet<usize> {
    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize> {
        self.iter().cloned()
    }
}

impl IndexProvider for Vec<usize> {
    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize> {
        self.iter().cloned()
    }
}

//--------------------------------------------------------------
// Immutable providers
//--------------------------------------------------------------

/// Trait for providing length of data
pub trait LenProvider {
    fn len(&self) -> usize;
}

/// Trait for providing iteration over positions
pub trait PosIterProvider {
    fn iter_pos(&self) -> impl PosIterator<'_>;
}

/// Trait for providing iteration over atomic masses 
pub trait MassIterProvider {
    fn iter_masses(&self) -> impl ExactSizeIterator<Item = f32>;
}

/// Trait for providing iteration over atoms
pub trait AtomIterProvider {
    fn iter_atoms(&self) -> impl AtomIterator<'_>;
}

/// Trait for providing access to periodic box
pub trait BoxProvider {
    /// Get reference to the periodic box or `None` if there is no box.
    fn get_box(&self) -> Option<&PeriodicBox>;
    
    /// Get reference to the periodic box or an error if there is no box.
    fn require_box(&self) -> Result<&PeriodicBox, PeriodicBoxError> {
        self.get_box().ok_or_else(|| PeriodicBoxError::NoPbc)
    }
}

/// Trait for providing iteration over particles
pub trait ParticleIterProvider: IndexProvider {
    fn iter_particle(&self) -> impl ExactSizeIterator<Item = Particle<'_>>;
}

/// Trait for providing random access to particles
pub trait RandomParticleProvider: RandomPosProvider+RandomAtomProvider {
    unsafe fn nth_particle_unchecked(&self, i: usize) -> Particle<'_>;

    fn first_particle(&self) -> Particle {
        unsafe { self.nth_particle_unchecked(0) }
    }

    fn last_particle(&self) -> Particle {
        unsafe { self.nth_particle_unchecked(self.num_atoms() - 1) }
    }

    fn nth_particle(&self, i: usize) -> Option<Particle> {
        if i < self.num_pos() {
            Some(unsafe { self.nth_particle_unchecked(i) })
        } else {
            None
        }
    }
}

/// Trait for providing random access to positions
pub trait RandomPosProvider: PosIterProvider {
    fn num_pos(&self) -> usize;

    unsafe fn nth_pos_unchecked(&self, i: usize) -> &Pos;

    fn nth_pos(&self, i: usize) -> Option<&Pos> {
        if i < self.num_pos() {
            Some(unsafe { self.nth_pos_unchecked(i) })
        } else {
            None
        }
    }

    fn first_pos(&self) -> &Pos {
        unsafe { self.nth_pos_unchecked(0) }
    }

    fn last_pos(&self) -> &Pos {
        unsafe { self.nth_pos_unchecked(self.num_pos()-1) }
    }
}

/// Trait for providing random access to atoms
pub trait RandomAtomProvider: AtomIterProvider {
    fn num_atoms(&self) -> usize;

    unsafe fn nth_atom_unchecked(&self, i: usize) -> &Atom;

    fn nth_atom(&self, i: usize) -> Option<&Atom> {
        if i < self.num_atoms() {
            Some(unsafe { self.nth_atom_unchecked(i) })
        } else {
            None
        }
    }

    fn first_atom(&self) -> &Atom {
        unsafe { self.nth_atom_unchecked(0) }
    }

    fn last_atom(&self) -> &Atom {
        unsafe { self.nth_atom_unchecked(self.num_atoms() - 1) }
    }
}

/// Trait for providing access to bonds
pub trait BondsProvider {
    fn num_bonds(&self) -> usize;

    fn iter_bonds(&self) -> impl Iterator<Item = &[usize; 2]>;

    unsafe fn nth_bond_unchecked(&self, i: usize) -> &[usize; 2];

    fn nth_bond(&self, i: usize) -> Option<&[usize; 2]> {
        if i < self.num_bonds() {
            Some(unsafe { self.nth_bond_unchecked(i) })
        } else {
            None
        }
    }
}

/// Trait for providing access to molecules (atoms subsets connected by bonds)
pub trait MoleculesProvider {
    fn num_molecules(&self) -> usize;

    fn iter_molecules(&self) -> impl Iterator<Item = &[usize; 2]>;

    unsafe fn nth_molecule_unchecked(&self, i: usize) -> &[usize; 2];

    fn nth_molecule(&self, i: usize) -> Option<&[usize; 2]> {
        if i < self.num_molecules() {
            Some(unsafe { self.nth_molecule_unchecked(i) })
        } else {
            None
        }
    }
}

//--------------------------------------------------------------
// Mutable providers
//--------------------------------------------------------------

/// Trait for providing mutable iteration over positions
pub trait PosIterMutProvider: PosIterProvider {
    fn iter_pos_mut(&self) -> impl PosMutIterator<'_>;
}

/// Trait for providing mutable random access to positions
pub trait RandomPosMutProvider: RandomPosProvider + PosIterMutProvider {
    unsafe fn nth_pos_mut_unchecked(&self, i: usize) -> &mut Pos;

    fn nth_pos_mut(&self, i: usize) -> Option<&mut Pos> {
        if i < self.num_pos() {
            Some(unsafe { self.nth_pos_mut_unchecked(i) })
        } else {
            None
        }
    }
    
    fn first_pos_mut(&self) -> &mut Pos {
        unsafe { self.nth_pos_mut_unchecked(0) }
    }
    
    fn last_pos_mut(&self) -> &mut Pos {
        unsafe { self.nth_pos_mut_unchecked(self.num_pos()-1) }
    }
}

/// Trait for providing mutable iteration over atoms
pub trait AtomsIterMutProvider: AtomIterProvider {
    fn iter_atoms_mut(&self) -> impl AtomMutIterator<'_>;
}

/// Trait for providing mutable iteration over particles
pub trait ParticleIterMutProvider: IndexProvider {
    fn iter_particle_mut(&self) -> impl ExactSizeIterator<Item = ParticleMut<'_>>;
}

/// Trait for providing mutable random access to atoms
pub trait RandomAtomMutProvider: RandomAtomProvider {
    fn nth_atom_mut(&self, i: usize) -> Option<&mut Atom> {
        if i < self.num_atoms() {
            Some(unsafe { self.nth_atom_mut_unchecked(i) })
        } else {
            None
        }
    }

    unsafe fn nth_atom_mut_unchecked(&self, i: usize) -> &mut Atom;

    fn first_atom_mut(&self) -> &mut Atom {
        unsafe { self.nth_atom_mut_unchecked(0) }
    }

    fn last_atom_mut(&self) -> &Atom {
        unsafe { self.nth_atom_mut_unchecked(self.num_atoms() - 1) }
    }
}

/// Trait for providing mutable access to periodic box
pub trait BoxMutProvider {
    fn get_box_mut(&self) -> Option<&mut PeriodicBox>;
}

/// Trait for providing mutable random access to particles
pub trait RandomParticleMutProvider: RandomPosMutProvider+RandomAtomMutProvider {
    unsafe fn nth_particle_mut_unchecked(&self, i: usize) -> ParticleMut;

    fn first_particle_mut(&self) -> ParticleMut {
        unsafe { self.nth_particle_mut_unchecked(0) }
    }

    fn last_particle_mut(&self) -> ParticleMut {
        unsafe { self.nth_particle_mut_unchecked(self.num_atoms() - 1) }
    }

    fn nth_particle_mut(&self, i: usize) -> Option<ParticleMut> {
        if i < self.num_pos() {
            Some(unsafe { self.nth_particle_mut_unchecked(i) })
        } else {
            None
        }
    }
}

//----------------------------------------------------
impl PosIterProvider for Vec<Pos> {
    fn iter_pos(&self) -> impl PosIterator<'_> + Clone {
        self.iter()
    }
}

impl PosIterProvider for Pos {
    fn iter_pos(&self) -> impl PosIterator<'_> + Clone {
        std::iter::once(self)
    }
}
