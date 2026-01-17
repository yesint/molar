use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator};

use crate::prelude::*;

//--------------------------------------------------------------
// Index
//--------------------------------------------------------------

/// Trait for selected indices
pub trait IndexProvider: LenProvider {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize;

    fn iter_index(&self) -> impl Iterator<Item = usize> {
        (0..self.len()).map(|i| unsafe { self.get_index_unchecked(i) })
    }
       
    fn first_index(&self) -> usize {
        unsafe { self.get_index_unchecked(0) }
    }

    fn last_index(&self) -> usize {
        unsafe { self.get_index_unchecked(self.len()-1) }
    }

    /// Creates a string in Gromacs index format representing self.
    fn as_gromacs_ndx_str(&self, name: impl AsRef<str>) -> String {
        use itertools::Itertools;
        let name = name.as_ref();
        let mut s = format!("[ {} ]\n", name);
        for chunk in &self.iter_index().chunks(15) {
            let line: String = chunk.map(|i| (i + 1).to_string()).join(" ");
            s.push_str(&line);
            s.push('\n');
        }
        s
    }
}

/// Trair for parallel iteration over indexes
pub trait IndexParProvider: LenProvider {
    fn par_iter_index(&self) -> impl IndexedParallelIterator<Item = usize>;
}

/// Trait for entities that contain continuous slices of indices (selections and such)
pub trait IndexSliceProvider {
    // Provides *sorted* indexes as slice
    fn get_index_slice(&self) -> &[usize];
}

impl<T: IndexSliceProvider> LenProvider for T {
    fn len(&self) -> usize {
        self.get_index_slice().len()
    }
}

impl<T: IndexSliceProvider> IndexProvider for T {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.get_index_slice().get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> {
        self.get_index_slice().into_iter().cloned()
    }
}

impl<T: IndexSliceProvider> IndexParProvider for T {
    fn par_iter_index(&self) -> impl IndexedParallelIterator<Item = usize> {
        self.get_index_slice().par_iter().cloned()
    }
}

impl IndexSliceProvider for SVec {
    fn get_index_slice(&self) -> &[usize] {
        self.as_slice()
    }
}

//--------------------------------------------------------------
// Immutable providers
//--------------------------------------------------------------

/// Trait for providing length of provided data
pub trait LenProvider {
    fn len(&self) -> usize;
}

/// Trait for providing iteration over positions
/// 
/// This trait is not dyn compatible, so it can't be combined
/// with dyn-compatible RandomPosProvider
pub trait PosIterProvider {
    fn iter_pos(&self) -> impl PosIterator<'_>;
}

/// Trait for providing parallel iteration over positions
pub trait PosParIterProvider {
    fn par_iter_pos(&self) -> impl IndexedParallelIterator<Item = &Pos>;
}

/// Trait for providing parallel mutable iteration over positions
pub trait PosParIterMutProvider {
    fn par_iter_pos_mut(&mut self) -> impl IndexedParallelIterator<Item = &mut Pos>;
}

/// Trait for providing iteration over atomic masses
pub trait MassIterProvider {
    fn iter_masses(&self) -> impl Iterator<Item = f32>;
}

/// Trait for providing iteration over atoms
pub trait AtomIterProvider {
    fn iter_atoms(&self) -> impl AtomIterator<'_>;
}

pub trait AtomParIterProvider {
    fn par_iter_atoms(&self) -> impl IndexedParallelIterator<Item = &Atom>;
}

/// Trait for providing parallel mutable iteration over atoms
pub trait AtomParIterMutProvider {
    fn par_iter_atoms_mut(&mut self) -> impl IndexedParallelIterator<Item = &mut Atom>;
}

impl<T: AtomIterProvider> MassIterProvider for T {
    fn iter_masses(&self) -> impl Iterator<Item = f32> {
        self.iter_atoms().map(|at| at.mass)
    }
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

/// Trait for getting time
pub trait TimeProvider {
    fn get_time(&self) -> f32;
}

/// Trait for providing iteration over particles
pub trait ParticleIterProvider {
    fn iter_particle(&self) -> impl Iterator<Item = Particle<'_>>;
}

pub trait ParticleParIterProvider {
    fn par_iter_particle(&self) -> impl IndexedParallelIterator<Item = Particle<'_>>;
}

/// Trait for providing parallel mutable iteration over particles
pub trait ParticleParIterMutProvider {
    fn par_iter_particle_mut(&mut self) -> impl IndexedParallelIterator<Item = ParticleMut<'_>>;
}

/// Trait for providing random access to particles
pub trait RandomParticleProvider: RandomPosProvider + RandomAtomProvider {
    unsafe fn get_particle_unchecked(&self, i: usize) -> Particle<'_>;

    fn first_particle(&self) -> Particle<'_> {
        unsafe { self.get_particle_unchecked(0) }
    }

    fn last_particle(&self) -> Particle<'_> {
        unsafe { self.get_particle_unchecked(self.len() - 1) }
    }

    fn get_particle(&self, i: usize) -> Option<Particle<'_>> {
        if i < self.len() {
            Some(unsafe { self.get_particle_unchecked(i) })
        } else {
            None
        }
    }
}

/// Trait for providing random access to positions
pub trait RandomPosProvider: LenProvider {
    unsafe fn get_pos_unchecked(&self, i: usize) -> &Pos;

    fn get_pos(&self, i: usize) -> Option<&Pos> {
        if i < self.len() {
            Some(unsafe { self.get_pos_unchecked(i) })
        } else {
            None
        }
    }

    fn first_pos(&self) -> &Pos {
        unsafe { self.get_pos_unchecked(0) }
    }

    fn last_pos(&self) -> &Pos {
        unsafe { self.get_pos_unchecked(self.len() - 1) }
    }
}

/// Trait for providing random access to atoms
pub trait RandomAtomProvider: LenProvider {
    unsafe fn get_atom_unchecked(&self, i: usize) -> &Atom;

    fn get_atom(&self, i: usize) -> Option<&Atom> {
        if i < self.len() {
            Some(unsafe { self.get_atom_unchecked(i) })
        } else {
            None
        }
    }

    fn first_atom(&self) -> &Atom {
        unsafe { self.get_atom_unchecked(0) }
    }

    fn last_atom(&self) -> &Atom {
        unsafe { self.get_atom_unchecked(self.len() - 1) }
    }
}

/// Trait for providing access to bonds
pub trait RandomBondProvider {
    fn num_bonds(&self) -> usize;

    unsafe fn get_bond_unchecked(&self, i: usize) -> &[usize; 2];

    fn get_bond(&self, i: usize) -> Option<&[usize; 2]> {
        if i < self.num_bonds() {
            Some(unsafe { self.get_bond_unchecked(i) })
        } else {
            None
        }
    }
}

pub trait BondIterProvider {
    fn iter_bonds(&self) -> impl Iterator<Item = &[usize; 2]>;
}

/// Trait for providing access to molecules (atoms subsets connected by bonds)
pub trait RandomMoleculeProvider {
    fn num_molecules(&self) -> usize;

    unsafe fn get_molecule_unchecked(&self, i: usize) -> &[usize; 2];

    fn get_molecule(&self, i: usize) -> Option<&[usize; 2]> {
        if i < self.num_molecules() {
            Some(unsafe { self.get_molecule_unchecked(i) })
        } else {
            None
        }
    }

    /// Splits by molecule and returns an iterator over them.
    /// If molecule is only partially contained in self then only this part is returned (molecules are clipped).
    /// If there are no molecules in [Topology] return an empty iterator.
    fn split_mol_iter(&self) -> impl Iterator<Item = Sel>
    where
        Self: Sized + IndexProvider,
    {
        // Iterate over molecules and find those inside selection
        let first = self.first_index();
        let last = self.last_index();

        let mut molid = 0;

        let next_fn = move || {
            if self.num_molecules() == 0 {
                return None;
            }

            let res = match self.get_molecule(molid) {
                Some(mol) => {
                    let b = mol[0];
                    let e = mol[1];
                    if b < first && e >= first && e <= last {
                        // molecule starts before Sel
                        Some(0..=e - first)
                    } else if b >= first && e <= last {
                        // molecule inside Sel
                        Some(b - first..=e - first)
                    } else if b >= first && b <= last && e > last {
                        // molecule ends after Sel
                        Some(b - first..=last - first)
                    } else {
                        None
                    }
                }
                None => None,
            }
            .map(|r| Sel(unsafe { SVec::from_sorted(r.into_iter().collect()) }));

            molid += 1;
            res
        };

        std::iter::from_fn(next_fn)
    }
}

pub trait MoleculeIterProvider {
    fn iter_molecules(&self) -> impl Iterator<Item = &[usize; 2]>;
}

//--------------------------------------------------------------
// Mutable providers
//--------------------------------------------------------------

/// Trait for providing mutable iteration over positions
pub trait PosIterMutProvider {
    fn iter_pos_mut(&mut self) -> impl PosMutIterator<'_>;
}

/// Trait for providing mutable random access to positions
pub trait RandomPosMutProvider: LenProvider {
    unsafe fn get_pos_mut_unchecked(&mut self, i: usize) -> &mut Pos;

    fn get_pos_mut(&mut self, i: usize) -> Option<&mut Pos> {
        if i < self.len() {
            Some(unsafe { self.get_pos_mut_unchecked(i) })
        } else {
            None
        }
    }

    fn first_pos_mut(&mut self) -> &mut Pos {
        unsafe { self.get_pos_mut_unchecked(0) }
    }

    fn last_pos_mut(&mut self) -> &mut Pos {
        unsafe { self.get_pos_mut_unchecked(self.len() - 1) }
    }
}

/// Trait for providing mutable iteration over atoms
pub trait AtomIterMutProvider {
    fn iter_atoms_mut(&mut self) -> impl AtomMutIterator<'_>;

    /// Sets same name to all selected atoms
    fn set_same_name(&mut self, val: &str)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.name = val.into();
        }
    }

    /// Sets same resname to all selected atoms
    fn set_same_resname(&mut self, val: &str)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.resname = val.into();
        }
    }

    /// Sets same resid to all selected atoms
    fn set_same_resid(&mut self, val: i32)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.resid = val;
        }
    }

    /// Sets same chain to all selected atoms
    fn set_same_chain(&mut self, val: char)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.chain = val;
        }
    }

    /// Sets same mass to all selected atoms
    fn set_same_mass(&mut self, val: f32)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.mass = val;
        }
    }

    /// Sets same B-factor to all selected atoms
    fn set_same_bfactor(&mut self, val: f32)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.bfactor = val;
        }
    }
}


/// Trait for providing mutable iteration over particles
pub trait ParticleIterMutProvider: IndexProvider {
    fn iter_particle_mut(&mut self) -> impl Iterator<Item = ParticleMut<'_>>;
}

/// Trait for providing mutable random access to atoms
pub trait RandomAtomMutProvider: LenProvider {
    fn get_atom_mut(&mut self, i: usize) -> Option<&mut Atom> {
        if i < self.len() {
            Some(unsafe { self.get_atom_mut_unchecked(i) })
        } else {
            None
        }
    }

    unsafe fn get_atom_mut_unchecked(&mut self, i: usize) -> &mut Atom;

    fn first_atom_mut(&mut self) -> &mut Atom {
        unsafe { self.get_atom_mut_unchecked(0) }
    }

    fn last_atom_mut(&mut self) -> &Atom {
        unsafe { self.get_atom_mut_unchecked(self.len() - 1) }
    }
}

/// Trait for providing mutable access to periodic box
pub trait BoxMutProvider {
    fn get_box_mut(&mut self) -> Option<&mut PeriodicBox>;
}

/// Trait for setting simulation time
pub trait TimeMutProvider {
    fn set_time(&mut self, t: f32);
}

/// Trait for providing mutable random access to particles
pub trait RandomParticleMutProvider: RandomPosMutProvider + RandomAtomMutProvider {
    unsafe fn get_particle_mut_unchecked(&mut self, i: usize) -> ParticleMut<'_>;

    fn first_particle_mut(&mut self) -> ParticleMut<'_> {
        unsafe { self.get_particle_mut_unchecked(0) }
    }

    fn last_particle_mut(&mut self) -> ParticleMut<'_> {
        unsafe { self.get_particle_mut_unchecked(self.len() - 1) }
    }

    fn get_particle_mut(&mut self, i: usize) -> Option<ParticleMut<'_>> {
        if i < self.len() {
            Some(unsafe { self.get_particle_mut_unchecked(i) })
        } else {
            None
        }
    }
}

//----------------------------------------------------
impl PosIterProvider for Vec<Pos> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.iter()
    }
}

impl PosIterProvider for Pos {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        std::iter::once(self)
    }
}
