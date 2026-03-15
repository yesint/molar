use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator};

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

/// Trait for parallel iteration over indexes
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
// LenProvider
//--------------------------------------------------------------

/// Trait for providing length of provided data
pub trait LenProvider {
    fn len(&self) -> usize;
}

//--------------------------------------------------------------
// New element traits (v5)
//--------------------------------------------------------------

/// Trait for random access to positions by absolute index.
/// `pos_unchecked(i)` takes the *absolute* storage index (not the selection-local index).
pub trait PosProvider: LenProvider {
    unsafe fn pos_unchecked(&self, i: usize) -> &Pos;
}

/// Trait for mutable random access to positions.
pub trait PosMutProvider: PosProvider {
    unsafe fn pos_mut_unchecked(&mut self, i: usize) -> &mut Pos;
}

/// Trait for random access to atoms by absolute index.
pub trait AtomProvider: LenProvider {
    unsafe fn atom_unchecked(&self, i: usize) -> &Atom;
}

/// Trait for mutable random access to atoms.
pub trait AtomMutProvider: AtomProvider {
    unsafe fn atom_mut_unchecked(&mut self, i: usize) -> &mut Atom;
}

/// Trait for access to periodic box
pub trait BoxProvider {
    fn get_box(&self) -> Option<&PeriodicBox>;

    fn require_box(&self) -> Result<&PeriodicBox, PeriodicBoxError> {
        self.get_box().ok_or(PeriodicBoxError::NoPbc)
    }
}

/// Trait for mutable access to periodic box
pub trait BoxMutProvider: BoxProvider {
    fn get_box_mut(&mut self) -> Option<&mut PeriodicBox>;
}

/// Trait for getting time
pub trait TimeProvider {
    fn get_time(&self) -> f32;
}

/// Trait for setting simulation time
pub trait TimeMutProvider: TimeProvider {
    fn set_time(&mut self, t: f32);
}

/// Trait for access to bonds
pub trait BondProvider {
    fn bonds_raw(&self) -> &[[usize; 2]];

    fn iter_bonds(&self) -> impl ExactSizeIterator<Item = &[usize; 2]> {
        self.bonds_raw().iter()
    }

    fn iter_bonds_dyn<'a>(&'a self) -> Box<dyn ExactSizeIterator<Item = &'a [usize; 2]> + 'a> {
        Box::new(self.iter_bonds())
    }

    fn num_bonds(&self) -> usize {
        self.bonds_raw().len()
    }

    // For backwards compatibility
    unsafe fn get_bond_unchecked(&self, i: usize) -> &[usize; 2] {
        self.bonds_raw().get_unchecked(i)
    }

    fn get_bond(&self, i: usize) -> Option<&[usize; 2]> {
        self.bonds_raw().get(i)
    }
}

/// Trait for access to molecules
pub trait MolProvider {
    fn molecules_raw(&self) -> &[[usize; 2]];

    fn iter_molecules(&self) -> impl ExactSizeIterator<Item = &[usize; 2]> {
        self.molecules_raw().iter()
    }

    fn num_molecules(&self) -> usize {
        self.molecules_raw().len()
    }

    // For backwards compatibility
    unsafe fn get_molecule_unchecked(&self, i: usize) -> &[usize; 2] {
        self.molecules_raw().get_unchecked(i)
    }

    fn get_molecule(&self, i: usize) -> Option<&[usize; 2]> {
        self.molecules_raw().get(i)
    }
}

//--------------------------------------------------------------
// Backward-compat aliases (used in distance_search.rs, sasa.rs etc.)
// These delegate to the new traits.
//--------------------------------------------------------------

// PosIterProvider - NOT a blanket from PosProvider to avoid conflicts.
// Types implementing SysProvider get it via PosProvider::iter_pos()
// and can implement PosIterProvider separately if needed.
pub trait PosIterProvider {
    fn iter_pos(&self) -> impl Iterator<Item = &Pos>;
}

// AtomIterProvider - NOT a blanket from AtomProvider.
pub trait AtomIterProvider {
    fn iter_atoms(&self) -> impl Iterator<Item = &Atom>;
}

// MassIterProvider - NOT a blanket.
pub trait MassIterProvider {
    fn iter_masses(&self) -> impl Iterator<Item = f32>;
}

// Default: AtomIterProvider implies MassIterProvider
impl<T: AtomIterProvider> MassIterProvider for T {
    fn iter_masses(&self) -> impl Iterator<Item = f32> {
        self.iter_atoms().map(|a| a.mass)
    }
}

// PosParIterProvider
pub trait PosParIterProvider {
    fn par_iter_pos(&self) -> impl IndexedParallelIterator<Item = &Pos>;
}

// AtomParIterProvider
pub trait AtomParIterProvider {
    fn par_iter_atoms(&self) -> impl IndexedParallelIterator<Item = &Atom>;
}

// Keep RandomPosProvider as alias
pub trait RandomPosProvider: LenProvider {
    unsafe fn get_pos_unchecked(&self, i: usize) -> &Pos;

    fn get_pos(&self, i: usize) -> Option<&Pos> {
        if i < self.len() { Some(unsafe { self.get_pos_unchecked(i) }) } else { None }
    }

    fn first_pos(&self) -> &Pos { unsafe { self.get_pos_unchecked(0) } }
    fn last_pos(&self) -> &Pos { unsafe { self.get_pos_unchecked(self.len()-1) } }
}

// Keep RandomAtomProvider as alias
pub trait RandomAtomProvider: LenProvider {
    unsafe fn get_atom_unchecked(&self, i: usize) -> &Atom;

    fn get_atom(&self, i: usize) -> Option<&Atom> {
        if i < self.len() { Some(unsafe { self.get_atom_unchecked(i) }) } else { None }
    }

    fn first_atom(&self) -> &Atom { unsafe { self.get_atom_unchecked(0) } }
    fn last_atom(&self) -> &Atom { unsafe { self.get_atom_unchecked(self.len()-1) } }
}

// Keep RandomPosMutProvider
pub trait RandomPosMutProvider: LenProvider {
    unsafe fn get_pos_mut_unchecked(&mut self, i: usize) -> &mut Pos;

    fn get_pos_mut(&mut self, i: usize) -> Option<&mut Pos> {
        if i < self.len() { Some(unsafe { self.get_pos_mut_unchecked(i) }) } else { None }
    }

    fn first_pos_mut(&mut self) -> &mut Pos { unsafe { self.get_pos_mut_unchecked(0) } }
    fn last_pos_mut(&mut self) -> &mut Pos {
        let last = self.len() - 1;
        unsafe { self.get_pos_mut_unchecked(last) }
    }
}

// Keep RandomAtomMutProvider
pub trait RandomAtomMutProvider: LenProvider {
    unsafe fn get_atom_mut_unchecked(&mut self, i: usize) -> &mut Atom;

    fn get_atom_mut(&mut self, i: usize) -> Option<&mut Atom> {
        if i < self.len() { Some(unsafe { self.get_atom_mut_unchecked(i) }) } else { None }
    }

    fn first_atom_mut(&mut self) -> &mut Atom { unsafe { self.get_atom_mut_unchecked(0) } }
}

// Blanket impls: PosProvider -> RandomPosProvider (via selection indexing)
impl<T: PosProvider + IndexProvider> RandomPosProvider for T {
    unsafe fn get_pos_unchecked(&self, i: usize) -> &Pos {
        let idx = self.get_index_unchecked(i);
        PosProvider::pos_unchecked(self, idx)
    }
}

// Blanket impls: PosMutProvider -> RandomPosMutProvider (via selection indexing)
impl<T: PosMutProvider + IndexProvider> RandomPosMutProvider for T {
    unsafe fn get_pos_mut_unchecked(&mut self, i: usize) -> &mut Pos {
        let idx = self.get_index_unchecked(i);
        PosMutProvider::pos_mut_unchecked(self, idx)
    }
}

// Blanket impls: AtomProvider -> RandomAtomProvider
impl<T: AtomProvider + IndexProvider> RandomAtomProvider for T {
    unsafe fn get_atom_unchecked(&self, i: usize) -> &Atom {
        let idx = self.get_index_unchecked(i);
        AtomProvider::atom_unchecked(self, idx)
    }
}

// Blanket impls: AtomMutProvider -> RandomAtomMutProvider
impl<T: AtomMutProvider + IndexProvider> RandomAtomMutProvider for T {
    unsafe fn get_atom_mut_unchecked(&mut self, i: usize) -> &mut Atom {
        let idx = self.get_index_unchecked(i);
        AtomMutProvider::atom_mut_unchecked(self, idx)
    }
}

// Keep PosIterMutProvider
pub trait PosIterMutProvider {
    fn iter_pos_mut(&mut self) -> impl Iterator<Item = &mut Pos>;
}

impl<T: PosMutProvider + IndexProvider> PosIterMutProvider for T {
    fn iter_pos_mut(&mut self) -> impl Iterator<Item = &mut Pos> {
        let len = self.len();
        let p = self as *mut T;
        (0..len).map(move |i| unsafe {
            let s = &mut *p;
            let idx = s.get_index_unchecked(i);
            s.pos_mut_unchecked(idx)
        })
    }
}

// Keep AtomIterMutProvider with the set_ methods
pub trait AtomIterMutProvider {
    fn iter_atoms_mut(&mut self) -> impl Iterator<Item = &mut Atom>;

    /// Sets same name to all selected atoms
    fn set_same_name(&mut self, val: &str)
    where
        Self: Sized,
    {
        use crate::atom::ATOM_NAME_EXPECT;
        let s = AtomStr::try_from_str(val).expect(ATOM_NAME_EXPECT);
        for a in self.iter_atoms_mut() {
            a.name = s;
        }
    }

    /// Sets same resname to all selected atoms
    fn set_same_resname(&mut self, val: &str)
    where
        Self: Sized,
    {
        use crate::atom::ATOM_RESNAME_EXPECT;
        let s = AtomStr::try_from_str(val).expect(ATOM_RESNAME_EXPECT);
        for a in self.iter_atoms_mut() {
            a.resname = s;
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

impl<T: AtomMutProvider + IndexProvider> AtomIterMutProvider for T {
    fn iter_atoms_mut(&mut self) -> impl Iterator<Item = &mut Atom> {
        let len = self.len();
        let p = self as *mut T;
        (0..len).map(move |i| unsafe {
            let s = &mut *p;
            let idx = s.get_index_unchecked(i);
            s.atom_mut_unchecked(idx)
        })
    }
}

// ParticleIterProvider
pub trait ParticleIterProvider {
    fn iter_particle(&self) -> impl Iterator<Item = Particle<'_>>;
}

// ParticleParIterProvider
pub trait ParticleParIterProvider {
    fn par_iter_particle(&self) -> impl IndexedParallelIterator<Item = Particle<'_>>;
}

// ParticleIterMutProvider
pub trait ParticleIterMutProvider {
    fn iter_particle_mut(&mut self) -> impl Iterator<Item = ParticleMut<'_>>;
}

// ParticleParIterMutProvider
pub trait ParticleParIterMutProvider {
    fn par_iter_particle_mut(&mut self) -> impl IndexedParallelIterator<Item = ParticleMut<'_>>;
}

// RandomParticleProvider
pub trait RandomParticleProvider: RandomPosProvider + RandomAtomProvider {
    unsafe fn get_particle_unchecked(&self, i: usize) -> Particle<'_>;

    fn first_particle(&self) -> Particle<'_> { unsafe { self.get_particle_unchecked(0) } }
    fn last_particle(&self) -> Particle<'_> {
        let last = self.len()-1;
        unsafe { self.get_particle_unchecked(last) }
    }
    fn get_particle(&self, i: usize) -> Option<Particle<'_>> {
        if i < self.len() { Some(unsafe { self.get_particle_unchecked(i) }) } else { None }
    }
}

// RandomParticleMutProvider
pub trait RandomParticleMutProvider: RandomPosMutProvider + RandomAtomMutProvider {
    unsafe fn get_particle_mut_unchecked(&mut self, i: usize) -> ParticleMut<'_>;

    fn first_particle_mut(&mut self) -> ParticleMut<'_> { unsafe { self.get_particle_mut_unchecked(0) } }
    fn last_particle_mut(&mut self) -> ParticleMut<'_> {
        let last = self.len()-1;
        unsafe { self.get_particle_mut_unchecked(last) }
    }
    fn get_particle_mut(&mut self, i: usize) -> Option<ParticleMut<'_>> {
        if i < self.len() { Some(unsafe { self.get_particle_mut_unchecked(i) }) } else { None }
    }
}

// Keep PosParIterMutProvider
pub trait PosParIterMutProvider {
    fn par_iter_pos_mut(&mut self) -> impl IndexedParallelIterator<Item = &mut Pos>;
}

// Keep AtomParIterProvider
pub trait AtomParIterMutProvider {
    fn par_iter_atoms_mut(&mut self) -> impl IndexedParallelIterator<Item = &mut Atom>;
}

// Keep RandomBondProvider for backward compat (maps to BondProvider)
pub trait RandomBondProvider {
    fn num_bonds(&self) -> usize;
    unsafe fn get_bond_unchecked(&self, i: usize) -> &[usize; 2];
    fn get_bond(&self, i: usize) -> Option<&[usize; 2]> {
        if i < self.num_bonds() { Some(unsafe { self.get_bond_unchecked(i) }) } else { None }
    }
}

// Blanket: BondProvider -> RandomBondProvider
impl<T: BondProvider> RandomBondProvider for T {
    fn num_bonds(&self) -> usize { self.bonds_raw().len() }
    unsafe fn get_bond_unchecked(&self, i: usize) -> &[usize; 2] {
        self.bonds_raw().get_unchecked(i)
    }
}

// Keep BondIterProvider for backward compat
pub trait BondIterProvider {
    fn iter_bonds(&self) -> impl Iterator<Item = &[usize; 2]>;
}

impl<T: BondProvider> BondIterProvider for T {
    fn iter_bonds(&self) -> impl Iterator<Item = &[usize; 2]> {
        BondProvider::iter_bonds(self)
    }
}

// Keep RandomMoleculeProvider for backward compat
pub trait RandomMoleculeProvider {
    fn num_molecules(&self) -> usize;
    unsafe fn get_molecule_unchecked(&self, i: usize) -> &[usize; 2];
    fn get_molecule(&self, i: usize) -> Option<&[usize; 2]> {
        if i < self.num_molecules() { Some(unsafe { self.get_molecule_unchecked(i) }) } else { None }
    }

    fn split_mol_iter(&self) -> impl Iterator<Item = Sel>
    where Self: Sized + IndexProvider,
    {
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
                        Some(0..=e - first)
                    } else if b >= first && e <= last {
                        Some(b - first..=e - first)
                    } else if b >= first && b <= last && e > last {
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

// Blanket: MolProvider -> RandomMoleculeProvider
impl<T: MolProvider> RandomMoleculeProvider for T {
    fn num_molecules(&self) -> usize { self.molecules_raw().len() }
    unsafe fn get_molecule_unchecked(&self, i: usize) -> &[usize; 2] {
        self.molecules_raw().get_unchecked(i)
    }
}

// Keep MoleculeIterProvider
pub trait MoleculeIterProvider {
    fn iter_molecules(&self) -> impl Iterator<Item = &[usize; 2]>;
}

impl<T: MolProvider> MoleculeIterProvider for T {
    fn iter_molecules(&self) -> impl Iterator<Item = &[usize; 2]> {
        MolProvider::iter_molecules(self)
    }
}

//--------------------------------------------------------------
// BoxMutProvider - also keep TimeMutProvider at module level
// (already defined as new trait above)
//--------------------------------------------------------------

//----------------------------------------------------
// Keep Vec<Pos> impls for State backward compat
impl PosIterProvider for Vec<Pos> {
    fn iter_pos(&self) -> impl Iterator<Item = &Pos> {
        self.iter()
    }
}

impl PosIterProvider for Pos {
    fn iter_pos(&self) -> impl Iterator<Item = &Pos> {
        std::iter::once(self)
    }
}
