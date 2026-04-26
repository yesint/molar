use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator};

use crate::atom::{ATOM_NAME_EXPECT, ATOM_RESNAME_EXPECT};
use crate::prelude::*;

//--------------------------------------------------------------
// Index
//--------------------------------------------------------------

/// Trait for selected indices
pub trait IndexProvider: LenProvider {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize;

    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize> {
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

    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize> {
        self.get_index_slice().iter().cloned()
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
// Core trait
//--------------------------------------------------------------

/// Trait for providing length of provided data
pub trait LenProvider {
    fn len(&self) -> usize;
}

//==============================================================
// Element traits (v5)
//==============================================================

/// Element trait providing immutable access to positions via index.
///
/// Implementors must supply a raw pointer to the contiguous positions array.
/// All other methods (iteration, random access, parallel iteration) are provided
/// as defaults using `coords_ptr()` together with `IndexProvider`.
pub trait PosProvider: LenProvider + IndexProvider {
    /// Raw pointer to the beginning of the positions array.
    ///
    /// # Safety
    /// Pointer must remain valid as long as `self` is alive.
    unsafe fn coords_ptr(&self) -> *const Pos;

    fn iter_pos(&self) -> impl PosIterator<'_> {
        let cp = unsafe { self.coords_ptr() };
        unsafe { self.iter_index().map(move |i| &*cp.add(i)) }
    }

    fn par_iter_pos(&self) -> impl IndexedParallelIterator<Item = &Pos>
    where
        Self: IndexParProvider,
    {
        let p = unsafe { self.coords_ptr() } as usize; // trick to make pointer Sync
        unsafe { self.par_iter_index().map(move |i| &*(p as *const Pos).add(i)) }
    }

    unsafe fn get_pos_unchecked(&self, i: usize) -> &Pos {
        let ind = self.get_index_unchecked(i);
        &*self.coords_ptr().add(ind)
    }

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

/// Element trait providing mutable access to positions.
///
/// Extends `PosProvider` with mutable iteration and random access.
pub trait PosMutProvider: PosProvider {
    unsafe fn coords_ptr_mut(&mut self) -> *mut Pos {
        self.coords_ptr() as *mut Pos
    }

    fn iter_pos_mut(&mut self) -> impl PosMutIterator<'_> {
        (0..self.len()).map(|i| {
            let ind = unsafe { self.get_index_unchecked(i) };
            unsafe { &mut *self.coords_ptr_mut().add(ind) }
        })
    }

    fn par_iter_pos_mut(&mut self) -> impl IndexedParallelIterator<Item = &mut Pos>
    where
        Self: IndexParProvider,
    {
        let p = unsafe { self.coords_ptr_mut() } as usize; // trick to make pointer Sync
        unsafe { self.par_iter_index().map(move |i| &mut *(p as *mut Pos).add(i)) }
    }

    unsafe fn get_pos_mut_unchecked(&mut self, i: usize) -> &mut Pos {
        let ind = self.get_index_unchecked(i);
        &mut *self.coords_ptr_mut().add(ind)
    }

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

/// Element trait providing immutable access to atoms via index.
///
/// Implementors must supply a raw pointer to the contiguous atoms array.
/// All other methods are provided as defaults.
pub trait AtomProvider: LenProvider + IndexProvider {
    /// Raw pointer to the beginning of the atoms array.
    ///
    /// # Safety
    /// Pointer must remain valid as long as `self` is alive.
    unsafe fn atoms_ptr(&self) -> *const Atom;

    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        let ap = unsafe { self.atoms_ptr() };
        unsafe { self.iter_index().map(move |i| &*ap.add(i)) }
    }

    fn par_iter_atoms(&self) -> impl IndexedParallelIterator<Item = &Atom>
    where
        Self: IndexParProvider,
    {
        let p = unsafe { self.atoms_ptr() } as usize; // trick to make pointer Sync
        unsafe { self.par_iter_index().map(move |i| &*(p as *const Atom).add(i)) }
    }

    fn iter_masses(&self) -> impl Iterator<Item = Float> {
        self.iter_atoms().map(|at| at.mass)
    }

    unsafe fn get_atom_unchecked(&self, i: usize) -> &Atom {
        let ind = self.get_index_unchecked(i);
        &*self.atoms_ptr().add(ind)
    }

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

/// Element trait providing mutable access to atoms.
///
/// Extends `AtomProvider` with mutable iteration and random access,
/// plus convenience setters.
pub trait AtomMutProvider: AtomProvider {
    unsafe fn atoms_ptr_mut(&mut self) -> *mut Atom {
        self.atoms_ptr() as *mut Atom
    }

    fn iter_atoms_mut(&mut self) -> impl AtomMutIterator<'_> {
        (0..self.len()).map(|i| {
            let ind = unsafe { self.get_index_unchecked(i) };
            unsafe { &mut *self.atoms_ptr_mut().add(ind) }
        })
    }

    fn par_iter_atoms_mut(&mut self) -> impl IndexedParallelIterator<Item = &mut Atom>
    where
        Self: IndexParProvider,
    {
        let p = unsafe { self.atoms_ptr_mut() } as usize; // trick to make pointer Sync
        unsafe { self.par_iter_index().map(move |i| &mut *(p as *mut Atom).add(i)) }
    }

    unsafe fn get_atom_mut_unchecked(&mut self, i: usize) -> &mut Atom {
        let ind = self.get_index_unchecked(i);
        &mut *self.atoms_ptr_mut().add(ind)
    }

    fn get_atom_mut(&mut self, i: usize) -> Option<&mut Atom> {
        if i < self.len() {
            Some(unsafe { self.get_atom_mut_unchecked(i) })
        } else {
            None
        }
    }

    fn first_atom_mut(&mut self) -> &mut Atom {
        unsafe { self.get_atom_mut_unchecked(0) }
    }

    fn last_atom_mut(&mut self) -> &mut Atom {
        unsafe { self.get_atom_mut_unchecked(self.len() - 1) }
    }

    /// Sets same name to all selected atoms
    fn set_same_name(&mut self, val: &str)
    where
        Self: Sized,
    {
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
    fn set_same_mass(&mut self, val: Float)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.mass = val;
        }
    }

    /// Sets same B-factor to all selected atoms
    fn set_same_bfactor(&mut self, val: Float)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.bfactor = val;
        }
    }
}

/// Element trait providing immutable access to bonds.
///
/// Default methods provide `iter_bonds`, `num_bonds`, `get_bond`.
pub trait BondProvider {
    fn num_bonds(&self) -> usize;

    unsafe fn get_bond_unchecked(&self, i: usize) -> &[usize; 2];

    fn get_bond(&self, i: usize) -> Option<&[usize; 2]> {
        if i < self.num_bonds() {
            Some(unsafe { self.get_bond_unchecked(i) })
        } else {
            None
        }
    }

    fn iter_bonds(&self) -> impl Iterator<Item = &[usize; 2]>;
}

/// Element trait providing immutable access to molecules.
///
/// Default methods provide `iter_molecules`, `num_molecules`, `get_molecule`,
/// and `split_mol_iter`.
pub trait MolProvider {
    fn num_molecules(&self) -> usize;

    unsafe fn get_molecule_unchecked(&self, i: usize) -> &[usize; 2];

    fn get_molecule(&self, i: usize) -> Option<&[usize; 2]> {
        if i < self.num_molecules() {
            Some(unsafe { self.get_molecule_unchecked(i) })
        } else {
            None
        }
    }

    fn iter_molecules(&self) -> impl Iterator<Item = &[usize; 2]>;

    /// Splits by molecule and returns an iterator over them.
    /// If molecule is only partially contained in self then only this part is returned (molecules are clipped).
    /// If there are no molecules in [Topology] return an empty iterator.
    fn split_mol_iter(&self) -> impl Iterator<Item = Sel>
    where
        Self: Sized + IndexProvider,
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

//--------------------------------------------------------------
// Box and time providers (unchanged)
//--------------------------------------------------------------

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
    fn get_time(&self) -> Float;
}

/// Trait for providing mutable access to periodic box
pub trait BoxMutProvider {
    fn get_box_mut(&mut self) -> Option<&mut PeriodicBox>;
}

/// Trait for setting simulation time
pub trait TimeMutProvider {
    fn set_time(&mut self, t: Float);
}

//--------------------------------------------------------------
// Vel / Force providers
//--------------------------------------------------------------

/// Element trait providing immutable access to velocities via index.
///
/// Implementors supply one `unsafe` raw pointer method. Returns a null pointer
/// when velocities are absent. All iteration and access methods are provided
/// as defaults using `vel_ptr()` and `IndexProvider`.
pub trait VelProvider: LenProvider + IndexProvider {
    /// Raw pointer to the beginning of the velocities array, or **null** if absent.
    ///
    /// # Safety
    /// If non-null, the pointer must remain valid as long as `self` is alive and
    /// must point to a contiguous array of at least `max_global_index + 1` elements.
    unsafe fn vel_ptr(&self) -> *const Vel;

    /// Returns `true` when velocities are present in the current state.
    fn has_vel(&self) -> bool {
        !unsafe { self.vel_ptr() }.is_null()
    }

    /// Iterates over velocities of selected atoms, or returns [`StateError::NoVelocities`] if absent.
    fn iter_vel(&self) -> Result<impl VelIterator<'_>, StateError> {
        let vp = unsafe { self.vel_ptr() };
        if vp.is_null() { return Err(StateError::NoVelocities); }
        Ok(unsafe { self.iter_index().map(move |i| &*vp.add(i)) })
    }

    /// Parallel iterator over velocities. Requires [`IndexParProvider`].
    fn par_iter_vel(&self) -> Result<impl IndexedParallelIterator<Item = &Vel>, StateError>
    where
        Self: IndexParProvider,
    {
        let vp = unsafe { self.vel_ptr() };
        if vp.is_null() { return Err(StateError::NoVelocities); }
        let p = vp as usize; // trick to make pointer Sync
        Ok(unsafe { self.par_iter_index().map(move |i| &*(p as *const Vel).add(i)) })
    }

    /// Random access to a velocity by local (selection) index. Returns `None` if out-of-bounds or absent.
    fn get_vel(&self, i: usize) -> Option<&Vel> {
        if i >= self.len() { return None; }
        let vp = unsafe { self.vel_ptr() };
        if vp.is_null() { return None; }
        Some(unsafe { &*vp.add(self.get_index_unchecked(i)) })
    }
}

/// Element trait providing mutable access to velocities.
///
/// Extends [`VelProvider`] with mutable iteration and random access.
pub trait VelMutProvider: VelProvider {
    /// Raw mutable pointer to the velocities array, or **null** if absent.
    ///
    /// # Safety
    /// Same contract as [`VelProvider::vel_ptr`].
    unsafe fn vel_ptr_mut(&mut self) -> *mut Vel {
        self.vel_ptr() as *mut Vel
    }

    /// Mutable iterator over velocities of selected atoms.
    fn iter_vel_mut(&mut self) -> Result<impl VelMutIterator<'_>, StateError> {
        let vp = unsafe { self.vel_ptr_mut() };
        if vp.is_null() { return Err(StateError::NoVelocities); }
        Ok((0..self.len()).map(move |i| {
            let ind = unsafe { self.get_index_unchecked(i) };
            unsafe { &mut *vp.add(ind) }
        }))
    }

    /// Parallel mutable iterator over velocities.
    fn par_iter_vel_mut(&mut self) -> Result<impl IndexedParallelIterator<Item = &mut Vel>, StateError>
    where
        Self: IndexParProvider,
    {
        let vp = unsafe { self.vel_ptr_mut() };
        if vp.is_null() { return Err(StateError::NoVelocities); }
        let p = vp as usize; // trick to make pointer Sync
        Ok(unsafe { self.par_iter_index().map(move |i| &mut *(p as *mut Vel).add(i)) })
    }

    /// Mutable random access to a velocity by local (selection) index.
    fn get_vel_mut(&mut self, i: usize) -> Option<&mut Vel> {
        if i >= self.len() { return None; }
        let vp = unsafe { self.vel_ptr_mut() };
        if vp.is_null() { return None; }
        Some(unsafe { &mut *vp.add(self.get_index_unchecked(i)) })
    }
}

/// Element trait providing immutable access to forces via index.
///
/// Mirror of [`VelProvider`] for forces.
pub trait ForceProvider: LenProvider + IndexProvider {
    /// Raw pointer to the beginning of the forces array, or **null** if absent.
    ///
    /// # Safety
    /// If non-null, the pointer must remain valid as long as `self` is alive.
    unsafe fn force_ptr(&self) -> *const Force;

    /// Returns `true` when forces are present in the current state.
    fn has_force(&self) -> bool {
        !unsafe { self.force_ptr() }.is_null()
    }

    /// Iterates over forces of selected atoms, or returns [`StateError::NoForces`] if absent.
    fn iter_force(&self) -> Result<impl ForceIterator<'_>, StateError> {
        let fp = unsafe { self.force_ptr() };
        if fp.is_null() { return Err(StateError::NoForces); }
        Ok(unsafe { self.iter_index().map(move |i| &*fp.add(i)) })
    }

    /// Parallel iterator over forces.
    fn par_iter_force(&self) -> Result<impl IndexedParallelIterator<Item = &Force>, StateError>
    where
        Self: IndexParProvider,
    {
        let fp = unsafe { self.force_ptr() };
        if fp.is_null() { return Err(StateError::NoForces); }
        let p = fp as usize; // trick to make pointer Sync
        Ok(unsafe { self.par_iter_index().map(move |i| &*(p as *const Force).add(i)) })
    }

    /// Random access to a force by local (selection) index.
    fn get_force(&self, i: usize) -> Option<&Force> {
        if i >= self.len() { return None; }
        let fp = unsafe { self.force_ptr() };
        if fp.is_null() { return None; }
        Some(unsafe { &*fp.add(self.get_index_unchecked(i)) })
    }
}

/// Element trait providing mutable access to forces.
///
/// Mirror of [`VelMutProvider`] for forces.
pub trait ForceMutProvider: ForceProvider {
    /// Raw mutable pointer to the forces array, or **null** if absent.
    ///
    /// # Safety
    /// Same contract as [`ForceProvider::force_ptr`].
    unsafe fn force_ptr_mut(&mut self) -> *mut Force {
        self.force_ptr() as *mut Force
    }

    /// Mutable iterator over forces of selected atoms.
    fn iter_force_mut(&mut self) -> Result<impl ForceMutIterator<'_>, StateError> {
        let fp = unsafe { self.force_ptr_mut() };
        if fp.is_null() { return Err(StateError::NoForces); }
        Ok((0..self.len()).map(move |i| {
            let ind = unsafe { self.get_index_unchecked(i) };
            unsafe { &mut *fp.add(ind) }
        }))
    }

    /// Parallel mutable iterator over forces.
    fn par_iter_force_mut(&mut self) -> Result<impl IndexedParallelIterator<Item = &mut Force>, StateError>
    where
        Self: IndexParProvider,
    {
        let fp = unsafe { self.force_ptr_mut() };
        if fp.is_null() { return Err(StateError::NoForces); }
        let p = fp as usize; // trick to make pointer Sync
        Ok(unsafe { self.par_iter_index().map(move |i| &mut *(p as *mut Force).add(i)) })
    }

    /// Mutable random access to a force by local (selection) index.
    fn get_force_mut(&mut self, i: usize) -> Option<&mut Force> {
        if i >= self.len() { return None; }
        let fp = unsafe { self.force_ptr_mut() };
        if fp.is_null() { return None; }
        Some(unsafe { &mut *fp.add(self.get_index_unchecked(i)) })
    }
}

//--------------------------------------------------------------
// Particle iterator (for iter_particle, par_iter_particle, etc.)
// These live here as concrete structs/methods are in Analysis trait
//--------------------------------------------------------------

/// Provides iteration over (id, atom, pos) particles
pub trait ParticleIterProvider {
    fn iter_particle(&self) -> impl Iterator<Item = Particle<'_>>;
}

pub trait ParticleParIterProvider {
    fn par_iter_particle(&self) -> impl IndexedParallelIterator<Item = Particle<'_>>;
}

/// Trait for providing mutable iteration over particles
pub trait ParticleIterMutProvider: IndexProvider {
    fn iter_particle_mut(&mut self) -> impl Iterator<Item = ParticleMut<'_>>;
}

//--------------------------------------------------------------
// Blanket impls for Vec<Pos> and Pos
// These objects don't use index indirection — they provide
// iter_pos() directly.
//--------------------------------------------------------------

/// Identity IndexProvider for Vec<Pos>
impl LenProvider for Vec<Pos> {
    fn len(&self) -> usize {
        Vec::len(self)
    }
}

impl IndexProvider for Vec<Pos> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        i
    }

    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize> {
        0..self.len()
    }
}

impl PosProvider for Vec<Pos> {
    unsafe fn coords_ptr(&self) -> *const Pos {
        self.as_ptr()
    }
}

impl LenProvider for Pos {
    fn len(&self) -> usize {
        1
    }
}

impl IndexProvider for Pos {
    unsafe fn get_index_unchecked(&self, _i: usize) -> usize {
        0
    }
    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize> {
        std::iter::once(0)
    }
}

impl PosProvider for Pos {
    unsafe fn coords_ptr(&self) -> *const Pos {
        self as *const Pos
    }
}
