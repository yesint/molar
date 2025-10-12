use std::path::Path;

use crate::prelude::*;

//================================================
/// Selection is just an index detached from the State.
/// Has to be bound to System before doing something with the subset of atoms.
/// It is guaranteed to be non-empty so no run-time checks for this
/// are needed when binding to system
//================================================
pub struct Sel(SVec);

impl Sel {
    pub fn from_vec(index: Vec<usize>) -> Result<Self, SelectionError> {
        if index.is_empty() {
            Err(SelectionError::EmptySlice)
        } else {
            Ok(Self(SVec::from_unsorted(index)))
        }
    }

    /// Binds Sel to System for read-only access
    pub fn bind<'a>(&'a self, sys: &'a System) -> Result<SubSystem<'a>, SelectionError> {
        let last_ind = self.0.len() - 1;
        if unsafe { *self.0.get_unchecked(last_ind) } < sys.top.len() {
            Ok(SubSystem {
                top: &sys.top,
                st: &sys.st,
                index: &self.0,
            })
        } else {
            Err(SelectionError::OutOfBounds(0, 0))
        }
    }

    /// Binds Sel to System for read-write access
    pub fn bind_mut<'a>(&'a self, sys: &'a mut System) -> Result<SubSystemMut<'a>, SelectionError> {
        // The cost calling is just one comparison
        let last_ind = self.0.len() - 1;
        if unsafe { *self.0.get_unchecked(last_ind) } < sys.top.len() {
            Ok(SubSystemMut {
                top: &mut sys.top,
                st: &mut sys.st,
                index: &self.0,
            })
        } else {
            Err(SelectionError::OutOfBounds(0, 0))
        }
    }
}

//================================================
/// System that stores Topology and State
//================================================
#[derive(Debug, Default)]
pub struct System {
    top: Topology,
    st: State,
}

impl System {
    pub fn new(top: Topology, st: State) -> Self {
        Self { top, st }
    }

    pub fn from_file(fname: impl AsRef<Path>) -> Result<Self, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        Ok(Self::new(top, st))
    }

    /// Create new selection based on provided definition.
    fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        Ok(Sel(def.into_sel_index(&self.top, &self.st, None)?))
    }
}

//================================================
/// Read only subsystem
/// Implements only read-only analysis traits
//================================================
pub struct SubSystem<'a> {
    top: &'a Topology,
    st: &'a State,
    index: &'a SVec,
}

impl SubSystem<'_> {
    /// Create new sub-selection based on provided definition.
    fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        Ok(Sel(def.into_sel_index(self.top, self.st, Some(self.index))?))
    }
}

impl LenProvider for SubSystem<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SubSystem<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl AtomPosAnalysis for SubSystem<'_> {
    fn atom_ptr(&self) -> *const Atom {
        self.top.atoms.as_ptr()
    }

    fn pos_ptr(&self) -> *const Pos {
        self.st.coords.as_ptr()
    }
}

impl NonAtomPosAnalysis for SubSystem<'_> {
    fn top_ref(&self) -> &Topology {
        &self.top
    }

    fn st_ref(&self) -> &State {
        &self.st
    }
}

//================================================
/// Read-write subsystem having access to all fields of Topology and State
//================================================
pub struct SubSystemMut<'a> {
    top: &'a mut Topology,
    st: &'a mut State,
    index: &'a SVec,
}

impl SubSystemMut<'_> {
    /// Create new sub-selection based on provided definition.
    fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        Ok(Sel(def.into_sel_index(self.top, self.st, Some(self.index))?))
    }
}

impl LenProvider for SubSystemMut<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SubSystemMut<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl AtomPosAnalysis for SubSystemMut<'_> {
    fn atom_ptr(&self) -> *const Atom {
        self.top.atoms.as_ptr()
    }

    fn pos_ptr(&self) -> *const Pos {
        self.st.coords.as_ptr()
    }
}


impl AtomPosAnalysisMut for SubSystemMut<'_> {
    fn atom_mut_ptr(&mut self) -> *mut Atom {
        self.top.atoms.as_mut_ptr()
    }

    fn pos_mut_ptr(&mut self) -> *mut Pos {
        self.st.coords.as_mut_ptr()
    }
}

impl NonAtomPosAnalysisMut for SubSystemMut<'_> {
    fn st_ref_mut(&mut self) -> &mut State {
        &mut self.st
    }

    fn top_ref_mut(&mut self) -> &mut Topology {
        &mut self.top
    }
}

//================================================
/// Read-write subsystem for non-blocking parallel access to atoms and posisitons
/// Doesn't have access to shared fields such as box and bonds.
//================================================
pub struct SubSystemParMut<'a> {
    pos_ptr: *mut Pos,
    atom_ptr: *mut Atom,
    index: &'a SVec,
}
unsafe impl Sync for SubSystemParMut<'_> {}
unsafe impl Send for SubSystemParMut<'_> {}

impl LenProvider for SubSystemParMut<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SubSystemParMut<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl AtomPosAnalysisMut for SubSystemParMut<'_> {
    fn atom_mut_ptr(&mut self) -> *mut Atom {
        self.atom_ptr
    }

    fn pos_mut_ptr(&mut self) -> *mut Pos {
        self.pos_ptr
    }
}

//================================================
/// Umbrella trait for implementing read-only analysis traits
trait AtomPosAnalysis: LenProvider + IndexProvider + Sized {
    // Raw pointer to the beginning of array of atoms
    fn atom_ptr(&self) -> *const Atom;
    // Raw pointer to the beginning of array of coords
    fn pos_ptr(&self) -> *const Pos;
}

trait NonAtomPosAnalysis: Sized {
    fn top_ref(&self) -> &Topology;
    fn st_ref(&self) -> &State;
}
//================================================


//================================================
/// Umbrella trait for implementing read-write analysis traits
trait AtomPosAnalysisMut: LenProvider + IndexProvider + Sized {
    // Raw pointer to the beginning of array of atoms
    fn atom_mut_ptr(&mut self) -> *mut Atom;
    // Raw pointer to the beginning of array of coords
    fn pos_mut_ptr(&mut self) -> *mut Pos;
}

trait NonAtomPosAnalysisMut: Sized {
    fn top_ref_mut(&mut self) -> &mut Topology;
    fn st_ref_mut(&mut self) -> &mut State;
}
//================================================



//═══════════════════════════════════════════════════════════
//  Blanket trait implementations 
//═══════════════════════════════════════════════════════════


//██████  Immutable analysis traits

impl<T: AtomPosAnalysis> PosIterProvider for T {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        unsafe {
            self.iter_index()
                .map(|i| &*self.pos_ptr().add(i))
        }
    }
}

impl<T: AtomPosAnalysis> AtomIterProvider for T {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        unsafe {
            self.iter_index()
                .map(|i| &*self.atom_ptr().add(i))
        }
    }
}

impl<T: AtomPosAnalysis> RandomPosProvider for T {
    unsafe fn get_pos_unchecked(&self, i: usize) -> &Pos {
        let ind = self.get_index_unchecked(i);
        &*self.pos_ptr().add(ind)
    }
}

impl<T: AtomPosAnalysis> RandomAtomProvider for T {
    unsafe fn get_atom_unchecked(&self, i: usize) -> &Atom {
        let ind = self.get_index_unchecked(i);
        &*self.atom_ptr().add(ind)
    }
}

impl<T: AtomPosAnalysis> RandomParticleProvider for T {
    unsafe fn get_particle_unchecked(&self, i: usize) -> Particle<'_> {
        let ind = self.get_index_unchecked(i);
        Particle {
            id: ind,
            atom: &*self.atom_ptr().add(ind),
            pos: &*self.pos_ptr().add(ind),
        }
    }
}

// Non atom-pos:

impl<T: NonAtomPosAnalysis> BoxProvider for T {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.st_ref().get_box()
    }
}

impl<T: NonAtomPosAnalysis> TimeProvider for T {
    fn get_time(&self) -> f32 {
        self.st_ref().time
    }
}

impl<T: NonAtomPosAnalysis> RandomMoleculeProvider for T {
    fn num_molecules(&self) -> usize {
        self.top_ref().num_molecules()
    }

    unsafe fn get_molecule_unchecked(&self, i: usize) -> &[usize; 2] {
        self.top_ref().get_molecule_unchecked(i)
    }
}

impl<T: NonAtomPosAnalysis> MoleculeIterProvider for T {
    fn iter_molecules(&self) -> impl Iterator<Item = &[usize; 2]> {
        self.top_ref().iter_molecules()
    }
}

impl<T: NonAtomPosAnalysis> RandomBondProvider for T {
    fn num_bonds(&self) -> usize {
        self.top_ref().num_bonds()
    }

    unsafe fn get_bond_unchecked(&self, i: usize) -> &[usize; 2] {
        self.top_ref().get_bond_unchecked(i)
    }
}

impl<T: NonAtomPosAnalysis> BondIterProvider for T {
    fn iter_bonds(&self) -> impl Iterator<Item = &[usize; 2]> {
        self.top_ref().iter_bonds()
    }
}

//██████  Measure traits

impl<T: AtomPosAnalysis> MeasurePos for T {}
impl<T: AtomPosAnalysis + NonAtomPosAnalysis> MeasurePeriodic for T {}
impl<T: AtomPosAnalysis> MeasureMasses for T {}
impl<T: AtomPosAnalysis> MeasureRandomAccess for T {}

//██████  Mutable analysis traits

impl<T: AtomPosAnalysisMut> RandomAtomMutProvider for T {
    unsafe fn get_atom_mut_unchecked(&mut self, i: usize) -> &mut Atom {
        let ind = self.get_index_unchecked(i);
        &mut *self.atom_mut_ptr().add(ind)
    }
}

impl<T: AtomPosAnalysisMut> AtomIterMutProvider for T {
    fn iter_atoms_mut(&mut self) -> impl AtomMutIterator<'_> {
        let mut i = 0;
        let func = move || {
            if i < self.len() {
                let ind = unsafe{self.get_index_unchecked(i)};
                let res = unsafe{&mut *self.atom_mut_ptr().add(ind)};
                i += 1;
                Some(res)
            } else {
                None
            }
        };
        std::iter::from_fn(func)
    }
}

impl<T: AtomPosAnalysisMut> PosIterMutProvider for T {
    fn iter_pos_mut(&mut self) -> impl PosMutIterator<'_> {
        let mut i = 0;
        let func = move || {
            if i < self.len() {
                let ind = unsafe{self.get_index_unchecked(i)};
                let res = unsafe{&mut *self.pos_mut_ptr().add(ind)};
                i += 1;
                Some(res)
            } else {
                None
            }
        };
        std::iter::from_fn(func)
    }
}

impl<T: AtomPosAnalysisMut> RandomPosMutProvider for T {
    unsafe fn get_pos_mut_unchecked(&mut self, i: usize) -> &mut Pos {
        let ind = self.get_index_unchecked(i);
        &mut *self.pos_mut_ptr().add(ind)
    }
}

impl<T: AtomPosAnalysisMut> RandomParticleMutProvider for T {
    unsafe fn get_particle_mut_unchecked(&mut self, i: usize) -> ParticleMut {
        let ind = self.get_index_unchecked(i);
        ParticleMut {
            id: ind,
            atom: unsafe{&mut *self.atom_mut_ptr().add(ind)},
            pos: unsafe{&mut *self.pos_mut_ptr().add(ind)},
        }
    }
}

// Non-atompos
impl<T: NonAtomPosAnalysisMut> TimeMutProvider for T {
    fn set_time(&mut self, t: f32) {
        self.st_ref_mut().time = t;
    }
}

impl<T: NonAtomPosAnalysisMut> BoxMutProvider for T {
    fn get_box_mut(&mut self) -> Option<&mut PeriodicBox> {
        self.st_ref_mut().pbox.as_mut()
    }
}

//██████  Modify traits

impl<T: AtomPosAnalysisMut> ModifyPos for T {}
impl<T: AtomPosAnalysisMut + NonAtomPosAnalysisMut> ModifyPeriodic for T {}
impl<T: AtomPosAnalysisMut + NonAtomPosAnalysisMut> ModifyRandomAccess for T {}

//====================================================================



#[cfg(test)]
mod tests {
    use super::System;
    use crate::{
        core::{system::Sel, AtomIterMutProvider, AtomIterProvider},
        io::FileHandler,
    };

    #[test]
    fn test1() -> anyhow::Result<()> {
        let mut h = FileHandler::open("tests/protein.pdb").unwrap();
        let (top, st) = h.read().unwrap();
        let mut sys = System::new(top, st);
        let sel1 = Sel::from_vec(vec![1, 2, 6, 7])?;
        
        for at in sel1.bind(&sys)?.iter_atoms() {
            println!("{} {}", at.name, at.resname);
        }

        let mut lock = sel1.bind_mut(&mut sys)?;
        for at in lock.iter_atoms_mut() {
            at.bfactor += 1.0;
        }

        let lock2 = sel1.bind(&sys)?;
        for at in lock2.iter_atoms() {
            println!("{} ", at.bfactor);
        }

        //drop(lock);

        Ok(())
    }
}
