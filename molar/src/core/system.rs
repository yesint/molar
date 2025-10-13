use std::path::Path;

use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator};

use crate::prelude::*;

//===========================================================================
/// Selection is just an index detached from the Topology and State.
/// Has to be bound to System before doing something with the subset of atoms.
/// It is guaranteed to be non-empty so no run-time checks for this
/// are needed when binding to system
//===========================================================================
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

    /// Creates a string in Gromacs index format representing self.
    fn as_gromacs_ndx_str(&self, name: impl AsRef<str>) -> String {
        use itertools::Itertools;
        let name = name.as_ref();
        let mut s = format!("[ {} ]\n", name);
        for chunk in &self.0.iter().chunks(15) {
            let line: String = chunk.map(|i| (i + 1).to_string()).join(" ");
            s.push_str(&line);
            s.push('\n');
        }
        s
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
        Ok(Sel(def.into_sel_index(
            self.top,
            self.st,
            Some(self.index),
        )?))
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
        Ok(Sel(def.into_sel_index(
            self.top,
            self.st,
            Some(self.index),
        )?))
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

//============================================================================
/// Collection of non-overlapping selections that could be mutated in parallel
/// Selections don't have access to shared fields such as box and bonds.
//============================================================================
pub struct ParSplit {
    selections: Vec<Sel>,
    max_index: usize,
}

/// Bound parallel split
pub struct ParSplitBound<'a>(Vec<SubSystemParMut<'a>>);

impl ParSplit {
    /// Binds split to System for read-write access
    pub fn bind_mut<'a>(
        &'a self,
        sys: &'a mut System,
    ) -> Result<ParSplitBound<'a>, SelectionError> {
        // The cost calling is just one comparison
        if self.max_index < sys.top.len() {
            // Bind each selection in mutable parallel way
            let bound_sels: Vec<_> = self
                .selections
                .iter()
                .map(|sel| SubSystemParMut {
                    index: &sel.0,
                    pos_ptr: sys.st.coords.as_mut_ptr(),
                    atom_ptr: sys.top.atoms.as_mut_ptr(),
                })
                .collect();
            Ok(ParSplitBound(bound_sels))
        } else {
            Err(SelectionError::OutOfBounds(0, 0))
        }
    }
}

impl<'a> ParSplitBound<'a> {
    /// Returns parallel iterator over stored parallel selections.
    pub fn par_iter(&'a self) -> rayon::slice::Iter<'a, SubSystemParMut<'a>> {
        self.0.par_iter()
    }

    /// Returns parallel mutable iterator over stored parallel selections.
    pub fn par_iter_mut(&'a mut self) -> rayon::slice::IterMut<'a, SubSystemParMut<'a>> {
        self.0.par_iter_mut()
    }

    /// Returns serial iterator over stored parallel selections.
    pub fn iter(&'a self) -> impl Iterator<Item = &'a SubSystemParMut<'a>> {
        self.0.iter()
    }

    /// Returns serial mutable iterator over stored parallel selections.
    pub fn iter_mut(&'a mut self) -> impl Iterator<Item = &'a mut SubSystemParMut<'a>> {
        self.0.iter_mut()
    }
}

//================================================
/// Umbrella trait for implementing read-only analysis traits involving atoms and positions
//================================================
trait AtomPosAnalysis: LenProvider + IndexProvider + Sized {
    // Raw pointer to the beginning of array of atoms
    fn atom_ptr(&self) -> *const Atom;
    // Raw pointer to the beginning of array of coords
    fn pos_ptr(&self) -> *const Pos;

    /// Creates a parallel split based on provided closure.
    /// A closure takes a [Particle] and returns a distinct value for each piece.
    /// New selection is created whenever new return value differes from the previous one.
    fn split_par<F, R>(&self, func: F) -> Result<ParSplit, SelectionError>
    where
        F: Fn(Particle) -> Option<R>,
        R: Default + PartialOrd,
        Self: Sized,
    {
        let selections: Vec<Sel> = self.split_iter(func).collect();

        if selections.is_empty() {
            return Err(SelectionError::EmptySplit);
        }

        // Compute last index
        let max_index = selections
            .iter()
            .map(|sel| *sel.0.last().unwrap())
            .max()
            .unwrap();

        Ok(ParSplit {
            selections,
            max_index,
        })
    }

    // Internal splitting function
    fn split_iter<RT, F>(&self, func: F) -> impl Iterator<Item = Sel>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(Particle) -> Option<RT>,
    {
        let mut cur_val = RT::default();
        let mut cur = 0usize;

        let next_fn = move || {
            let mut index = Vec::<usize>::new();

            while cur < self.len() {
                let p = unsafe { self.get_particle_unchecked(cur) };
                let id = p.id;
                let val = func(p);

                if let Some(val) = val {
                    if val == cur_val {
                        // Current selection continues. Add current index
                        index.push(id);
                    } else if index.is_empty() {
                        // The very first id is not default, this is Ok, add index
                        // and update self.id
                        cur_val = val;
                        index.push(id);
                    } else {
                        // The end of current selection
                        cur_val = val; // Update val for the next selection
                        return Some(Sel(unsafe { SVec::from_sorted(index) }));
                    }
                }
                // Next particle
                cur += 1;
            }

            // Return any remaining index as last selection
            if !index.is_empty() {
                return Some(Sel(unsafe { SVec::from_sorted(index) }));
            }

            // If we are here stop iterating
            None
        };

        std::iter::from_fn(next_fn)
    }

    fn split_resindex_iter(&self) -> impl Iterator<Item = Sel> {
        self.split_iter(|p| Some(p.atom.resindex))
    }

    /// Creates an "expanded" selection that includes all atoms with the same attributes as encountered in current selection.
    /// Provided closure takes an [Atom] and returns an "attribute" value (for example resid, resindex, chain).
    /// All atoms in the parent [Topology] that have the same attribute as any atom of the current selection will be selected.
    /// Functionally this method is equivalent to "same <whatever> as" selections in VMD or Gromacs.
    /// ## Example - select whole residues
    /// ```ignore
    /// let sel = system.select("name CA and resid 1:10")?;
    /// let whole_res = sel.whole_attr(|atom| &atom.resid);
    /// ```
    fn whole_attr<T>(&self, attr_fn: fn(&Atom) -> &T) -> Sel
    where
        T: Eq + std::hash::Hash + Copy,
    {
        // Collect all properties from the inner
        let mut properties = std::collections::HashSet::<T>::new();
        for at in self.iter_atoms() {
            properties.insert(*attr_fn(at));
        }

        let mut ind = vec![];
        // Loop over all atoms with the same property
        for (i, at) in self.iter_atoms().enumerate() {
            let cur_prop = attr_fn(at);
            if properties.contains(cur_prop) {
                ind.push(i);
            }
        }

        Sel(unsafe { SVec::from_sorted(ind) })
    }

    /// Selects whole residiues present in the current selection (in terms of resindex)
    fn whole_residues(&self) -> Sel {
        self.whole_attr(|at| &at.resindex)
    }

    /// Selects whole chains present in the current selection
    fn whole_chains(&self) -> Sel {
        self.whole_attr(|at| &at.chain)
    }

    /// Computes the Solvet Accessible Surface Area (SASA).
    fn sasa(&self) -> SasaResults {
        molar_powersasa::compute_sasa(
            self.len(),
            0.14,
            |i| unsafe {
                let ind = self.get_index_unchecked(i);
                self.pos_ptr().add(ind) as *mut f32
            },
            |i: usize| self.get_particle(i).unwrap().atom.vdw(),
        )
    }
}

//================================================
/// Umbrella trait for implementing read-only analysis traits NOT involving atoms and positions
//================================================
trait NonAtomPosAnalysis: LenProvider + IndexProvider + Sized {
    fn top_ref(&self) -> &Topology;
    fn st_ref(&self) -> &State;

    /// Splits by molecule and returns an iterator over them. 
    /// If molecule is only partially contained in self then only this part is returned (molecules are clipped).
    /// If there are no molecules in [Topology] return an empty iterator.
    fn split_mol_iter(&self) -> impl Iterator<Item = Sel>
    where
        Self: Sized,
    {
        let n = self.len();
        // Iterate over molecules and find those inside selection
        let first = unsafe{self.get_index_unchecked(0)};
        let last = unsafe{self.get_index_unchecked(n-1)};

        let mut molid = 0;

        let next_fn = move || {
            if self.top_ref().num_molecules() == 0 {
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
            }.map(|r| Sel(unsafe{SVec::from_sorted(r.into_iter().collect())}));

            molid += 1;
            res
        };

        std::iter::from_fn(next_fn)
    }
}
//================================================

//================================================
/// Umbrella trait for implementing read-write analysis traits involving atoms and positions
//================================================
trait AtomPosAnalysisMut: LenProvider + IndexProvider + Sized {
    // Raw pointer to the beginning of array of atoms
    fn atom_mut_ptr(&mut self) -> *mut Atom;
    // Raw pointer to the beginning of array of coords
    fn pos_mut_ptr(&mut self) -> *mut Pos;
}

//================================================
/// Umbrella trait for implementing read-write analysis traits NOT involving atoms and positions
//================================================
trait NonAtomPosAnalysisMut: Sized {
    fn top_ref_mut(&mut self) -> &mut Topology;
    fn st_ref_mut(&mut self) -> &mut State;
}

//═══════════════════════════════════════════════════════════
//  Blanket trait implementations
//═══════════════════════════════════════════════════════════

//██████  Immutable analysis traits involving atoms and positions only

impl<T: AtomPosAnalysis> PosIterProvider for T {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        unsafe { self.iter_index().map(|i| &*self.pos_ptr().add(i)) }
    }
}

impl<T: AtomPosAnalysis> AtomIterProvider for T {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        unsafe { self.iter_index().map(|i| &*self.atom_ptr().add(i)) }
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

//██████  Immutable analysis traits NOT involving atoms and positions

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

//██████  Mutable analysis traits involving only atoms and positions

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
                let ind = unsafe { self.get_index_unchecked(i) };
                let res = unsafe { &mut *self.atom_mut_ptr().add(ind) };
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
                let ind = unsafe { self.get_index_unchecked(i) };
                let res = unsafe { &mut *self.pos_mut_ptr().add(ind) };
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
            atom: unsafe { &mut *self.atom_mut_ptr().add(ind) },
            pos: unsafe { &mut *self.pos_mut_ptr().add(ind) },
        }
    }
}

//██████  Mutable analysis traits NOT involving atoms and positions

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

//██████  Measure traits

impl<T: AtomPosAnalysis> MeasurePos for T {}
impl<T: AtomPosAnalysis + NonAtomPosAnalysis> MeasurePeriodic for T {}
impl<T: AtomPosAnalysis> MeasureMasses for T {}
impl<T: AtomPosAnalysis> MeasureRandomAccess for T {}

//██████  Modify traits

impl<T: AtomPosAnalysisMut> ModifyPos for T {}
impl<T: AtomPosAnalysisMut + NonAtomPosAnalysisMut> ModifyPeriodic for T {}
impl<T: AtomPosAnalysisMut + NonAtomPosAnalysisMut> ModifyRandomAccess for T {}

//====================================================================================

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
