use crate::prelude::*;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator};
use std::path::Path;
use thiserror::Error;

//===========================================================================
/// Selection index detached from any [System].
/// It is guaranteed to be non-empty.
//===========================================================================
pub struct SelIndex(pub(crate) SVec);

impl SelIndex {
    pub fn from_vec(index: Vec<usize>) -> Result<Self, SelectionError> {
        if index.is_empty() {
            Err(SelectionError::EmptySlice)
        } else {
            Ok(Self(SVec::from_unsorted(index)))
        }
    }

    pub fn into_svec(self) -> SVec {
        self.0
    }

    pub fn from_svec(index: SVec) -> Result<Self, SelectionError> {
        if index.is_empty() {
            Err(SelectionError::EmptySlice)
        } else {
            Ok(Self(index))
        }
    }

    fn from_iter(iter: impl Iterator<Item = usize>) -> Result<Self, SelectionError> {
        let index: Vec<_> = iter.collect();
        Ok(Self::from_vec(index)?)
    }
}

impl LenProvider for SelIndex {
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl IndexProvider for SelIndex {
    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.0.iter().cloned()
    }

    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.0.get_unchecked(i)
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
    pub fn new(top: Topology, st: State) -> Result<Self, SelectionError> {
        check_topology_state_sizes(&top, &st)?;
        Ok(Self { top, st })
    }

    pub fn from_file(fname: impl AsRef<Path>) -> Result<Self, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        Ok(Self::new(top, st)?)
    }

    /// Create new detached selection index based on provided definition.
    pub fn select_as_index(&self, def: impl SelectionDef) -> Result<SelIndex, SelectionError> {
        Ok(SelIndex(def.into_sel_index(&self.top, &self.st, None)?))
    }

    /// Create all detached
    pub fn select_all_as_index(&self) -> SelIndex {
        SelIndex(unsafe { SVec::from_sorted((0..self.len()).into_iter().collect()) })
    }

    /// Create new selection based on provided definition.
    pub fn select(&self, def: impl SelectionDef) -> Result<Sel<'_>, SelectionError> {
        Ok(Sel {
            sys: self,
            index: def.into_sel_index(&self.top, &self.st, None)?,
        })
    }

    pub fn select_all(&self) -> Sel<'_> {
        Sel {
            sys: self,
            index: unsafe { SVec::from_sorted((0..self.len()).into_iter().collect()) },
        }
    }

    pub fn select_mut(&mut self, def: impl SelectionDef) -> Result<SelMut<'_>, SelectionError> {
        let index = def.into_sel_index(&self.top, &self.st, None)?;
        Ok(SelMut { sys: self, index })
    }

    pub fn select_all_mut(&mut self) -> SelMut<'_> {
        let index = unsafe { SVec::from_sorted((0..self.len()).into_iter().collect()) };
        SelMut { sys: self, index }
    }

    /// Binds detached selection index to make borrowed selection.
    /// `ind`  is not consumed.
    pub fn bind<'a>(&'a self, ind: &'a SelIndex) -> Result<SelBorrowing<'a>, SelectionError> {
        // No need to check for empty index since it's guaranteed to be non-empty
        let last = unsafe { *ind.0.get_unchecked(ind.0.len()) };
        if last >= self.top.len() {
            Err(SelectionError::IndexValidation(
                *ind.0.first().unwrap(),
                last,
                self.top.len() - 1,
            ))
        } else {
            Ok(SelBorrowing {
                sys: self,
                index: ind.0.as_slice(),
            })
        }
    }

    /// Mutably binds detached selection index to make borrowed selection.
    /// `ind`  is not consumed.
    pub fn bind_mut<'a>(
        &'a mut self,
        ind: &'a SelIndex,
    ) -> Result<SelBorrowingMut<'a>, SelectionError> {
        // No need to check for empty index since it's guaranteed to be non-empty
        let last = unsafe { *ind.0.get_unchecked(ind.0.len()) };
        if last >= self.top.len() {
            Err(SelectionError::IndexValidation(
                *ind.0.first().unwrap(),
                last,
                self.top.len() - 1,
            ))
        } else {
            Ok(SelBorrowingMut {
                sys: self,
                index: ind.0.as_slice(),
            })
        }
    }

    // TODO: This function allocated a vector each time when called!
    // This is inefficient. Think how to pre-allocate.
    pub fn bind_par<'a>(
        &'a mut self,
        par: &'a ParSplit,
    ) -> Result<ParSplitBound<'a>, SelectionError> {
        if par.max_index < self.top.len() {
            Ok(ParSplitBound {
                sys: self,
                indexes: &par.selections,
            })
        } else {
            Err(SelectionError::OutOfBounds(par.max_index, self.top.len()))
        }
    }

    pub fn split<'a, RT, F>(&'a self, func: F) -> impl Iterator<Item = Sel<'a>>
    where
        RT: Default + std::cmp::PartialEq + 'a,
        F: Fn(Particle) -> Option<RT> + 'a,
    {
        self.split_as_index(func).map(|sel| Sel {
            sys: self,
            index: sel.0,
        })
    }

    pub fn split_resindex(&self) -> impl Iterator<Item = Sel<'_>> {
        self.split(|p| Some(p.atom.resindex))
    }

    pub fn set_state(&mut self, st: State) -> Result<State, SelectionError> {
        if self.len() != st.len() {
            Err(SelectionError::IncompatibleState)
        } else {
            Ok(std::mem::replace(&mut self.st, st))
        }
    }

    pub fn set_topology(&mut self, top: Topology) -> Result<Topology, SelectionError> {
        if self.len() != top.len() {
            Err(SelectionError::IncompatibleTopology)
        } else {
            Ok(std::mem::replace(&mut self.top, top))
        }
    }

    pub fn release(self) -> (Topology, State) {
        (self.top, self.st)
    }

    //===============
    // Modifying
    //===============

    /// Append selection derived from self
    pub fn append_self_index(&mut self, sel: &SelIndex) -> Result<SelIndex, SelectionError> {
        let old_last = self.len() - 1;
        let last_ind = unsafe { *sel.0.get_unchecked(sel.0.len() - 1) };
        if last_ind >= self.top.len() {
            return Err(SelectionError::IndexValidation(
                sel.0[0],
                last_ind,
                self.top.len(),
            ))?;
        }
        let pos: Vec<_> = sel.0.iter().map(|i| &self.st.coords[*i]).cloned().collect();
        let atoms: Vec<_> = sel.0.iter().map(|i| &self.top.atoms[*i]).cloned().collect();
        self.st.add_coords(pos.into_iter());
        self.top.add_atoms(atoms.into_iter());
        Ok(SelIndex::from_iter(old_last + 1..self.len())?)
    }

    pub fn append_atoms_coords<'a>(
        &mut self,
        atoms: impl AtomIterator<'a>,
        coords: impl PosIterator<'a>,
    ) -> Result<SelIndex, SelectionError> {
        let old_last = self.len() - 1;
        self.st.add_coords(coords.cloned());
        self.top.add_atoms(atoms.cloned());
        check_topology_state_sizes(&self.top, &self.st)?;
        Ok(SelIndex::from_iter(old_last + 1..self.len())?)
    }

    pub fn append(
        &mut self,
        data: &(impl AtomIterProvider + PosIterProvider),
    ) -> Result<SelIndex, SelectionError> {
        let old_last = self.len() - 1;
        self.st.add_coords(data.iter_pos().cloned());
        self.top.add_atoms(data.iter_atoms().cloned());
        check_topology_state_sizes(&self.top, &self.st)?;
        Ok(SelIndex::from_iter(old_last + 1..self.len())?)
    }

    pub fn remove(
        &mut self,
        removed: impl Iterator<Item = usize> + Clone,
    ) -> Result<(), BuilderError> {
        self.st.remove_coords(removed.clone())?;
        self.top.remove_atoms(removed)?;
        Ok(())
    }

    pub fn set_box_from(&mut self, src: &impl BoxProvider) {
        self.st.pbox = src.get_box().cloned();
    }

    pub fn multiply_periodically(&mut self, nbox: [usize; 3]) -> Result<(), SelectionError> {
        if self.get_box().is_none() {
            return Err(PeriodicBoxError::NoPbc)?;
        }
        let m = self.require_box()?.get_matrix();
        let all = self.select_all_as_index();
        for x in 0..=nbox[0] {
            for y in 0..=nbox[1] {
                for z in 0..=nbox[2] {
                    if x == 0 && y == 0 && z == 0 {
                        continue;
                    }
                    let added = self.append_self_index(&all)?;
                    let shift =
                        m.column(0) * x as f32 + m.column(1) * y as f32 + m.column(2) * z as f32;
                    self.select_mut(added)?.translate(&shift);
                }
            }
        }
        // Scale the box
        self.get_box_mut().unwrap().scale_vectors([
            nbox[0] as f32,
            nbox[1] as f32,
            nbox[2] as f32,
        ])?;

        // Re-assign resindex
        self.top.assign_resindex();
        Ok(())
    }

    pub fn assign_resindex(&mut self) {
        self.top.assign_resindex();
    }
}

impl TopologyWrite for System {}
impl StateWrite for System {}
impl TopologyStateWrite for System {}

impl LenProvider for System {
    fn len(&self) -> usize {
        self.top.atoms.len()
    }
}

impl IndexProvider for System {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        i
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        (0..self.len()).into_iter()
    }
}

impl AtomPosAnalysis for System {
    fn atoms_ptr(&self) -> *const Atom {
        self.top.atoms.as_ptr()
    }

    fn coords_ptr(&self) -> *const Pos {
        self.st.coords.as_ptr()
    }
}

impl AtomPosAnalysisMut for System {
    fn atoms_ptr_mut(&mut self) -> *mut Atom {
        self.top.atoms.as_mut_ptr()
    }

    fn coords_ptr_mut(&mut self) -> *mut Pos {
        self.st.coords.as_mut_ptr()
    }
}

impl NonAtomPosAnalysis for System {
    fn top_ptr(&self) -> *const Topology {
        &self.top
    }

    fn st_ptr(&self) -> *const State {
        &self.st
    }
}

impl NonAtomPosAnalysisMut for System {
    fn top_ptr_mut(&mut self) -> *mut Topology {
        &mut self.top
    }

    fn st_ptr_mut(&mut self) -> *mut State {
        &mut self.st
    }
}

//================================================
/// Read only subsystem that owns its index
/// Implements only read-only analysis traits
//================================================

#[derive(Clone, Debug)]
pub struct Sel<'a> {
    sys: &'a System,
    index: SVec,
}

impl Sel<'_> {
    /// Create new unbound sub-selection based on provided definition.
    pub fn select_as_index(&self, def: impl SelectionDef) -> Result<SelIndex, SelectionError> {
        Ok(SelIndex(def.into_sel_index(
            &self.sys.top,
            &self.sys.st,
            Some(self.index.as_slice()),
        )?))
    }

    /// Create new bound sub-selection based on provided definition.
    pub fn select(&self, def: impl SelectionDef) -> Result<Self, SelectionError> {
        Ok(Self {
            sys: self.sys,
            index: def.into_sel_index(&self.sys.top, &self.sys.st, Some(self.index.as_slice()))?,
        })
    }

    pub fn into_index(self) -> SelIndex {
        SelIndex(self.index)
    }

    pub fn clone_index(&self) -> SelIndex {
        SelIndex(self.index.clone())
    }

    pub fn split<'a, RT, F>(&'a self, func: F) -> impl Iterator<Item = Self> + 'a
    where
        RT: Default + std::cmp::PartialEq + 'a,
        F: Fn(Particle) -> Option<RT> + 'a,
    {
        self.split_as_index(func).map(|sel| Self {
            sys: &self.sys,
            index: sel.0,
        })
    }

    pub fn split_resindex(&self) -> impl Iterator<Item = Self> + '_ {
        self.split(|p| Some(p.atom.resindex))
    }
}

impl LenProvider for Sel<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for Sel<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl AtomPosAnalysis for Sel<'_> {
    fn atoms_ptr(&self) -> *const Atom {
        self.sys.top.atoms.as_ptr()
    }

    fn coords_ptr(&self) -> *const Pos {
        self.sys.st.coords.as_ptr()
    }
}

impl NonAtomPosAnalysis for Sel<'_> {
    fn top_ptr(&self) -> *const Topology {
        &self.sys.top
    }

    fn st_ptr(&self) -> *const State {
        &self.sys.st
    }
}

impl TopologyWrite for Sel<'_> {}
impl StateWrite for Sel<'_> {}
impl TopologyStateWrite for Sel<'_> {}

//================================================
/// Read-write bound subsystem having access to
/// all fields of Topology and State
//================================================
pub struct SelMut<'a> {
    sys: &'a mut System,
    index: SVec,
}

impl SelMut<'_> {
    /// Create new sub-selection based on provided definition.
    pub fn unbind(self) -> SelIndex {
        SelIndex(self.index)
    }

    pub fn into_index(self) -> SelIndex {
        SelIndex(self.index)
    }

    pub fn clone_index(&self) -> SelIndex {
        SelIndex(self.index.clone())
    }
}

impl LenProvider for SelMut<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SelMut<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl AtomPosAnalysis for SelMut<'_> {
    fn atoms_ptr(&self) -> *const Atom {
        self.sys.top.atoms.as_ptr()
    }

    fn coords_ptr(&self) -> *const Pos {
        self.sys.st.coords.as_ptr()
    }
}

impl AtomPosAnalysisMut for SelMut<'_> {
    fn atoms_ptr_mut(&mut self) -> *mut Atom {
        self.sys.top.atoms.as_mut_ptr()
    }

    fn coords_ptr_mut(&mut self) -> *mut Pos {
        self.sys.st.coords.as_mut_ptr()
    }
}

impl NonAtomPosAnalysis for SelMut<'_> {
    fn st_ptr(&self) -> *const State {
        &self.sys.st
    }

    fn top_ptr(&self) -> *const Topology {
        &self.sys.top
    }
}

impl NonAtomPosAnalysisMut for SelMut<'_> {
    fn st_ptr_mut(&mut self) -> *mut State {
        &mut self.sys.st
    }

    fn top_ptr_mut(&mut self) -> *mut Topology {
        &mut self.sys.top
    }
}

//================================================
/// Read only selection that borrows its index
/// Implements only read-only analysis traits
//================================================

#[derive(Clone, Debug)]
pub struct SelBorrowing<'a> {
    sys: &'a System,
    index: &'a [usize],
}

impl SelBorrowing<'_> {
    /// Create new owned sub-selection based on provided definition.
    pub fn select(&self, def: impl SelectionDef) -> Result<Sel<'_>, SelectionError> {
        let index = def.into_sel_index(&self.sys.top, &self.sys.st, Some(self.index))?;
        Ok(Sel {
            index,
            sys: &self.sys,
        })
    }

    pub fn clone_index(&self) -> SelIndex {
        SelIndex(SVec::from_iter(self.index.iter().cloned()))
    }
}

impl LenProvider for SelBorrowing<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SelBorrowing<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl AtomPosAnalysis for SelBorrowing<'_> {
    fn atoms_ptr(&self) -> *const Atom {
        self.sys.top.atoms.as_ptr()
    }

    fn coords_ptr(&self) -> *const Pos {
        self.sys.st.coords.as_ptr()
    }
}

impl NonAtomPosAnalysis for SelBorrowing<'_> {
    fn top_ptr(&self) -> *const Topology {
        &self.sys.top
    }

    fn st_ptr(&self) -> *const State {
        &self.sys.st
    }
}

impl TopologyWrite for SelBorrowing<'_> {}
impl StateWrite for SelBorrowing<'_> {}
impl TopologyStateWrite for SelBorrowing<'_> {}

//================================================
/// Read-write selection that borrows its index
//================================================
pub struct SelBorrowingMut<'a> {
    sys: &'a mut System,
    index: &'a [usize],
}

impl SelBorrowingMut<'_> {
    /// Create new sub-selection based on provided definition.
    pub fn select(&self, def: impl SelectionDef) -> Result<Sel<'_>, SelectionError> {
        let index = def.into_sel_index(&self.sys.top, &self.sys.st, Some(self.index))?;
        Ok(Sel {
            index,
            sys: &self.sys,
        })
    }

    pub fn clone_index(&self) -> SelIndex {
        SelIndex(SVec::from_iter(self.index.iter().cloned()))
    }
}

impl LenProvider for SelBorrowingMut<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SelBorrowingMut<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl AtomPosAnalysis for SelBorrowingMut<'_> {
    fn atoms_ptr(&self) -> *const Atom {
        self.sys.top.atoms.as_ptr()
    }

    fn coords_ptr(&self) -> *const Pos {
        self.sys.st.coords.as_ptr()
    }
}

impl AtomPosAnalysisMut for SelBorrowingMut<'_> {
    fn atoms_ptr_mut(&mut self) -> *mut Atom {
        self.sys.top.atoms.as_mut_ptr()
    }

    fn coords_ptr_mut(&mut self) -> *mut Pos {
        self.sys.st.coords.as_mut_ptr()
    }
}

impl NonAtomPosAnalysis for SelBorrowingMut<'_> {
    fn st_ptr(&self) -> *const State {
        &self.sys.st
    }

    fn top_ptr(&self) -> *const Topology {
        &self.sys.top
    }
}

impl NonAtomPosAnalysisMut for SelBorrowingMut<'_> {
    fn st_ptr_mut(&mut self) -> *mut State {
        &mut self.sys.st
    }

    fn top_ptr_mut(&mut self) -> *mut Topology {
        &mut self.sys.top
    }
}

//================================================
/// Read-write subsystem for non-blocking parallel access to atoms and posisitons
/// Doesn't have access to shared fields such as box and bonds.
//================================================
pub struct SelPar<'a> {
    coords_ptr: *mut Pos,
    atoms_ptr: *mut Atom,
    index: &'a [usize],
}

unsafe impl Sync for SelPar<'_> {}
unsafe impl Send for SelPar<'_> {}

impl LenProvider for SelPar<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SelPar<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl AtomPosAnalysis for SelPar<'_> {
    fn atoms_ptr(&self) -> *const Atom {
        self.atoms_ptr
    }

    fn coords_ptr(&self) -> *const Pos {
        self.coords_ptr
    }
}

impl AtomPosAnalysisMut for SelPar<'_> {
    fn atoms_ptr_mut(&mut self) -> *mut Atom {
        self.atoms_ptr
    }

    fn coords_ptr_mut(&mut self) -> *mut Pos {
        self.coords_ptr
    }
}

//============================================================================
/// Collection of non-overlapping selections that could be mutated in parallel
/// Selections don't have access to shared fields such as box and bonds.
//============================================================================
pub struct ParSplit {
    pub(crate) selections: Vec<SelIndex>,
    pub(crate) max_index: usize,
}

/// Bound parallel split
pub struct ParSplitBound<'a> {
    sys: &'a mut System,
    indexes: &'a Vec<SelIndex>,
}

impl<'a> ParSplitBound<'a> {
    /// Returns parallel iterator over parallel selections.
    pub fn par_iter(&'a mut self) -> impl IndexedParallelIterator<Item = SelPar<'a>> {
        self.indexes.par_iter().map(|ind| SelPar {
            index: &ind.0,
            atoms_ptr: self.sys.atoms_ptr() as *mut Atom,
            coords_ptr: self.sys.coords_ptr() as *mut Pos,
        })
    }

    /// Returns serial iterator over parallel selections.
    pub fn iter(&'a self) -> impl Iterator<Item = SelPar<'a>> {
        self.indexes.iter().map(|ind| SelPar {
            index: &ind.0,
            atoms_ptr: self.sys.atoms_ptr() as *mut Atom,
            coords_ptr: self.sys.coords_ptr() as *mut Pos,
        })
    }
}

#[macro_export]
macro_rules! with_sels {
    ($inp:expr, $($sel:ident),+ , $body:block) => {{
        $(let $sel = $sel.bind(&$inp);)+
        $body
    }};
}

//====================================================================================

#[cfg(test)]
mod tests {
    use super::*;
    //use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator};

    #[test]
    fn test1() -> anyhow::Result<()> {
        let mut h = FileHandler::open("tests/protein.pdb").unwrap();
        let (top, st) = h.read().unwrap();
        let mut sys = System::new(top, st)?;
        let sel1 = sys.select(vec![1, 2, 6, 7])?;

        for at in sel1.iter_atoms() {
            println!("{} {}", at.name, at.resname);
        }

        let ind = sel1.into_index();
        for at in sys.bind_mut(&ind)?.iter_atoms_mut() {
            at.bfactor += 1.0;
        }

        // let lock2 = sel1.bind(&sys);
        // for at in lock2.iter_atoms() {
        //     println!("{} ", at.bfactor);
        // }

        //drop(lock);

        Ok(())
    }

    #[test]
    fn test1_1() -> anyhow::Result<()> {
        let mut sys = System::from_file("tests/albumin.pdb")?;

        let par = sys.split_par(|p| {
            if p.atom.resname != "SOL" {
                Some(p.atom.resindex)
            } else {
                None
            }
        })?;

        sys.bind_par(&par)?.par_iter().try_for_each(|mut sel| {
            //println!("{} {}", sel.len(), sel.first_atom().resname);
            if sel.first_atom().resname == "ALA" {
                sel.first_pos_mut().coords[1] = 555.5;
                println!("{}", sel.first_pos());
            } else {
                sel.first_pos_mut().coords[1] = 777.7;
            }
            Ok::<_, SelectionError>(())
        })?;

        // Add serial selection
        let ca = sys.select("name CA")?;
        let cb = sys.select("name CB")?;
        println!("#ca: {} {}", ca.len(), cb.center_of_mass()?);

        for a in cb.iter_atoms().take(10) {
            println!("{}", a.name);
        }

        //Iter serial views
        let b = sys.bind_par(&par)?;
        let serials: Vec<_> = b.iter().collect();
        println!("serial #5: {}", serials[5].first_particle().pos);

        Ok(())
    }
}
