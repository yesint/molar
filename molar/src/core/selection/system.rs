use crate::prelude::*;
use std::path::Path;

//================================================
/// System that stores Topology and State
//================================================
#[derive(Debug, Default)]
pub struct System {
    pub(super) top: Topology,
    pub(super) st: State,
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
        let last = unsafe { *ind.0.get_unchecked(ind.0.len() - 1) };
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
        let last = unsafe { *ind.0.get_unchecked(ind.0.len() - 1) };
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

    pub fn bind_par<'a>(&'a self, par: &'a ParSplit) -> Result<ParSplitBound<'a>, SelectionError> {
        if par.max_index < self.top.len() {
            Ok(ParSplitBound {
                sys: self,
                indexes: &par.selections,
            })
        } else {
            Err(SelectionError::OutOfBounds(par.max_index, self.top.len()))
        }
    }

    pub fn bind_par_mut<'a>(
        &'a mut self,
        par: &'a ParSplit,
    ) -> Result<ParSplitBoundMut<'a>, SelectionError> {
        if par.max_index < self.top.len() {
            Ok(ParSplitBoundMut {
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
