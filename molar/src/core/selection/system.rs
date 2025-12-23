use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator};

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

    /// Create unbound sub-selection
    pub fn sub_select(&self, ind: &Sel, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        Ok(self.try_bind(ind)?.select_as_index(def)?)
    }

    /// Create all detached
    pub fn select_all(&self) -> Sel {
        Sel(unsafe { SVec::from_sorted((0..self.len()).into_iter().collect()) })
    }

    /// Create new bound selection based on provided definition.
    pub fn select_bound(&self, def: impl SelectionDef) -> Result<SelOwnBound<'_>, SelectionError> {
        Ok(SelOwnBound {
            sys: self,
            index: def.into_sel_index(&self.top, &self.st, None)?,
        })
    }

    pub fn select_all_bound(&self) -> SelOwnBound<'_> {
        SelOwnBound {
            sys: self,
            index: unsafe { SVec::from_sorted((0..self.len()).into_iter().collect()) },
        }
    }

    pub fn select_bound_mut(
        &mut self,
        def: impl SelectionDef,
    ) -> Result<SelOwnBoundMut<'_>, SelectionError> {
        let index = def.into_sel_index(&self.top, &self.st, None)?;
        Ok(SelOwnBoundMut { sys: self, index })
    }

    pub fn select_all_bound_mut(&mut self) -> SelOwnBoundMut<'_> {
        let index = unsafe { SVec::from_sorted((0..self.len()).into_iter().collect()) };
        SelOwnBoundMut { sys: self, index }
    }

    /// Binds detached selection index to make borrowed selection.
    /// `sel`  is not consumed.
    pub fn bind<'a>(&'a self, sel: &'a Sel) -> SelBound<'a> {
        // No need to check for empty index since it's guaranteed to be non-empty
        let last = unsafe { *sel.0.get_unchecked(sel.0.len() - 1) };
        if last >= self.top.len() {
            panic!("selection is out of bounds");
        } else {
            SelBound {
                sys: self,
                index: sel.0.as_slice(),
            }
        }
    }

    /// Binds detached selection index to make borrowed selection.
    /// `sel`  is not consumed.
    pub fn try_bind<'a>(&'a self, sel: &'a Sel) -> Result<SelBound<'a>, SelectionError> {
        // No need to check for empty index since it's guaranteed to be non-empty
        let last = unsafe { *sel.0.get_unchecked(sel.0.len() - 1) };
        if last >= self.top.len() {
            Err(SelectionError::IndexValidation(
                *sel.0.first().unwrap(),
                last,
                self.top.len() - 1,
            ))
        } else {
            Ok(SelBound {
                sys: self,
                index: sel.0.as_slice(),
            })
        }
    }

    /// Mutably binds detached selection index to make borrowed selection.
    /// `sel`  is not consumed.
    pub fn bind_mut<'a>(&'a mut self, sel: &'a Sel) -> SelBoundMut<'a> {
        // No need to check for empty index since it's guaranteed to be non-empty
        let last = unsafe { *sel.0.get_unchecked(sel.0.len() - 1) };
        if last >= self.top.len() {
            panic!("selection is out of bounds");
        } else {
            SelBoundMut {
                sys: self,
                index: sel.0.as_slice(),
            }
        }
    }

    /// Mutably binds detached selection index to make borrowed selection.
    /// `sel`  is not consumed.
    pub fn try_bind_mut<'a>(&'a mut self, sel: &'a Sel) -> Result<SelBoundMut<'a>, SelectionError> {
        // No need to check for empty index since it's guaranteed to be non-empty
        let last = unsafe { *sel.0.get_unchecked(sel.0.len() - 1) };
        if last >= self.top.len() {
            Err(SelectionError::IndexValidation(
                *sel.0.first().unwrap(),
                last,
                self.top.len() - 1,
            ))
        } else {
            Ok(SelBoundMut {
                sys: self,
                index: sel.0.as_slice(),
            })
        }
    }

    /// Returns mutable parallel iterator over parallel selections.
    pub fn iter_par_split_mut<'a>(&'a mut self, par: &'a ParSplit) -> impl IndexedParallelIterator<Item = SelParMut<'a>> {
        par.check_bounds(self);
        par.selections
            .par_iter()
            .map(|sel| SelParMut::new(self, &sel.0))
    }

    /// Returns parallel iterator over parallel selections.
    pub fn iter_par_split<'a>(&'a mut self, par: &'a ParSplit) -> impl IndexedParallelIterator<Item = SelPar<'a>> {
        par.check_bounds(self);
        par.selections
            .par_iter()
            .map(|sel| SelPar::new(self, &sel.0))
    }

    pub fn split_bound<'a, RT, F>(&'a self, func: F) -> impl Iterator<Item = SelOwnBound<'a>>
    where
        RT: Default + std::cmp::PartialEq + 'a,
        F: Fn(Particle) -> Option<RT> + 'a,
    {
        self.split(func).map(|sel| SelOwnBound {
            sys: self,
            index: sel.0,
        })
    }

    pub fn split_resindex_bound(&self) -> impl Iterator<Item = SelOwnBound<'_>> {
        self.split_bound(|p| Some(p.atom.resindex))
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
    pub fn append_from_self(&mut self, sel: &Sel) -> Result<Sel, SelectionError> {
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
        Ok(Sel::from_iter(old_last + 1..self.len())?)
    }

    pub fn append_atoms_coords<'a>(
        &mut self,
        atoms: impl AtomIterator<'a>,
        coords: impl PosIterator<'a>,
    ) -> Result<Sel, SelectionError> {
        let old_last = self.len() - 1;
        self.st.add_coords(coords.cloned());
        self.top.add_atoms(atoms.cloned());
        check_topology_state_sizes(&self.top, &self.st)?;
        Ok(Sel::from_iter(old_last + 1..self.len())?)
    }

    pub fn append(
        &mut self,
        data: &(impl AtomIterProvider + PosIterProvider),
    ) -> Result<Sel, SelectionError> {
        let old_last = self.len() - 1;
        self.st.add_coords(data.iter_pos().cloned());
        self.top.add_atoms(data.iter_atoms().cloned());
        check_topology_state_sizes(&self.top, &self.st)?;
        Ok(Sel::from_iter(old_last + 1..self.len())?)
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
        let all = self.select_all();
        for x in 0..=nbox[0] {
            for y in 0..=nbox[1] {
                for z in 0..=nbox[2] {
                    if x == 0 && y == 0 && z == 0 {
                        continue;
                    }
                    let added = self.append_from_self(&all)?;
                    let shift =
                        m.column(0) * x as f32 + m.column(1) * y as f32 + m.column(2) * z as f32;
                    self.select_bound_mut(added)?.translate(&shift);
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

impl Selectable for System {
    fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        Ok(Sel(def.into_sel_index(&self.top, &self.st, None)?))
    }
}

impl SaveTopology for System {}
impl SaveState for System {}
impl SaveTopologyState for System {}

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

// Analog of `bind` with `&sel >> &system` syntax
impl<'a> std::ops::Shr<&'a System> for &'a Sel {
    type Output = SelBound<'a>;
    fn shr(self, rhs: &'a System) -> Self::Output {
        rhs.bind(self)
    }
}

// Analog of `bind_mut` with `&sel >> &mut system` syntax
impl<'a> std::ops::Shr<&'a mut System> for &'a Sel {
    type Output = SelBoundMut<'a>;
    fn shr(self, rhs: &'a mut System) -> Self::Output {
        rhs.bind_mut(self)
    }
}
