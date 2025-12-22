use crate::prelude::*;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator};

//================================================
/// Read-only subsystem for non-blocking parallel access to atoms and posisitons
/// Doesn't have access to shared fields such as box and bonds.
//================================================
pub struct SelPar<'a> {
    coords_ptr: *const Pos,
    atoms_ptr: *const Atom,
    index: &'a [usize],
    sys: &'a System,
}

impl<'a> SelPar<'a> {
    pub(crate) fn new(sys: &'a System, index: &'a [usize]) -> SelPar<'a> {
        Self {
            index,
            atoms_ptr: sys.atoms_ptr(),
            coords_ptr: sys.coords_ptr(),
            sys,
        }
    }
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

impl NonAtomPosAnalysis for SelPar<'_> {
    fn top_ptr(&self) -> *const Topology {
        &self.sys.top
    }

    fn st_ptr(&self) -> *const State {
        &self.sys.st
    }
}

//================================================
/// Read-write subsystem for non-blocking parallel access to atoms and posisitons
/// Doesn't have access to shared fields such as box and bonds.
//================================================
pub struct SelParMut<'a> {
    coords_ptr: *mut Pos,
    atoms_ptr: *mut Atom,
    index: &'a [usize],
    sys: &'a System,
}

impl<'a> SelParMut<'a> {
    pub(crate) fn new(sys: &'a System, index: &'a [usize]) -> SelParMut<'a> {
        Self {
            index,
            atoms_ptr: sys.atoms_ptr() as *mut Atom,
            coords_ptr: sys.coords_ptr() as *mut Pos,
            sys,
        }
    }
}

unsafe impl Sync for SelParMut<'_> {}
unsafe impl Send for SelParMut<'_> {}

impl LenProvider for SelParMut<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SelParMut<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl AtomPosAnalysis for SelParMut<'_> {
    fn atoms_ptr(&self) -> *const Atom {
        self.atoms_ptr
    }

    fn coords_ptr(&self) -> *const Pos {
        self.coords_ptr
    }
}

impl AtomPosAnalysisMut for SelParMut<'_> {
    fn atoms_ptr_mut(&mut self) -> *mut Atom {
        self.atoms_ptr
    }

    fn coords_ptr_mut(&mut self) -> *mut Pos {
        self.coords_ptr
    }
}

impl NonAtomPosAnalysis for SelParMut<'_> {
    fn top_ptr(&self) -> *const Topology {
        &self.sys.top
    }

    fn st_ptr(&self) -> *const State {
        &self.sys.st
    }
}
//============================================================================
/// Collection of non-overlapping selections that could be mutated in parallel
/// Selections don't have access to shared fields such as box and bonds.
//============================================================================
pub struct ParSplit {
    pub(crate) selections: Vec<Sel>,
    pub(crate) max_index: usize,
}

/// Bound immutable parallel split
pub struct ParSplitBound<'a> {
    pub(super) sys: &'a System,
    pub(super) indexes: &'a Vec<Sel>,
}

impl<'a> ParSplitBound<'a> {
    /// Returns parallel iterator over parallel selections.
    pub fn par_iter(&'a self) -> impl IndexedParallelIterator<Item = SelPar<'a>> {
        self.indexes
            .par_iter()
            .map(|ind| SelPar::new(self.sys, &ind.0))
    }

    /// Returns serial iterator over parallel selections.
    pub fn iter(&'a self) -> impl Iterator<Item = SelPar<'a>> {
        self.indexes.iter().map(|ind| SelPar::new(self.sys, &ind.0))
    }

    pub fn get(&self, i: usize) -> SelPar<'_> {
        SelPar::new(self.sys, &self.indexes[i].0)
    }
}

/// Bound mutable parallel split
pub struct ParSplitBoundMut<'a> {
    pub(super) sys: &'a mut System,
    pub(super) indexes: &'a Vec<Sel>,
}

impl<'a> ParSplitBoundMut<'a> {
    /// Returns parallel iterator over parallel selections.
    pub fn par_iter_mut(&'a mut self) -> impl IndexedParallelIterator<Item = SelParMut<'a>> {
        self.indexes
            .par_iter()
            .map(|ind| SelParMut::new(self.sys, &ind.0))
    }

    /// Returns serial iterator over parallel selections.
    pub fn iter_mut(&'a mut self) -> impl Iterator<Item = SelParMut<'a>> {
        self.indexes
            .iter()
            .map(|ind| SelParMut::new(self.sys, &ind.0))
    }

    pub fn get_mut(&mut self, i: usize) -> SelParMut<'_> {
        SelParMut::new(self.sys, &self.indexes[i].0)
    }
}

#[cfg(test)]
mod tests {
    use rayon::iter::ParallelIterator;

    use crate::core::modify::ModifyPeriodic;
    use crate::core::{AtomPosAnalysis, System};

    #[test]
    fn par_unwrap() -> anyhow::Result<()> {
        // Load file
        let mut sys = System::from_file("tests/membr.gro")?;

        // Make a parallel split for each POPG lipid molecule
        let par = sys.split_par(|p| {
            if p.atom.resname == "POPG" {
                // Whenever new distinct result is returned form this closure
                // new selection is created, so each distinct POPG residue
                // becomes a separate selection.
                Some(p.atom.resindex)
            } else {
                // All other atoms are ignored
                None
            }
        })?;

        // Bind split to a system
        sys.bind_par_mut(&par)? 
            // Get rayon parallel iterator over selections
            .par_iter_mut() 
            // Run unwrap on each selection in parallel
            .try_for_each(|mut sel| sel.unwrap_simple())?; 
        Ok(())
    }
}
