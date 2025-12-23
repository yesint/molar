use crate::prelude::*;

//================================================
/// Read-only subsystem for non-blocking parallel access to atoms and posisitons
/// Doesn't have access to shared fields such as box and bonds.
//================================================
pub struct SelPar<'a> {
    index: &'a [usize],
    sys: *const System,
}

impl<'a> SelPar<'a> {
    pub(crate) fn new(sys: &'a System, index: &'a [usize]) -> SelPar<'a> {
        Self {
            index,
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
        unsafe{ (*self.sys).top.atoms.as_ptr()}
    }

    fn coords_ptr(&self) -> *const Pos {
        unsafe{ (*self.sys).st.coords.as_ptr()}
    }
}

impl NonAtomPosAnalysis for SelPar<'_> {
    fn top_ptr(&self) -> *const Topology {
        unsafe{ &raw const (*self.sys).top }
    }

    fn st_ptr(&self) -> *const State {
        unsafe{ &raw const (*self.sys).st }
    }
}

//================================================
/// Read-write subsystem for non-blocking parallel access to atoms and posisitons
/// Doesn't have access to shared fields such as box and bonds.
//================================================
pub struct SelParMut<'a> {
    index: &'a [usize],
    sys: *const System,
}

impl<'a> SelParMut<'a> {
    pub(crate) fn new(sys: &'a System, index: &'a [usize]) -> SelParMut<'a> {
        Self {
            index,
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
        unsafe{ (*self.sys).top.atoms.as_ptr()}
    }

    fn coords_ptr(&self) -> *const Pos {
        unsafe{ (*self.sys).st.coords.as_ptr()}
    }
}

impl AtomPosAnalysisMut for SelParMut<'_> {
    fn atoms_ptr_mut(&mut self) -> *mut Atom {
        // This creates temp &sys for each parallel selection
        // since no &mut sys is ever created this is fine
        unsafe{ (*self.sys).top.atoms.as_ptr() as *mut Atom }
    }

    fn coords_ptr_mut(&mut self) -> *mut Pos {
        // This creates temp &sys for each parallel selection
        // since no &mut sys is ever created this is fine
        unsafe{ (*self.sys).st.coords.as_ptr() as *mut Pos }
    }
}

impl NonAtomPosAnalysis for SelParMut<'_> {
    fn top_ptr(&self) -> *const Topology {
        unsafe{ &raw const (*self.sys).top }
    }

    fn st_ptr(&self) -> *const State {
        unsafe{ &raw const (*self.sys).st }
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

impl ParSplit {
    pub(crate) fn check_bounds(&self, sys: &System) {
        if self.max_index >= sys.len() {
            panic!("max index of ParSplit is out of bounds");
        }
    }

    pub fn get_bound_mut<'a>(&'a self, sys: &'a mut System, i: usize) -> SelParMut<'a> {
        self.check_bounds(sys);
        SelParMut::new(sys, &self.selections[i].0)
    }

    pub fn get_bound<'a>(&'a self, sys: &'a System, i: usize) -> SelPar<'a> {
        self.check_bounds(sys);
        SelPar::new(sys, &self.selections[i].0)
    }

    pub fn into_selections(self) -> Vec<Sel> {
        self.selections
    }
}


#[cfg(test)]
mod tests {
    use rayon::iter::ParallelIterator;

    use crate::prelude::*;

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
        sys.iter_par_split_mut(&par) 
            // Run unwrap on each selection in parallel
            .try_for_each(|mut sel| sel.unwrap_simple())?; 
        Ok(())
    }
}
