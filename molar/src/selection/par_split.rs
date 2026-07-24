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

impl IndexSliceProvider for SelPar<'_> {
    fn get_index_slice(&self) -> &[usize] {
        self.index
    }
}

impl SystemProvider for SelPar<'_> {
    fn get_system_ptr(&self) -> *const System {
        self.sys
    }
}

//================================================
/// Read-write subsystem for non-blocking parallel access to atoms and posisitons
/// Doesn't have access to shared fields such as box and bonds.
//================================================
pub struct SelParMut<'a> {
    index: &'a [usize],
    sys: *mut System,
}

impl<'a> SelParMut<'a> {
    pub(crate) fn new(sys: *mut System, index: &'a [usize]) -> SelParMut<'a> {
        Self {
            index,
            sys,
        }
    }
}

unsafe impl Sync for SelParMut<'_> {}
unsafe impl Send for SelParMut<'_> {}

impl IndexSliceProvider for SelParMut<'_> {
    fn get_index_slice(&self) -> &[usize] {
        self.index
    }
}

impl SystemProvider for SelParMut<'_> {
    fn get_system_ptr(&self) -> *const System {
        self.sys as *const System
    }
}

// PosMutProvider and AtomMutProvider are implemented directly (not via SysMutProvider)
// so that shared fields (box, bonds, molecules) cannot be mutated during parallel use.
impl PosMutProvider for SelParMut<'_> {
    unsafe fn coords_ptr_mut(&mut self) -> *mut Pos { unsafe {
        (*self.sys).st.coords.as_mut_ptr()
    }}
}
impl AtomMutProvider for SelParMut<'_> {
    fn atom_storage_mut(&mut self) -> &mut AtomStorage {
        unsafe { &mut (*self.sys).top.atoms }
    }
}
// VelMutProvider and ForceMutProvider are implemented directly (same reason as above).
// VelProvider and ForceProvider are covered by blanket impls via SystemProvider.
impl VelMutProvider for SelParMut<'_> {
    unsafe fn vel_ptr_mut(&mut self) -> *mut Vel { unsafe {
        let v = &mut (*self.sys).st.velocities;
        if v.is_empty() { std::ptr::null_mut() } else { v.as_mut_ptr() }
    }}
}
impl ForceMutProvider for SelParMut<'_> {
    unsafe fn force_ptr_mut(&mut self) -> *mut Force { unsafe {
        let v = &mut (*self.sys).st.forces;
        if v.is_empty() { std::ptr::null_mut() } else { v.as_mut_ptr() }
    }}
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
        SelParMut::new(sys as *mut System, &self.selections[i].0)
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
            if p.atom.get_resname() == "POPG" {
                // Whenever new distinct result is returned form this closure
                // new selection is created, so each distinct POPG residue
                // becomes a separate selection.
                Some(p.atom.get_resindex())
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

    /// Mutate an *atom* column (not just positions) in parallel over disjoint index sets and
    /// verify the writes land correctly — exercises `AtomRefMut` core-column setters across
    /// threads (the `SelParMut` → `iter_atoms_mut` path under SoA storage).
    #[test]
    fn par_atom_column_write() -> anyhow::Result<()> {
        // In-memory system: 20 residues × 5 atoms.
        let mut top = Topology::default();
        for res in 0..20i32 {
            for _ in 0..5 {
                top.atoms
                    .push_row(&Atom::new().with_name("C").with_resname("RES").with_resid(res));
            }
        }
        top.assign_resindex();
        let n = top.atoms.len();
        let mut sys = System::new(top, State::new_fake(n))?;

        // One parallel selection per residue (disjoint index sets).
        let par = sys.split_par(|p| Some(p.atom.get_resindex()))?;

        // In parallel, set each atom's bfactor to its resid.
        sys.iter_par_split_mut(&par).try_for_each(|mut sel| {
            for mut a in sel.iter_atoms_mut() {
                let r = a.get_resid() as Float;
                a.set_bfactor(r);
            }
            Ok::<_, SelectionError>(())
        })?;

        // Every atom's bfactor must now equal its resid.
        for a in sys.iter_atoms() {
            assert_eq!(a.get_bfactor(), a.get_resid() as Float);
        }
        Ok(())
    }
}
