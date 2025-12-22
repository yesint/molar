use crate::prelude::*;

/// Trait for objects that support selecting
pub trait Selectable {
    fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError>;
}

//============================================================
/// Umbrella trait for implementing read-only analysis traits
/// involving only atoms and positions.
///
/// It assumes that atoms and coordinates are stored in 
/// contigous arrays.
//============================================================
pub trait AtomPosAnalysis: LenProvider + IndexProvider + Sized {
    // Using raw pointers here allows to easily implement fast unsfe accessors
    // without working around borrow checker issues.

    // Raw pointer to the beginning of array of atoms
    fn atoms_ptr(&self) -> *const Atom;
    // Raw pointer to the beginning of array of coords
    fn coords_ptr(&self) -> *const Pos;

    /// Creates a parallel split based on provided closure.
    /// A closure takes a [Particle] and returns a distinct value for each piece.
    /// New selection is created whenever new return value differes from the previous one.
    fn split_par<F, R>(&self, func: F) -> Result<ParSplit, SelectionError>
    where
        F: Fn(Particle) -> Option<R>,
        R: Default + PartialOrd,
        Self: Sized,
    {
        let selections: Vec<Sel> = self.split(func).collect();

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
    fn split<RT, F>(&self, func: F) -> impl Iterator<Item = Sel>
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

    fn split_resindex(&self) -> impl Iterator<Item = Sel> {
        self.split(|p| Some(p.atom.resindex))
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
                self.coords_ptr().add(ind) as *mut f32
            },
            |i: usize| self.get_particle(i).unwrap().atom.vdw(),
        )
    }
}

//============================================================
/// Umbrella trait for implementing read-only analysis traits
///  NOT involving atoms and positions
//============================================================
pub trait NonAtomPosAnalysis: LenProvider + IndexProvider + Sized {
    fn top_ptr(&self) -> *const Topology;
    fn st_ptr(&self) -> *const State;

    /// Splits by molecule and returns an iterator over them.
    /// If molecule is only partially contained in self then only this part is returned (molecules are clipped).
    /// If there are no molecules in [Topology] return an empty iterator.
    fn split_mol_iter(&self) -> impl Iterator<Item = Sel>
    where
        Self: Sized,
    {
        // Iterate over molecules and find those inside selection
        let first = self.first_index();
        let last =  self.last_index();

        let mut molid = 0;

        let next_fn = move || {
            if unsafe{&*self.top_ptr()}.num_molecules() == 0 {
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
            }
            .map(|r| Sel(unsafe { SVec::from_sorted(r.into_iter().collect()) }));

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
pub trait AtomPosAnalysisMut: LenProvider + IndexProvider + Sized {
    // Raw pointer to the beginning of array of atoms
    fn atoms_ptr_mut(&mut self) -> *mut Atom;
    // Raw pointer to the beginning of array of coords
    fn coords_ptr_mut(&mut self) -> *mut Pos;
}

//================================================
/// Umbrella trait for implementing read-write analysis traits NOT involving atoms and positions
//================================================
pub trait NonAtomPosAnalysisMut: Sized {
    fn top_ptr_mut(&mut self) -> *mut Topology;
    fn st_ptr_mut(&mut self) -> *mut State;
}

//═══════════════════════════════════════════════════════════
//  Blanket trait implementations
//═══════════════════════════════════════════════════════════

//██████  Immutable analysis traits involving atoms and positions only

impl<T: AtomPosAnalysis> PosIterProvider for T {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        unsafe { self.iter_index().map(|i| &*self.coords_ptr().add(i)) }
    }
}

impl<T: AtomPosAnalysis> AtomIterProvider for T {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        unsafe { self.iter_index().map(|i| &*self.atoms_ptr().add(i)) }
    }
}

impl<T: AtomPosAnalysis> RandomPosProvider for T {
    unsafe fn get_pos_unchecked(&self, i: usize) -> &Pos {
        let ind = self.get_index_unchecked(i);
        &*self.coords_ptr().add(ind)
    }
}

impl<T: AtomPosAnalysis> RandomAtomProvider for T {
    unsafe fn get_atom_unchecked(&self, i: usize) -> &Atom {
        let ind = self.get_index_unchecked(i);
        &*self.atoms_ptr().add(ind)
    }
}

impl<T: AtomPosAnalysis> RandomParticleProvider for T {
    unsafe fn get_particle_unchecked(&self, i: usize) -> Particle<'_> {
        let ind = self.get_index_unchecked(i);
        Particle {
            id: ind,
            atom: &*self.atoms_ptr().add(ind),
            pos: &*self.coords_ptr().add(ind),
        }
    }
}

//██████  Immutable analysis traits NOT involving atoms and positions

impl<T: NonAtomPosAnalysis> BoxProvider for T {
    fn get_box(&self) -> Option<&PeriodicBox> {
        unsafe{&*self.st_ptr()}.get_box()
    }
}

impl<T: NonAtomPosAnalysis> TimeProvider for T {
    fn get_time(&self) -> f32 {
        unsafe{&*self.st_ptr()}.time
    }
}

impl<T: NonAtomPosAnalysis> RandomMoleculeProvider for T {
    fn num_molecules(&self) -> usize {
        unsafe{&*self.top_ptr()}.num_molecules()
    }

    unsafe fn get_molecule_unchecked(&self, i: usize) -> &[usize; 2] {
        unsafe{&*self.top_ptr()}.get_molecule_unchecked(i)
    }
}

impl<T: NonAtomPosAnalysis> MoleculeIterProvider for T {
    fn iter_molecules(&self) -> impl Iterator<Item = &[usize; 2]> {
        unsafe{&*self.top_ptr()}.iter_molecules()
    }
}

impl<T: NonAtomPosAnalysis> RandomBondProvider for T {
    fn num_bonds(&self) -> usize {
        unsafe{&*self.top_ptr()}.num_bonds()
    }

    unsafe fn get_bond_unchecked(&self, i: usize) -> &[usize; 2] {
        unsafe{&*self.top_ptr()}.get_bond_unchecked(i)
    }
}

impl<T: NonAtomPosAnalysis> BondIterProvider for T {
    fn iter_bonds(&self) -> impl Iterator<Item = &[usize; 2]> {
        unsafe{&*self.top_ptr()}.iter_bonds()
    }
}

//██████  Mutable analysis traits involving only atoms and positions

impl<T: AtomPosAnalysisMut> RandomAtomMutProvider for T {
    unsafe fn get_atom_mut_unchecked(&mut self, i: usize) -> &mut Atom {
        let ind = self.get_index_unchecked(i);
        &mut *self.atoms_ptr_mut().add(ind)
    }
}

impl<T: AtomPosAnalysisMut> AtomIterMutProvider for T {
    fn iter_atoms_mut(&mut self) -> impl AtomMutIterator<'_> {
        let mut i = 0;
        let func = move || {
            if i < self.len() {
                let ind = unsafe { self.get_index_unchecked(i) };
                let res = unsafe { &mut *self.atoms_ptr_mut().add(ind) };
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
                let res = unsafe { &mut *self.coords_ptr_mut().add(ind) };
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
        &mut *self.coords_ptr_mut().add(ind)
    }
}

impl<T: AtomPosAnalysisMut> RandomParticleMutProvider for T {
    unsafe fn get_particle_mut_unchecked(&mut self, i: usize) -> ParticleMut<'_> {
        let ind = self.get_index_unchecked(i);
        ParticleMut {
            id: ind,
            atom: unsafe { &mut *self.atoms_ptr_mut().add(ind) },
            pos: unsafe { &mut *self.coords_ptr_mut().add(ind) },
        }
    }
}

//██████  Mutable analysis traits NOT involving atoms and positions

impl<T: NonAtomPosAnalysisMut> TimeMutProvider for T {
    fn set_time(&mut self, t: f32) {
        unsafe{&mut *self.st_ptr_mut()}.time = t;
    }
}

impl<T: NonAtomPosAnalysisMut> BoxMutProvider for T {
    fn get_box_mut(&mut self) -> Option<&mut PeriodicBox> {
        unsafe{&mut *self.st_ptr_mut()}.pbox.as_mut()
    }
}

//██████  Measure traits

impl<T: AtomPosAnalysis> MeasurePos for T {}
impl<T: AtomPosAnalysis + NonAtomPosAnalysis> MeasurePeriodic for T {}
impl<T: AtomPosAnalysis> MeasureMasses for T {}
impl<T: AtomPosAnalysis> MeasureRandomAccess for T {}

//██████  Modify traits

impl<T: AtomPosAnalysisMut> ModifyPos for T {}
impl<T: AtomPosAnalysisMut + NonAtomPosAnalysis> ModifyPeriodic for T {}
impl<T: AtomPosAnalysisMut + AtomPosAnalysis + NonAtomPosAnalysis> ModifyRandomAccess for T {}
