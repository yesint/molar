use rayon::iter::IndexedParallelIterator;
use crate::{prelude::*, selection::utils::{difference_sorted, intersection_sorted, union_sorted}};

/// Trait for objects that support creating selections
pub trait Selectable {
    // Make unbound immutable selection
    fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError>;
}

/// Trait for objects that support creating bound selections
pub trait SelectableBound: SysProvider + Selectable {
    // Make bound immutable selection
    fn select_bound(&self, def: impl SelectionDef) -> Result<SelOwnBound<'_>, SelectionError>;
}

//=============================================================================
/// Trait for things containing reference to System (mostly selections)
/// (previously called SystemProvider)
pub trait SysProvider {
    fn get_system_ptr(&self) -> *const System;
}

/// Trait for things containing mut reference to System (mostly selections)
/// (previously called SystemMutProvider)
pub trait SysMutProvider: SysProvider {
    fn get_system_mut(&mut self) -> *mut System {
        self.get_system_ptr() as *mut System
    }
}

/// Backward-compatible alias for [`SysProvider`].
pub use SysProvider as SystemProvider;
/// Backward-compatible alias for [`SysMutProvider`].
pub use SysMutProvider as SystemMutProvider;

//=============================================================================
/// Blanket impls: anything that has a SysProvider + IndexProvider gets all
/// the element traits for free.
//=============================================================================

impl<T: SysProvider + IndexProvider> PosProvider for T {
    unsafe fn coords_ptr(&self) -> *const Pos {
        unsafe { (*self.get_system_ptr()).st.coords.as_ptr() }
    }
}

impl<T: SysMutProvider + IndexProvider> PosMutProvider for T {}

impl<T: SysProvider + IndexProvider> AtomProvider for T {
    unsafe fn atoms_ptr(&self) -> *const Atom {
        unsafe { (*self.get_system_ptr()).top.atoms.as_ptr() }
    }
}

impl<T: SysMutProvider + IndexProvider> AtomMutProvider for T {}

impl<T: SysProvider + IndexProvider> BoxProvider for T {
    fn get_box(&self) -> Option<&PeriodicBox> {
        unsafe { (*self.get_system_ptr()).st.get_box() }
    }
}

impl<T: SysMutProvider + IndexProvider> BoxMutProvider for T {
    fn get_box_mut(&mut self) -> Option<&mut PeriodicBox> {
        unsafe { (*self.get_system_mut()).st.pbox.as_mut() }
    }
}

impl<T: SysProvider + IndexProvider> TimeProvider for T {
    fn get_time(&self) -> f32 {
        unsafe { (*self.get_system_ptr()).st.time }
    }
}

impl<T: SysMutProvider + IndexProvider> TimeMutProvider for T {
    fn set_time(&mut self, t: f32) {
        unsafe { (*self.get_system_mut()).st.time = t; }
    }
}

impl<T: SysProvider + IndexProvider> BondProvider for T {
    fn num_bonds(&self) -> usize {
        unsafe { BondProvider::num_bonds(&(*self.get_system_ptr()).top) }
    }

    unsafe fn get_bond_unchecked(&self, i: usize) -> &[usize; 2] {
        unsafe { (*self.get_system_ptr()).top.get_bond_unchecked(i) }
    }

    fn iter_bonds(&self) -> impl Iterator<Item = &[usize; 2]> {
        unsafe { (*self.get_system_ptr()).top.iter_bonds() }
    }
}

impl<T: SysProvider + IndexProvider> MolProvider for T {
    fn num_molecules(&self) -> usize {
        unsafe { (*self.get_system_ptr()).top.num_molecules() }
    }

    unsafe fn get_molecule_unchecked(&self, i: usize) -> &[usize; 2] {
        unsafe { (*self.get_system_ptr()).top.get_molecule_unchecked(i) }
    }

    fn iter_molecules(&self) -> impl Iterator<Item = &[usize; 2]> {
        unsafe { (*self.get_system_ptr()).top.iter_molecules() }
    }
}

//████████  Particle iteration
// Blanket for anything with PosProvider + AtomProvider + IndexProvider

impl<T: PosProvider + AtomProvider + IndexProvider> ParticleIterProvider for T {
    fn iter_particle(&self) -> impl Iterator<Item = Particle<'_>> {
        let cp = unsafe { self.coords_ptr() };
        let ap = unsafe { self.atoms_ptr() };
        unsafe {
            self.iter_index().map(move |i| Particle {
                id: i,
                atom: &*ap.add(i),
                pos: &*cp.add(i),
            })
        }
    }
}

impl<T: PosProvider + AtomProvider + IndexProvider + IndexParProvider> ParticleParIterProvider for T {
    fn par_iter_particle(&self) -> impl IndexedParallelIterator<Item = Particle<'_>> {
        let cp = unsafe { self.coords_ptr() } as usize;
        let ap = unsafe { self.atoms_ptr() } as usize;
        unsafe {
            self.par_iter_index().map(move |i| Particle {
                id: i,
                atom: &*(ap as *const Atom).add(i),
                pos: &*(cp as *const Pos).add(i),
            })
        }
    }
}

impl<T: PosMutProvider + AtomMutProvider + IndexProvider> ParticleIterMutProvider for T {
    fn iter_particle_mut(&mut self) -> impl Iterator<Item = ParticleMut<'_>> {
        let cp = unsafe { self.coords_ptr_mut() };
        let ap = unsafe { self.atoms_ptr_mut() };
        unsafe {
            self.iter_index().map(move |i| ParticleMut {
                id: i,
                atom: &mut *ap.add(i),
                pos: &mut *cp.add(i),
            })
        }
    }
}

//████████  Measure blanket impls

impl<T: PosProvider + LenProvider> MeasurePos for T {}
impl<T: PosProvider + AtomProvider + BoxProvider + LenProvider> MeasurePeriodic for T {}
impl<T: PosProvider + AtomProvider + LenProvider> MeasureMasses for T {}
impl<T: PosProvider> MeasureRandomAccess for T {}
impl<T: PosProvider + AtomProvider> MeasureAtomPos for T {}

//████████  Modify blanket impls

impl<T: PosMutProvider> ModifyPos for T {}
impl<T: PosMutProvider + BoxProvider + LenProvider> ModifyPeriodic for T {}
impl<T: PosMutProvider + BoxProvider + Sized> ModifyRandomAccess for T {}
impl<T: AtomMutProvider + LenProvider> ModifyAtoms for T {}

//=============================================================================
/// Trait for logical operations in any (borrowed or own) bound selections
pub trait SelectionLogic: IndexSliceProvider {
    type DerivedSel;
    fn clone_with_index(&self, index: SVec) -> Self::DerivedSel;

    fn or(&self, rhs: &impl IndexSliceProvider) -> Self::DerivedSel {
        let index = unsafe {union_sorted(self.get_index_slice(), rhs.get_index_slice())};
        self.clone_with_index(index)
    }

    fn and(&self, rhs: &impl IndexSliceProvider) -> Result<Self::DerivedSel,SelectionError> {
        let index = unsafe {intersection_sorted(self.get_index_slice(), rhs.get_index_slice())};
        if index.is_empty() {
            return Err(SelectionError::EmptyIntersection)
        }
        Ok(self.clone_with_index(index))
    }

    fn minus(&self, rhs: &impl IndexSliceProvider) -> Result<Self::DerivedSel,SelectionError> {
        let index = unsafe {difference_sorted(self.get_index_slice(), rhs.get_index_slice())};
        if index.is_empty() {
            return Err(SelectionError::EmptyDifference)
        }
        Ok(self.clone_with_index(index))
    }

    fn invert(&self, rhs: &impl IndexSliceProvider) -> Result<Self::DerivedSel,SelectionError>
    where Self: SysProvider,
    {
        let all = (0..unsafe{&*self.get_system_ptr()}.len()).into_iter().collect::<Vec<_>>();
        let index = unsafe {difference_sorted(&all, rhs.get_index_slice())};
        if index.is_empty() {
            return Err(SelectionError::EmptyDifference)
        }
        Ok(self.clone_with_index(index))
    }
}

//=============================================================================
// split / split_par  (previously on AtomPosAnalysis)
//=============================================================================

/// Trait providing `split` and `split_par` for selections that have full particle access.
///
/// Anything implementing `PosProvider + AtomProvider + IndexProvider` gets
/// `split`, `split_par`, `split_resindex`, `whole_attr`, `whole_residues`, `whole_chains`.
pub trait Analysis: PosProvider + AtomProvider + IndexProvider + Sized {
    /// Creates a parallel split based on provided closure.
    fn split_par<F, R>(&self, func: F) -> Result<ParSplit, SelectionError>
    where
        F: Fn(Particle) -> Option<R>,
        R: Default + PartialOrd,
    {
        let selections: Vec<Sel> = self.split(func).collect();

        if selections.is_empty() {
            return Err(SelectionError::EmptySplit);
        }

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
                let global_i = unsafe { self.get_index_unchecked(cur) };
                let p = Particle {
                    id: global_i,
                    atom: unsafe { self.get_atom_unchecked(cur) },
                    pos: unsafe { self.get_pos_unchecked(cur) },
                };
                let id = p.id;
                let val = func(p);

                if let Some(val) = val {
                    if val == cur_val {
                        index.push(id);
                    } else if index.is_empty() {
                        cur_val = val;
                        index.push(id);
                    } else {
                        cur_val = val;
                        return Some(Sel(unsafe { SVec::from_sorted(index) }));
                    }
                }
                cur += 1;
            }

            if !index.is_empty() {
                return Some(Sel(unsafe { SVec::from_sorted(index) }));
            }

            None
        };

        std::iter::from_fn(next_fn)
    }

    fn split_resindex(&self) -> impl Iterator<Item = Sel> {
        self.split(|p| Some(p.atom.resindex))
    }

    /// Creates an "expanded" selection that includes all atoms with the same attributes.
    fn whole_attr<T>(&self, attr_fn: fn(&Atom) -> &T) -> Sel
    where
        T: Eq + std::hash::Hash + Copy,
    {
        let mut properties = std::collections::HashSet::<T>::new();
        for at in self.iter_atoms() {
            properties.insert(*attr_fn(at));
        }

        let mut ind = vec![];
        for (i, at) in self.iter_atoms().enumerate() {
            let cur_prop = attr_fn(at);
            if properties.contains(cur_prop) {
                ind.push(i);
            }
        }

        Sel(unsafe { SVec::from_sorted(ind) })
    }

    /// Selects whole residues present in the current selection (in terms of resindex)
    fn whole_residues(&self) -> Sel {
        self.whole_attr(|at| &at.resindex)
    }

    /// Selects whole chains present in the current selection
    fn whole_chains(&self) -> Sel {
        self.whole_attr(|at| &at.chain)
    }
}

/// Blanket impl: anything with PosProvider + AtomProvider + IndexProvider gets Analysis
impl<T: PosProvider + AtomProvider + IndexProvider> Analysis for T {}
