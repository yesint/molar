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
/// Internal trait for things containing reference to System (mostly selections)
/// This is pub(crate) so that it's not part of public API.
pub(crate) trait SysProvider: LenProvider + IndexProvider {
    fn sys_ptr(&self) -> *const System;
}

/// Internal trait for things containing mutable reference to System
pub(crate) trait SysMutProvider: SysProvider {
    fn sys_ptr_mut(&self) -> *mut System;
}

//=============================================================================
/// Public trait for things containing reference to System (kept for backwards compat)
pub trait SystemProvider {
    fn get_system_ptr(&self) -> *const System;
}

/// Public trait for things containing mutable reference to System (kept for backwards compat)
pub trait SystemMutProvider: SystemProvider {
    fn get_system_mut(&mut self) -> *mut System {
        self.get_system_ptr() as *mut System
    }
}

// SysProvider <-> SystemProvider bridge
impl<T: SysProvider> SystemProvider for T {
    fn get_system_ptr(&self) -> *const System {
        self.sys_ptr()
    }
}

//=============================================================================
// Blanket impls: SysProvider → all read-only element traits

impl<T: SysProvider> PosProvider for T {
    unsafe fn pos_unchecked(&self, i: usize) -> &Pos {
        let sys = &*self.sys_ptr();
        &*sys.st.coords.as_ptr().add(i)
    }
}

impl<T: SysProvider> AtomProvider for T {
    unsafe fn atom_unchecked(&self, i: usize) -> &Atom {
        let sys = &*self.sys_ptr();
        &*sys.top.atoms.as_ptr().add(i)
    }
}

impl<T: SysProvider> BoxProvider for T {
    fn get_box(&self) -> Option<&PeriodicBox> {
        unsafe { (*self.sys_ptr()).st.pbox.as_ref() }
    }
}

impl<T: SysProvider> TimeProvider for T {
    fn get_time(&self) -> f32 {
        unsafe { (*self.sys_ptr()).st.time }
    }
}

impl<T: SysProvider> BondProvider for T {
    fn bonds_raw(&self) -> &[[usize; 2]] {
        unsafe { &(*self.sys_ptr()).top.bonds }
    }
}

impl<T: SysProvider> MolProvider for T {
    fn molecules_raw(&self) -> &[[usize; 2]] {
        unsafe { &(*self.sys_ptr()).top.molecules }
    }
}

//=============================================================================
// Blanket impls: SysMutProvider → all mutable element traits

impl<T: SysMutProvider> PosMutProvider for T {
    unsafe fn pos_mut_unchecked(&mut self, i: usize) -> &mut Pos {
        let sys = &mut *self.sys_ptr_mut();
        &mut *sys.st.coords.as_mut_ptr().add(i)
    }
}

impl<T: SysMutProvider> AtomMutProvider for T {
    unsafe fn atom_mut_unchecked(&mut self, i: usize) -> &mut Atom {
        let sys = &mut *self.sys_ptr_mut();
        &mut *sys.top.atoms.as_mut_ptr().add(i)
    }
}

impl<T: SysMutProvider> BoxMutProvider for T {
    fn get_box_mut(&mut self) -> Option<&mut PeriodicBox> {
        unsafe { (*self.sys_ptr_mut()).st.pbox.as_mut() }
    }
}

impl<T: SysMutProvider> TimeMutProvider for T {
    fn set_time(&mut self, t: f32) {
        unsafe { (*self.sys_ptr_mut()).st.time = t; }
    }
}

//=============================================================================
// Particle iteration blanket impls (for anything with SysProvider + IndexProvider)

impl<T: SysProvider + IndexProvider> ParticleIterProvider for T {
    fn iter_particle(&self) -> impl Iterator<Item = Particle<'_>> {
        let ap = unsafe { (*self.sys_ptr()).top.atoms.as_ptr() } as usize;
        let cp = unsafe { (*self.sys_ptr()).st.coords.as_ptr() } as usize;
        self.iter_index().map(move |i| unsafe {
            Particle {
                id: i,
                atom: &*(ap as *const Atom).add(i),
                pos: &*(cp as *const Pos).add(i),
            }
        })
    }
}

impl<T: SysProvider + IndexParProvider> ParticleParIterProvider for T {
    fn par_iter_particle(&self) -> impl IndexedParallelIterator<Item = Particle<'_>> {
        let ap = unsafe { (*self.sys_ptr()).top.atoms.as_ptr() } as usize;
        let cp = unsafe { (*self.sys_ptr()).st.coords.as_ptr() } as usize;
        self.par_iter_index().map(move |i| unsafe {
            Particle {
                id: i,
                atom: &*(ap as *const Atom).add(i),
                pos: &*(cp as *const Pos).add(i),
            }
        })
    }
}

impl<T: SysMutProvider + IndexProvider> ParticleIterMutProvider for T {
    fn iter_particle_mut(&mut self) -> impl Iterator<Item = ParticleMut<'_>> {
        let ap = unsafe { (*self.sys_ptr_mut()).top.atoms.as_mut_ptr() } as usize;
        let cp = unsafe { (*self.sys_ptr_mut()).st.coords.as_mut_ptr() } as usize;
        self.iter_index().map(move |i| unsafe {
            ParticleMut {
                id: i,
                atom: &mut *(ap as *mut Atom).add(i),
                pos: &mut *(cp as *mut Pos).add(i),
            }
        })
    }
}

impl<T: SysMutProvider + IndexParProvider> ParticleParIterMutProvider for T {
    fn par_iter_particle_mut(&mut self) -> impl IndexedParallelIterator<Item = ParticleMut<'_>> {
        let ap = unsafe { (*self.sys_ptr_mut()).top.atoms.as_mut_ptr() } as usize;
        let cp = unsafe { (*self.sys_ptr_mut()).st.coords.as_mut_ptr() } as usize;
        self.par_iter_index().map(move |i| unsafe {
            ParticleMut {
                id: i,
                atom: &mut *(ap as *mut Atom).add(i),
                pos: &mut *(cp as *mut Pos).add(i),
            }
        })
    }
}

impl<T: SysMutProvider + IndexParProvider> PosParIterMutProvider for T {
    fn par_iter_pos_mut(&mut self) -> impl IndexedParallelIterator<Item = &mut Pos> {
        let p = unsafe { (*self.sys_ptr_mut()).st.coords.as_mut_ptr() } as usize;
        self.par_iter_index().map(move |i| unsafe {
            &mut *(p as *mut Pos).add(i)
        })
    }
}

impl<T: SysMutProvider + IndexParProvider> AtomParIterMutProvider for T {
    fn par_iter_atoms_mut(&mut self) -> impl IndexedParallelIterator<Item = &mut Atom> {
        let p = unsafe { (*self.sys_ptr_mut()).top.atoms.as_mut_ptr() } as usize;
        self.par_iter_index().map(move |i| unsafe {
            &mut *(p as *mut Atom).add(i)
        })
    }
}

// RandomParticleProvider for SysProvider types
impl<T: SysProvider + IndexProvider> RandomParticleProvider for T {
    unsafe fn get_particle_unchecked(&self, i: usize) -> Particle<'_> {
        let ind = self.get_index_unchecked(i);
        Particle {
            id: ind,
            atom: &*(*self.sys_ptr()).top.atoms.as_ptr().add(ind),
            pos: &*(*self.sys_ptr()).st.coords.as_ptr().add(ind),
        }
    }
}

// RandomParticleMutProvider for SysMutProvider types
impl<T: SysMutProvider + IndexProvider> RandomParticleMutProvider for T {
    unsafe fn get_particle_mut_unchecked(&mut self, i: usize) -> ParticleMut<'_> {
        let ind = self.get_index_unchecked(i);
        ParticleMut {
            id: ind,
            atom: &mut *(*self.sys_ptr_mut()).top.atoms.as_mut_ptr().add(ind),
            pos: &mut *(*self.sys_ptr_mut()).st.coords.as_mut_ptr().add(ind),
        }
    }
}

//=============================================================================
// Blanket measure/modify trait impls for SysProvider types

impl<T: SysProvider + IndexProvider> MeasurePos for T {}
impl<T: SysProvider + IndexProvider> MeasurePeriodic for T {}
impl<T: SysProvider + IndexProvider> MeasureMasses for T {}
impl<T: SysProvider + IndexProvider> MeasureRandomAccess for T {}
impl<T: SysProvider + IndexProvider> MeasureAtomPos for T {}

impl<T: SysMutProvider + IndexProvider> ModifyPos for T {}
impl<T: SysMutProvider + IndexProvider> ModifyPeriodic for T {}
impl<T: SysMutProvider + IndexProvider> ModifyRandomAccess for T {}
impl<T: SysMutProvider + IndexProvider> ModifyAtoms for T {}

//=============================================================================
/// Blanket split methods for types that have particle access.
/// Available for all SysProvider + IndexProvider types (they get RandomParticleProvider via blanket).
pub trait SplittableByParticle: RandomParticleProvider + IndexProvider + LenProvider + Sized {
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

    fn split_resindex(&self) -> impl Iterator<Item = Sel> {
        self.split(|p| Some(p.atom.resindex))
    }
}

// Blanket impl for all SysProvider + IndexProvider types
impl<T: SysProvider + IndexProvider> SplittableByParticle for T {}

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
        let all = (0..unsafe{&*self.sys_ptr()}.len()).into_iter().collect::<Vec<_>>();
        let index = unsafe {difference_sorted(&all, rhs.get_index_slice())};
        if index.is_empty() {
            return Err(SelectionError::EmptyDifference)
        }
        Ok(self.clone_with_index(index))
    }
}

//=============================================================================
// AtomPosAnalysis / NonAtomPosAnalysis - kept for backwards compat with molar_python

/// Umbrella trait for implementing read-only analysis traits involving only atoms and positions.
pub trait AtomPosAnalysis: LenProvider + IndexProvider + Sized {
    fn atoms_ptr(&self) -> *const Atom;
    fn coords_ptr(&self) -> *const Pos;

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

    fn split<RT, F>(&self, func: F) -> impl Iterator<Item = Sel>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(Particle) -> Option<RT>,
    {
        let mut cur_val = RT::default();
        let mut cur = 0usize;
        let ap = self.atoms_ptr() as usize;
        let cp = self.coords_ptr() as usize;

        let next_fn = move || {
            let mut index = Vec::<usize>::new();

            while cur < self.len() {
                let id = unsafe { self.get_index_unchecked(cur) };
                let p = unsafe {
                    Particle {
                        id,
                        atom: &*(ap as *const Atom).add(id),
                        pos: &*(cp as *const Pos).add(id),
                    }
                };
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

    fn whole_residues(&self) -> Sel {
        self.whole_attr(|at| &at.resindex)
    }

    fn whole_chains(&self) -> Sel {
        self.whole_attr(|at| &at.chain)
    }

    fn iter_atoms(&self) -> impl Iterator<Item = &Atom> {
        self.iter_index().map(|i| unsafe { &*self.atoms_ptr().add(i) })
    }

    fn iter_pos(&self) -> impl Iterator<Item = &Pos> {
        self.iter_index().map(|i| unsafe { &*self.coords_ptr().add(i) })
    }
}

/// Umbrella trait for implementing read-only analysis traits NOT involving atoms and positions
pub trait NonAtomPosAnalysis: LenProvider + IndexProvider + Sized {
    fn top_ptr(&self) -> *const Topology;
    fn st_ptr(&self) -> *const State;
}

/// Umbrella trait for implementing read-write analysis traits involving atoms and positions
pub trait AtomPosAnalysisMut: AtomPosAnalysis {
    fn atoms_ptr_mut(&mut self) -> *mut Atom {
        self.atoms_ptr() as *mut Atom
    }
    fn coords_ptr_mut(&mut self) -> *mut Pos {
        self.coords_ptr() as *mut Pos
    }
}

/// Umbrella trait for implementing read-write analysis traits NOT involving atoms and positions
pub trait NonAtomPosAnalysisMut: NonAtomPosAnalysis {
    fn top_ptr_mut(&mut self) -> *mut Topology {
        self.top_ptr() as *mut Topology
    }
    fn st_ptr_mut(&mut self) -> *mut State {
        self.st_ptr() as *mut State
    }
}

// PosIterProvider / AtomIterProvider / PosParIterProvider / AtomParIterProvider for SysProvider + IndexProvider types
impl<T: SysProvider + IndexProvider> PosIterProvider for T {
    fn iter_pos(&self) -> impl Iterator<Item = &Pos> {
        let len = self.len();
        (0..len).map(move |i| unsafe { self.pos_unchecked(self.get_index_unchecked(i)) })
    }
}

impl<T: SysProvider + IndexParProvider + Sync> PosParIterProvider for T {
    fn par_iter_pos(&self) -> impl IndexedParallelIterator<Item = &Pos> {
        let p = self as *const T as usize;
        self.par_iter_index().map(move |i| unsafe {
            let s = &*(p as *const T);
            s.pos_unchecked(i)
        })
    }
}

impl<T: SysProvider + IndexProvider> AtomIterProvider for T {
    fn iter_atoms(&self) -> impl Iterator<Item = &Atom> {
        let len = self.len();
        (0..len).map(move |i| unsafe { self.atom_unchecked(self.get_index_unchecked(i)) })
    }
}

impl<T: SysProvider + IndexParProvider + Sync> AtomParIterProvider for T {
    fn par_iter_atoms(&self) -> impl IndexedParallelIterator<Item = &Atom> {
        let p = self as *const T as usize;
        self.par_iter_index().map(move |i| unsafe {
            let s = &*(p as *const T);
            s.atom_unchecked(i)
        })
    }
}

// Note: AtomPosAnalysis types (e.g. molar_python SelPy, TmpSel) do NOT get PosIterProvider /
// AtomIterProvider via blanket (would conflict with SysProvider blanket above). Instead, they must
// implement PosIterProvider and AtomIterProvider explicitly in their crate.
// They have iter_pos / iter_atoms via AtomPosAnalysis directly, but for MeasurePos etc.
// they need to explicitly impl PosIterProvider using AtomPosAnalysis::iter_pos/iter_atoms.
//
// Note: MeasurePos / MeasureMasses / MeasurePeriodic / MeasureRandomAccess / MeasureAtomPos
// and ModifyPos / ModifyPeriodic / ModifyRandomAccess / ModifyAtoms for AtomPosAnalysis types
// must be implemented directly in each concrete type (e.g., SelPy, SystemPy in molar_python).
// Blanket impls for these marker traits are only provided for SysProvider+IndexProvider types
// (see the blanket impls block above).
