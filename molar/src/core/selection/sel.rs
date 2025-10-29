use std::{
    cell::{Ref, RefCell, RefMut},
    path::Path,
    rc::Rc,
};

use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator};

use crate::prelude::*;

//===========================================================================

pub trait Guarded {
    type Guard<'a>
    where
        Self: 'a;
    fn bind(&self) -> Result<Self::Guard<'_>, SelectionError>;
}

pub trait GuardedMut {
    type GuardMut<'a>
    where
        Self: 'a;
    fn bind_mut(&self) -> Result<Self::GuardMut<'_>, SelectionError>;
}

//===========================================================================
/// Selection is just an index detached from the Topology and State.
/// Has to be bound to System before doing something with the subset of atoms.
/// It is guaranteed to be non-empty so no run-time checks for this
/// are needed when binding to system
//===========================================================================
pub struct Sel {
    sys: Rc<RefCell<SystemStorage>>,
    index: SVec,
}

impl Sel {
    pub fn get_system(&self) -> SystemBound<'_> {
        SystemBound {
            guard: self.sys.borrow(),
            sys: &self.sys,
        }
    }

    pub fn get_system_mut(&mut self) -> SystemBoundMut<'_> {
        SystemBoundMut {
            guard: self.sys.borrow_mut(),
        }
    }
    
    pub fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        Ok(self.bind()?.select(def)?)
    }

    pub fn select_bound(&self, def: impl SelectionDef) -> Result<SelBoundOwn<'_>, SelectionError> {
        let g = self.sys.borrow();
        Ok(SelBoundOwn {
            guard: self.sys.borrow(),
            index: def.into_sel_index(&g.top, &g.st, Some(&self.index))?,
            sys: Rc::clone(&self.sys),
        })
    }

    /// Creates a string in Gromacs index format.
    pub fn as_gromacs_ndx_str(&self, name: impl AsRef<str>) -> String {
        use itertools::Itertools;
        let name = name.as_ref();
        let mut s = format!("[ {} ]\n", name);
        for chunk in &self.index.iter().chunks(15) {
            let line: String = chunk.map(|i| (i + 1).to_string()).join(" ");
            s.push_str(&line);
            s.push('\n');
        }
        s
    }

    pub fn save(&self, fname: &str) -> Result<(), FileIoError> {
        self.bind().unwrap().save(fname)
    }
}

impl Guarded for Sel {
    type Guard<'a> = SelBound<'a>;

    fn bind(&self) -> Result<Self::Guard<'_>, SelectionError> {
        let guard = self.sys.try_borrow()?;
        let last_ind = self.index.len() - 1;
        if unsafe { *self.index.get_unchecked(last_ind) } < guard.top.len() {
            Ok(SelBound {
                guard,
                index: &self.index,
                sys: &self.sys,
            })
        } else {
            Err(SelectionError::OutOfBounds(0, 0))
        }
    }
}

impl GuardedMut for Sel {
    type GuardMut<'a> = SelBoundMut<'a>;

    fn bind_mut(&self) -> Result<Self::GuardMut<'_>, SelectionError> {
        let guard = self.sys.try_borrow_mut()?;
        let last_ind = self.index.len() - 1;
        if unsafe { *self.index.get_unchecked(last_ind) } < guard.top.len() {
            Ok(SelBoundMut {
                guard,
                index: &self.index,
                _sys: &self.sys,
            })
        } else {
            Err(SelectionError::OutOfBounds(0, 0))
        }
    }
}

impl LenProvider for Sel {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for Sel {
    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }

    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }
}

macro_rules! impl_measure_and_modify_for_guarded {
    ( $t:ty ) => {
        impl MeasurePos for $t {
            fn min_max(&self) -> (Pos, Pos) {
                self.bind().unwrap().min_max()
            }

            fn center_of_geometry(&self) -> Pos {
                self.bind().unwrap().center_of_geometry()
            }

            fn rmsd(&self, other: &Self) -> Result<f32, MeasureError> {
                let g1 = self.bind().unwrap();
                let g2 = other.bind().unwrap();
                g1.rmsd(&g2)
            }
        }

        impl MeasureMasses for $t {
            fn center_of_mass(&self) -> Result<Pos, MeasureError> {
                self.bind().unwrap().center_of_mass()
            }

            fn gyration(&self) -> Result<f32, MeasureError> {
                self.bind().unwrap().gyration()
            }

            fn inertia(&self) -> Result<(Vector3f, Matrix3f), MeasureError> {
                self.bind().unwrap().inertia()
            }

            fn principal_transform(&self) -> Result<nalgebra::IsometryMatrix3<f32>, MeasureError> {
                self.bind().unwrap().principal_transform()
            }
        }

        impl MeasurePeriodic for $t {
            fn center_of_geometry_pbc(&self) -> Result<Pos, MeasureError> {
                self.bind().unwrap().center_of_geometry_pbc()
            }

            fn center_of_geometry_pbc_dims(&self, dims: PbcDims) -> Result<Pos, MeasureError> {
                self.bind().unwrap().center_of_geometry_pbc_dims(dims)
            }

            fn center_of_mass_pbc(&self) -> Result<Pos, MeasureError> {
                self.bind().unwrap().center_of_mass_pbc()
            }

            fn center_of_mass_pbc_dims(&self, dims: PbcDims) -> Result<Pos, MeasureError> {
                self.bind().unwrap().center_of_mass_pbc_dims(dims)
            }

            fn gyration_pbc(&self) -> Result<f32, MeasureError> {
                self.bind().unwrap().gyration_pbc()
            }

            fn inertia_pbc(&self) -> Result<(Vector3f, Matrix3f), MeasureError> {
                self.bind().unwrap().inertia_pbc()
            }

            fn principal_transform_pbc(
                &self,
            ) -> Result<nalgebra::IsometryMatrix3<f32>, MeasureError> {
                self.bind().unwrap().principal_transform_pbc()
            }
        }

        impl MeasureRandomAccess for $t {
            fn lipid_tail_order(
                &self,
                order_type: OrderType,
                normals: &Vec<Vector3f>,
                bond_orders: &Vec<u8>,
            ) -> Result<nalgebra::DVector<f32>, LipidOrderError> {
                self.bind()
                    .unwrap()
                    .lipid_tail_order(order_type, normals, bond_orders)
            }
        }

        impl ModifyPos for $t {
            fn apply_transform(&mut self, tr: &nalgebra::IsometryMatrix3<f32>) {
                self.bind_mut().unwrap().apply_transform(tr);
            }

            fn rotate(&mut self, ax: &nalgebra::Unit<Vector3f>, ang: f32) {
                self.bind_mut().unwrap().rotate(ax, ang);
            }

            fn translate<S>(
                &mut self,
                shift: &nalgebra::Matrix<f32, nalgebra::Const<3>, nalgebra::Const<1>, S>,
            ) where
                S: nalgebra::storage::Storage<f32, nalgebra::Const<3>, nalgebra::Const<1>>,
            {
                self.bind_mut().unwrap().translate(shift);
            }
        }

        impl ModifyAtoms for $t {
            fn assign_resindex(&mut self) {
                self.bind_mut().unwrap().assign_resindex();
            }
        }

        impl ModifyPeriodic for $t {
            fn unwrap_simple_dim(&mut self, dims: PbcDims) -> Result<(), MeasureError> {
                self.bind_mut().unwrap().unwrap_simple_dim(dims)
            }
        }

        impl ModifyRandomAccess for $t {
            fn unwrap_connectivity_dim(
                &mut self,
                cutoff: f32,
                dims: PbcDims,
            ) -> Result<(), MeasureError> {
                self.bind_mut()
                    .unwrap()
                    .unwrap_connectivity_dim(cutoff, dims)
            }
        }
    };
}

impl_measure_and_modify_for_guarded!(Sel);
impl_measure_and_modify_for_guarded!(System);
//=================================================

/// Internal System storage behind the Rc
#[derive(Debug, Default)]
pub struct SystemStorage {
    top: Topology,
    st: State,
}

//================================================
/// System that stores Topology and State
//================================================
#[derive(Debug, Default)]
pub struct System(Rc<RefCell<SystemStorage>>);

impl System {
    pub fn new(top: Topology, st: State) -> Result<Self, SelectionError> {
        check_topology_state_sizes(&top, &st)?;
        Ok(Self(Rc::new(RefCell::new(SystemStorage { top, st }))))
    }

    pub fn new_empty() -> Self {
        Self(Rc::new(RefCell::new(SystemStorage::default())))
    }

    pub fn from_file(fname: impl AsRef<Path>) -> Result<Self, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        Ok(Self::new(top, st)?)
    }

    /// Create new selection based on provided definition.
    pub fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        let g = self.0.try_borrow()?;
        Ok(Sel {
            index: def.into_sel_index(&g.top, &g.st, None)?,
            sys: Rc::clone(&self.0),
        })
    }

    pub fn select_bound(&self, def: impl SelectionDef) -> Result<SelBoundOwn<'_>, SelectionError> {
        let g = self.0.borrow();
        let index = def.into_sel_index(&g.top, &g.st, None)?;
        Ok(SelBoundOwn {
            guard: g,
            index,
            sys: Rc::clone(&self.0),
        })
    }

    pub fn select_all(&self) -> Result<Sel, SelectionError> {
        let v: Vec<usize> = (0..self.len() - 1).into_iter().collect();
        let index = unsafe { SVec::from_sorted(v) };
        Ok(Sel {
            index,
            sys: Rc::clone(&self.0),
        })
    }

    pub fn set_state(&mut self, st: State) -> Result<State, SelectionError> {
        let mut g = self.0.try_borrow_mut()?;
        if !g.st.interchangeable(&st) {
            return Err(SelectionError::IncompatibleState);
        }
        Ok(std::mem::replace(&mut g.st, st))
    }

    pub fn set_topology(&mut self, top: Topology) -> Result<Topology, SelectionError> {
        let mut g = self.0.try_borrow_mut()?;
        if !g.top.interchangeable(&top) {
            return Err(SelectionError::IncompatibleTopology);
        }
        Ok(std::mem::replace(&mut g.top, top))
    }

    pub fn release(self) -> Result<(Topology, State), SelectionError> {
        let c = Rc::try_unwrap(self.0).map_err(|_| SelectionError::Release)?;
        let s = c.into_inner();
        Ok((s.top, s.st))
    }

    pub fn save(&self, fname: &str) -> Result<(), FileIoError> {
        self.bind().unwrap().save(fname)
    }

    pub fn append(&mut self, sel: &Sel) {
        if Rc::ptr_eq(&self.0, &sel.sys) {
            // The same system.
            // We can't extend atoms and coords while reading from them at the same time
            // because the buffers could be reallocated in the process
            // so we collect selection atoms and coords first and then add them
            // iterators over pos and atoms manually (unsafely)
            let sel_g = sel.bind().unwrap();
            let sel_pos: Vec<_> = sel_g.iter_pos().collect();
            let sel_atoms: Vec<_> = sel_g.iter_atoms().collect();
            // sel_g is not used below and thus is dropped here
            // Now obtain an exclusive lock and do appending
            let mut sys_g = self.bind_mut().unwrap();
            sys_g.guard.st.add_coords(sel_pos.into_iter().cloned());
            sys_g.guard.top.add_atoms(sel_atoms.into_iter().cloned());
        } else {
            // Different systems. No problem to hold two locks
            let mut sys_g = self.bind_mut().unwrap();
            let sel_g = sel.bind().unwrap();
            sys_g.guard.st.add_coords(sel_g.iter_pos().cloned());
            sys_g.guard.top.add_atoms(sel_g.iter_atoms().cloned());
        }
    }

    pub fn has_box(&self) -> bool {
        self.bind().unwrap().get_box().is_some()
    }
}

impl LenProvider for System {
    fn len(&self) -> usize {
        self.0.borrow().st.len()
    }
}

impl Guarded for System {
    type Guard<'a> = SystemBound<'a>;
    fn bind(&self) -> Result<Self::Guard<'_>, SelectionError> {
        Ok(SystemBound {
            guard: self.0.try_borrow()?,
            sys: &self.0,
        })
    }
}

impl GuardedMut for System {
    type GuardMut<'a> = SystemBoundMut<'a>;
    fn bind_mut(&self) -> Result<Self::GuardMut<'_>, SelectionError> {
        Ok(SystemBoundMut {
            guard: self.0.try_borrow_mut()?,
        })
    }
}

//================================================
pub struct SystemBound<'a> {
    guard: Ref<'a, SystemStorage>,
    sys: &'a Rc<RefCell<SystemStorage>>,
}

impl Guarded for SystemBound<'_> {
    type Guard<'a> = &'a Self
    where
        Self: 'a;
    
    fn bind(&self) -> Result<Self::Guard<'_>, SelectionError> {
        Ok(self)
    }
}

impl LenProvider for SystemBound<'_> {
    fn len(&self) -> usize {
        self.guard.st.len()
    }
}

impl IndexProvider for SystemBound<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        i
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        (0..self.len()).into_iter()
    }
}

impl SelectableGuard for SystemBound<'_> {
    fn sys_ref(&self) -> &Rc<RefCell<SystemStorage>> {
        self.sys
    }

    fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        Ok(Sel {
            index: def.into_sel_index(&self.guard.top, &self.guard.st, None)?,
            sys: Rc::clone(&self.sys),
        })
    }

    fn select_bound(&self, def: impl SelectionDef) -> Result<SelBoundOwn<'_>, SelectionError> {
        Ok(SelBoundOwn {
            guard: Ref::clone(&self.guard),
            index: def.into_sel_index(&self.guard.top, &self.guard.st, None)?,
            sys: Rc::clone(&self.sys),
        })
    }
}

impl AtomPosAnalysis for SystemBound<'_> {
    fn atom_ptr(&self) -> *const Atom {
        self.guard.top.atoms.as_ptr()
    }

    fn pos_ptr(&self) -> *const Pos {
        self.guard.st.coords.as_ptr()
    }
}

impl NonAtomPosAnalysis for SystemBound<'_> {
    fn top_ref(&self) -> &Topology {
        &self.guard.top
    }

    fn st_ref(&self) -> &State {
        &self.guard.st
    }
}

impl TopologyWrite for SystemBound<'_> {}
impl StateWrite for SystemBound<'_> {}
impl TopologyStateWrite for SystemBound<'_> {}

//================================================
pub struct SystemBoundMut<'a> {
    guard: RefMut<'a, SystemStorage>,
}

impl SystemBoundMut<'_> {
    pub fn set_state(&mut self, st: State) -> Result<State, SelectionError> {
        if !self.guard.st.interchangeable(&st) {
            return Err(SelectionError::IncompatibleState);
        }
        Ok(std::mem::replace(&mut self.guard.st, st))
    }

    pub fn set_topology(&mut self, top: Topology) -> Result<Topology, SelectionError> {
        if !self.guard.top.interchangeable(&top) {
            return Err(SelectionError::IncompatibleTopology);
        }
        Ok(std::mem::replace(&mut self.guard.top, top))
    }
}

impl GuardedMut for SystemBoundMut<'_> {
    type GuardMut<'a> = &'a Self
    where
        Self: 'a;

    fn bind_mut(&self) -> Result<Self::GuardMut<'_>, SelectionError> {
        Ok(self)
    }
}

impl LenProvider for SystemBoundMut<'_> {
    fn len(&self) -> usize {
        self.guard.st.len()
    }
}

impl IndexProvider for SystemBoundMut<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        i
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        (0..self.len()).into_iter()
    }
}

impl AtomPosAnalysis for SystemBoundMut<'_> {
    fn atom_ptr(&self) -> *const Atom {
        self.guard.top.atoms.as_ptr()
    }

    fn pos_ptr(&self) -> *const Pos {
        self.guard.st.coords.as_ptr()
    }
}

impl AtomPosAnalysisMut for SystemBoundMut<'_> {
    fn atom_mut_ptr(&mut self) -> *mut Atom {
        self.guard.top.atoms.as_mut_ptr()
    }

    fn pos_mut_ptr(&mut self) -> *mut Pos {
        self.guard.st.coords.as_mut_ptr()
    }
}

impl NonAtomPosAnalysis for SystemBoundMut<'_> {
    fn st_ref(&self) -> &State {
        &self.guard.st
    }

    fn top_ref(&self) -> &Topology {
        &self.guard.top
    }
}

impl NonAtomPosAnalysisMut for SystemBoundMut<'_> {
    fn st_ref_mut(&mut self) -> &mut State {
        &mut self.guard.st
    }

    fn top_ref_mut(&mut self) -> &mut Topology {
        &mut self.guard.top
    }
}

//================================================
/// Read only subsystem
/// Implements only read-only analysis traits
//================================================
pub struct SelBound<'a> {
    guard: Ref<'a, SystemStorage>,
    index: &'a SVec,
    sys: &'a Rc<RefCell<SystemStorage>>,
}

impl Guarded for SelBound<'_> {
    type Guard<'a> = &'a Self
    where
        Self: 'a;

    fn bind(&self) -> Result<Self::Guard<'_>, SelectionError> {
        Ok(self)
    }
}

impl<'a> SelBound<'a> {
    pub fn into_mut(self) -> Result<SelBoundMut<'a>, SelectionError> {
        drop(self.guard);
        Ok(SelBoundMut {
            index: self.index,
            _sys: self.sys,
            guard: self.sys.try_borrow_mut()?,
        })
    }
}

impl LenProvider for SelBound<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SelBound<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl SelectableGuard for SelBound<'_> {
    fn sys_ref(&self) -> &Rc<RefCell<SystemStorage>> {
        self.sys
    }

    fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        Ok(Sel {
            index: def.into_sel_index(&self.guard.top, &self.guard.st, Some(self.index))?,
            sys: Rc::clone(&self.sys),
        })
    }

    fn select_bound(&self, def: impl SelectionDef) -> Result<SelBoundOwn<'_>, SelectionError> {
        Ok(SelBoundOwn {
            guard: Ref::clone(&self.guard),
            index: def.into_sel_index(&self.guard.top, &self.guard.st, Some(self.index))?,
            sys: Rc::clone(&self.sys),
        })
    }
}

impl AtomPosAnalysis for SelBound<'_> {
    fn atom_ptr(&self) -> *const Atom {
        self.guard.top.atoms.as_ptr()
    }

    fn pos_ptr(&self) -> *const Pos {
        self.guard.st.coords.as_ptr()
    }
}

impl NonAtomPosAnalysis for SelBound<'_> {
    fn top_ref(&self) -> &Topology {
        &self.guard.top
    }

    fn st_ref(&self) -> &State {
        &self.guard.st
    }
}

impl TopologyWrite for SelBound<'_> {}
impl StateWrite for SelBound<'_> {}
impl TopologyStateWrite for SelBound<'_> {}

//================================================
/// Bound selection which owns its index.
/// Created by `select_bound()` methods of existing bound selections.
pub struct SelBoundOwn<'a> {
    guard: Ref<'a, SystemStorage>,
    index: SVec,
    sys: Rc<RefCell<SystemStorage>>,
}

impl Guarded for SelBoundOwn<'_> {
    type Guard<'a> = &'a Self
    where
        Self: 'a;

    fn bind(&self) -> Result<Self::Guard<'_>, SelectionError> {
        Ok(self)
    }
}

impl SelBoundOwn<'_> {
    /// Convert to unbound [Sel]
    pub fn into_unbound(self) -> Sel {
        Sel {
            index: self.index,
            sys: self.sys,
        }
    }
}

impl LenProvider for SelBoundOwn<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SelBoundOwn<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl SelectableGuard for SelBoundOwn<'_> {
    fn sys_ref(&self) -> &Rc<RefCell<SystemStorage>> {
        &self.sys
    }

    fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        Ok(Sel {
            index: def.into_sel_index(&self.guard.top, &self.guard.st, Some(&self.index))?,
            sys: Rc::clone(&self.sys),
        })
    }

    fn select_bound(&self, def: impl SelectionDef) -> Result<SelBoundOwn<'_>, SelectionError> {
        Ok(SelBoundOwn {
            guard: Ref::clone(&self.guard),
            index: def.into_sel_index(&self.guard.top, &self.guard.st, Some(&self.index))?,
            sys: Rc::clone(&self.sys),
        })
    }
}

impl AtomPosAnalysis for SelBoundOwn<'_> {
    fn atom_ptr(&self) -> *const Atom {
        self.guard.top.atoms.as_ptr()
    }

    fn pos_ptr(&self) -> *const Pos {
        self.guard.st.coords.as_ptr()
    }
}

impl NonAtomPosAnalysis for SelBoundOwn<'_> {
    fn top_ref(&self) -> &Topology {
        &self.guard.top
    }

    fn st_ref(&self) -> &State {
        &self.guard.st
    }
}

impl TopologyWrite for SelBoundOwn<'_> {}
impl StateWrite for SelBoundOwn<'_> {}
impl TopologyStateWrite for SelBoundOwn<'_> {}

//================================================
/// Read-write subsystem having access to all fields of Topology and State
//================================================
pub struct SelBoundMut<'a> {
    guard: RefMut<'a, SystemStorage>,
    index: &'a SVec,
    _sys: &'a Rc<RefCell<SystemStorage>>,
}

impl GuardedMut for SelBoundMut<'_> {
    type GuardMut<'a>
        = &'a Self
    where
        Self: 'a;
    fn bind_mut(&self) -> Result<Self::GuardMut<'_>, SelectionError> {
        Ok(self)
    }
}

impl<'a> SelBoundMut<'a> {
    pub fn into_immut(self) -> Result<SelBound<'a>, SelectionError> {
        drop(self.guard);
        Ok(SelBound {
            index: self.index,
            sys: self._sys,
            guard: self._sys.try_borrow()?,
        })
    }
}

impl LenProvider for SelBoundMut<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SelBoundMut<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl AtomPosAnalysis for SelBoundMut<'_> {
    fn atom_ptr(&self) -> *const Atom {
        self.guard.top.atoms.as_ptr()
    }

    fn pos_ptr(&self) -> *const Pos {
        self.guard.st.coords.as_ptr()
    }
}

impl AtomPosAnalysisMut for SelBoundMut<'_> {
    fn atom_mut_ptr(&mut self) -> *mut Atom {
        self.guard.top.atoms.as_mut_ptr()
    }

    fn pos_mut_ptr(&mut self) -> *mut Pos {
        self.guard.st.coords.as_mut_ptr()
    }
}

impl NonAtomPosAnalysis for SelBoundMut<'_> {
    fn st_ref(&self) -> &State {
        &self.guard.st
    }

    fn top_ref(&self) -> &Topology {
        &self.guard.top
    }
}

impl NonAtomPosAnalysisMut for SelBoundMut<'_> {
    fn st_ref_mut(&mut self) -> &mut State {
        &mut self.guard.st
    }

    fn top_ref_mut(&mut self) -> &mut Topology {
        &mut self.guard.top
    }
}

// impl TopologyWrite for SubSystemMut<'_> {}
// impl StateWrite for SubSystemMut<'_> {}
// impl TopologyStateWrite for SubSystemMut<'_> {}

//================================================
/// Read-write subsystem for non-blocking parallel access to atoms and posisitons
/// Doesn't have access to shared fields such as box and bonds.
//================================================
pub struct SubSystemParMut {
    pos_ptr: *mut Pos,
    atom_ptr: *mut Atom,
    index: SVec,
}
unsafe impl Sync for SubSystemParMut {}
unsafe impl Send for SubSystemParMut {}

impl LenProvider for SubSystemParMut {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SubSystemParMut {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl AtomPosAnalysis for SubSystemParMut {
    fn atom_ptr(&self) -> *const Atom {
        self.atom_ptr
    }

    fn pos_ptr(&self) -> *const Pos {
        self.pos_ptr
    }
}

impl AtomPosAnalysisMut for SubSystemParMut {
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
    selections: Vec<SubSystemParMut>,
    max_index: usize,
    sys: Rc<RefCell<SystemStorage>>,
}

/// Bound parallel split
pub struct ParSplitBound<'a> {
    pub selections: &'a mut Vec<SubSystemParMut>,
    _guard: RefMut<'a, SystemStorage>,
}

impl ParSplit {
    /// Binds split to System for read-write access
    pub fn bind_mut<'a>(&'a mut self) -> Result<ParSplitBound<'a>, SelectionError> {
        let mut g = self.sys.try_borrow_mut()?;
        // The cost calling is just one comparison
        if self.max_index < g.top.len() {
            // Set correct pointers in each selection
            self.selections.iter_mut().for_each(|sel| {
                sel.pos_ptr = g.st.coords.as_mut_ptr();
                sel.atom_ptr = g.top.atoms.as_mut_ptr();
            });
            Ok(ParSplitBound {
                selections: &mut self.selections,
                _guard: g,
            })
        } else {
            Err(SelectionError::OutOfBounds(0, 0))
        }
    }
}

impl ParSplitBound<'_> {
    /// Returns parallel iterator over stored parallel selections.
    pub fn iter_par(&self) -> rayon::slice::Iter<'_, SubSystemParMut> {
        self.selections.par_iter()
    }

    /// Returns parallel mutable iterator over stored parallel selections.
    pub fn iter_par_mut(&mut self) -> rayon::slice::IterMut<'_, SubSystemParMut> {
        self.selections.par_iter_mut()
    }

    /// Returns serial iterator over stored parallel selections.
    pub fn iter(&self) -> impl Iterator<Item = &SubSystemParMut> {
        self.selections.iter()
    }

    /// Returns serial mutable iterator over stored parallel selections.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut SubSystemParMut> {
        self.selections.iter_mut()
    }
}

//================================================
/// Umbrella trait for implementing read-only analysis traits involving atoms and positions
//================================================
pub trait SelectableGuard: AtomPosAnalysis {
    // Reference to the system
    fn sys_ref(&self) -> &Rc<RefCell<SystemStorage>>;

    fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError>;

    fn select_bound(&self, def: impl SelectionDef) -> Result<SelBoundOwn<'_>, SelectionError>;

    /// Creates a parallel split based on provided closure.
    /// A closure takes a [Particle] and returns a distinct value for each piece.
    /// New selection is created whenever new return value differes from the previous one.
    fn split_par<F, R>(&self, func: F) -> Result<ParSplit, SelectionError>
    where
        F: Fn(Particle) -> Option<R>,
        R: Default + PartialOrd,
        Self: Sized,
    {
        let svecs: Vec<SVec> = self.split_iter(func).collect();

        if svecs.is_empty() {
            return Err(SelectionError::EmptySplit);
        }

        // Compute last index
        let max_index = *svecs.iter().map(|v| v.last().unwrap()).max().unwrap();

        Ok(ParSplit {
            selections: svecs
                .into_iter()
                .map(|v| SubSystemParMut {
                    index: v,
                    atom_ptr: std::ptr::null_mut(),
                    pos_ptr: std::ptr::null_mut(),
                })
                .collect(),
            sys: Rc::clone(self.sys_ref()),
            max_index,
        })
    }

    // Internal splitting function
    fn split_iter<RT, F>(&self, func: F) -> impl Iterator<Item = SVec>
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
                        return Some(unsafe { SVec::from_sorted(index) });
                    }
                }
                // Next particle
                cur += 1;
            }

            // Return any remaining index as last selection
            if !index.is_empty() {
                return Some(unsafe { SVec::from_sorted(index) });
            }

            // If we are here stop iterating
            None
        };

        std::iter::from_fn(next_fn)
    }

    fn split_resindex_iter(&self) -> impl Iterator<Item = SVec> {
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

        Sel {
            index: unsafe { SVec::from_sorted(ind) },
            sys: Rc::clone(self.sys_ref()),
        }
    }

    /// Selects whole residiues present in the current selection (in terms of resindex)
    fn whole_residues(&self) -> Sel {
        self.whole_attr(|at| &at.resindex)
    }

    /// Selects whole chains present in the current selection
    fn whole_chains(&self) -> Sel {
        self.whole_attr(|at| &at.chain)
    }
}

pub trait AtomPosAnalysis: LenProvider + IndexProvider + Sized {
    // Raw pointer to the beginning of array of atoms
    fn atom_ptr(&self) -> *const Atom;
    // Raw pointer to the beginning of array of coords
    fn pos_ptr(&self) -> *const Pos;

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
pub trait NonAtomPosAnalysis: LenProvider + IndexProvider + Sized {
    fn top_ref(&self) -> &Topology;
    fn st_ref(&self) -> &State;

    /// Splits by molecule and returns an iterator over them.
    /// If molecule is only partially contained in self then only this part is returned (molecules are clipped).
    /// If there are no molecules in [Topology] return an empty iterator.
    fn split_mol_iter(&self) -> impl Iterator<Item = Sel>
    where
        Self: Sized + SelectableGuard,
    {
        let n = self.len();
        // Iterate over molecules and find those inside selection
        let first = unsafe { self.get_index_unchecked(0) };
        let last = unsafe { self.get_index_unchecked(n - 1) };

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
            }
            .map(|r| Sel {
                index: unsafe { SVec::from_sorted(r.into_iter().collect()) },
                sys: Rc::clone(&self.sys_ref()),
            });

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
    fn atom_mut_ptr(&mut self) -> *mut Atom;
    // Raw pointer to the beginning of array of coords
    fn pos_mut_ptr(&mut self) -> *mut Pos;
}

//================================================
/// Umbrella trait for implementing read-write analysis traits NOT involving atoms and positions
//================================================
pub trait NonAtomPosAnalysisMut: Sized {
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
    unsafe fn get_particle_mut_unchecked(&mut self, i: usize) -> ParticleMut<'_> {
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

//---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    //use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator};

    #[test]
    fn test1() -> anyhow::Result<()> {
        let sys = System::from_file("tests/protein.pdb")?;
        let sel1 = sys.select([1usize, 2, 6, 7])?;

        for at in sel1.bind()?.iter_atoms() {
            println!("{} {}", at.name, at.resname);
        }

        let mut lock = sel1.bind_mut()?;
        for at in lock.iter_atoms_mut() {
            at.bfactor += 1.0;
        }

        let lock2 = sel1.bind()?;
        for at in lock2.iter_atoms() {
            println!("{} ", at.bfactor);
        }

        //drop(lock);

        Ok(())
    }

    #[test]
    fn test1_1() -> anyhow::Result<()> {
        let sys = System::from_file("tests/albumin.pdb")?;

        let mut par = sys.bind()?.split_par(|p| {
            if p.atom.resname != "SOL" {
                Some(p.atom.resindex)
            } else {
                None
            }
        })?;

        par.bind_mut()?.iter_par_mut().try_for_each(|sel| {
            println!("{} {}", sel.len(), sel.first_atom().resname);
            if sel.first_atom().resname == "ALA" {
                sel.first_pos_mut().coords.add_scalar_mut(1.0);
                println!("{}", sel.first_pos());
            }
            Ok::<_, SelectionError>(())
        })?;

        // Add serial selection
        let ca = sys.select("name CA")?;
        let ca_b = ca.bind()?;
        let cb = sys.select("name CB")?;
        let cb_b = cb.bind()?;
        println!("#ca: {} {}", ca_b.len(), ca_b.center_of_mass()?);

        for a in cb_b.iter_atoms() {
            println!("{}", a.name);
        }

        //Iter serial views
        let b = par.bind_mut()?;
        let serials: Vec<_> = b.iter().collect();
        println!("serial #5: {}", serials[5].first_atom().resname);

        Ok(())
    }
}
