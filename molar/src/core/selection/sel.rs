use std::path::Path;

use nalgebra::{Const, IsometryMatrix3, Unit};
use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator};
use thiserror::Error;

use crate::prelude::*;

#[derive(Error,Debug)]
pub enum BindError {
    #[error("can't bind selection immutably: last index {0} is out of bounds (0:{1})")]
    Immut(usize,usize),
    #[error("can't bind selection mutably: last index {0} is out of bounds (0:{1})")]
    Mut(usize,usize),
}
//===========================================================================
/// Selection is just an index detached from the Topology and State.
/// Has to be bound to System before doing something with the subset of atoms.
/// It is guaranteed to be non-empty so no run-time checks for this
/// are needed when binding to system
//===========================================================================
pub struct Sel(pub(crate) SVec);

impl Sel {
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

    fn from_iter(iter: impl Iterator<Item=usize>) -> Self {
        let v = iter.collect();
        Self(unsafe{SVec::from_sorted(v)})
    }

    /// Binds Sel to System for read-only access
    pub fn try_bind<'a>(&'a self, sys: &'a System) -> Result<SubSystem<'a>, BindError> {
        let ind = self.0.len() - 1;
        let last_ind = unsafe { *self.0.get_unchecked(ind) };
        if last_ind < sys.top.len() {
            Ok(SubSystem {
                sys: &sys,
                index: &self.0,
            })
        } else {
            Err(BindError::Immut(last_ind, sys.top.len()))
        }
    }

    pub fn bind<'a>(&'a self, sys: &'a System) -> SubSystem<'a> {
        self.try_bind(sys).expect("binding selection immutably should not fail")
    }

    /// Binds Sel to System for read-write access
    pub fn try_bind_mut<'a>(&'a self, sys: &'a mut System) -> Result<SubSystemMut<'a>, BindError> {
        // The cost calling is just one comparison
        let ind = self.0.len() - 1;
        let last_ind = unsafe { *self.0.get_unchecked(ind) };
        if last_ind < sys.top.len() {
            Ok(SubSystemMut {
                sys,
                index: &self.0,
            })
        } else {
            Err(BindError::Mut(last_ind, sys.top.len()))
        }
    }

    pub fn bind_mut<'a>(&'a self, sys: &'a mut System) -> SubSystemMut<'a> {
        self.try_bind_mut(sys).expect("binding selection mutably should not fail")
    }
}

impl LenProvider for Sel {
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl IndexProvider for Sel {
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

    /// Create new selection based on provided definition.
    pub fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        Ok(Sel(def.into_sel_index(&self.top, &self.st, None)?))
    }

    pub fn select_all(&self) -> Sel {
        Sel(unsafe {SVec::from_sorted( (0..self.len()).into_iter().collect() )})
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

    /// Binds Sel to System for read-only access
    pub fn with<'a,T>(&'a self, sel: &'a Sel, op: T) -> Result<(), SelectionError> 
    where 
        T: Fn(SubSystem<'a>),
    {
        op(sel.try_bind(self)?);
        Ok(())
    }

    /// Binds Sel to System for read-write access
    pub fn with_mut<'a,T>(&'a mut self, sel: &'a Sel, mut op: T) -> Result<(), SelectionError>
    where 
        T: FnMut(SubSystemMut<'a>),
    {
        op(sel.try_bind_mut(self)?);
        Ok(())
    }

    //===============
    // Measuring
    //===============

    pub fn min_max(&self, sel: &Sel) -> Result<(Pos, Pos), MeasureError> {
        Ok(sel.try_bind(self)?.min_max())
    }

    pub fn center_of_geometry(&self, sel: &Sel) -> Result<Pos,MeasureError> {
        Ok(sel.try_bind(self)?.center_of_geometry())
    }

    pub fn center_of_geometry_pbc(&self, sel: &Sel) -> Result<Pos,MeasureError> {
        Ok(sel.try_bind(self)?.center_of_geometry_pbc()?)
    }

    pub fn center_of_geometry_pbc_dims(&self, sel: &Sel, dims: PbcDims) -> Result<Pos,MeasureError> {
        Ok(sel.try_bind(self)?.center_of_geometry_pbc_dims(dims)?)
    }

    pub fn rmsd(&self, sel1: &Sel, sel2: &Sel) -> Result<f32, MeasureError> {
        let b1 = sel1.try_bind(self)?;
        let b2 = sel2.try_bind(self)?;
        Ok(MeasurePos::rmsd(&b1, &b2)?)
    }

    pub fn rmsd_mw(&self, sel1: &Sel, sel2: &Sel) -> Result<f32, MeasureError> {
        let b1 = sel1.try_bind(self)?;
        let b2 = sel2.try_bind(self)?;
        Ok(rmsd_mw(&b1, &b2)?)
    }

    pub fn center_of_mass(&self, sel: &Sel) -> Result<Pos,MeasureError> {
        Ok(sel.try_bind(self)?.center_of_mass()?)
    }

    pub fn center_of_mass_pbc(&self, sel: &Sel) -> Result<Pos,MeasureError> {
        Ok(sel.try_bind(self)?.center_of_mass_pbc()?)
    }

    pub fn center_of_mass_pbc_dims(&self, sel: &Sel, dims: PbcDims) -> Result<Pos,MeasureError> {
        Ok(sel.try_bind(self)?.center_of_mass_pbc_dims(dims)?)
    }

    pub fn gyration(&self, sel: &Sel) -> Result<f32,MeasureError> {
        Ok(sel.try_bind(self)?.gyration()?)
    }

    pub fn gyration_pbc(&self, sel: &Sel) -> Result<f32,MeasureError> {
        Ok(sel.try_bind(self)?.gyration_pbc()?)
    }

    pub fn inertia(&self, sel: &Sel) -> Result<(Vector3f, Matrix3f), MeasureError> {
        Ok(sel.try_bind(self)?.inertia()?)
    }

    pub fn inertia_pbc(&self, sel: &Sel) -> Result<(Vector3f, Matrix3f), MeasureError> {
        Ok(sel.try_bind(self)?.inertia_pbc()?)
    }

    pub fn principal_transform(&self, sel: &Sel) -> Result<IsometryMatrix3<f32>, MeasureError> {
        Ok(sel.try_bind(self)?.principal_transform()?)
    }

    pub fn principal_transform_pbc(&self, sel: &Sel) -> Result<IsometryMatrix3<f32>, MeasureError> {
        Ok(sel.try_bind(self)?.principal_transform_pbc()?)
    }

    pub fn fit_transform(&self, sel1: &Sel, sel2: &Sel) -> Result<IsometryMatrix3<f32>, MeasureError> {
        let b1 = sel1.try_bind(self)?;
        let b2 = sel2.try_bind(self)?;
        Ok(fit_transform(&b1, &b2)?)
    }

    pub fn fit_transform_at_origin(&self, sel1: &Sel, sel2: &Sel) -> Result<IsometryMatrix3<f32>, MeasureError> {
        let b1 = sel1.try_bind(self)?;
        let b2 = sel2.try_bind(self)?;
        Ok(fit_transform_at_origin(&b1, &b2)?)
    }

    pub fn lipid_tail_order(
        &self,
        sel: &Sel,
        order_type: OrderType,
        normals: &Vec<Vector3f>,
        bond_orders: &Vec<u8>,
    ) -> Result<nalgebra::DVector<f32>, MeasureError> {
        Ok(sel.try_bind(self)?.lipid_tail_order(order_type,normals,bond_orders)?)
    }

    //===============
    // Modifying
    //===============

    pub fn translate<S>(&mut self, sel: &Sel, shift: &nalgebra::Matrix<f32, Const<3>, Const<1>, S>) -> Result<(),MeasureError>
    where
        S: nalgebra::storage::Storage<f32, Const<3>, Const<1>>,
    {
        Ok(sel.try_bind_mut(self)?.translate(shift))
    }

    pub fn rotate(&mut self, sel: &Sel, ax: &Unit<Vector3f>, ang: f32) -> Result<(),MeasureError> {
        Ok(sel.try_bind_mut(self)?.rotate(ax,ang))
    }

    pub fn apply_transform(&mut self, sel: &Sel, tr: &nalgebra::IsometryMatrix3<f32>) -> Result<(),MeasureError> {
        Ok(sel.try_bind_mut(self)?.apply_transform(tr))
    }

    pub fn unwrap_simple_dim(&mut self, sel: &Sel, dims: PbcDims) -> Result<(), MeasureError> {
        Ok(sel.try_bind_mut(self)?.unwrap_simple_dim(dims)?)
    }

    pub fn unwrap_simple(&mut self, sel: &Sel) -> Result<(), MeasureError> {
        Ok(sel.try_bind_mut(self)?.unwrap_simple()?)
    }

    pub fn unwrap_connectivity(&mut self, sel: &Sel, cutoff: f32) -> Result<(), MeasureError> {
        Ok(sel.try_bind_mut(self)?.unwrap_connectivity(cutoff)?)
    }

    pub fn unwrap_connectivity_dim(&mut self, sel: &Sel, cutoff: f32, dims: PbcDims) -> Result<(), MeasureError> {
        Ok(sel.try_bind_mut(self)?.unwrap_connectivity_dim(cutoff,dims)?)
    }

    // Saving
    pub fn save_sel(&self, sel: &Sel, fname: impl AsRef<Path>) -> Result<(),FileIoError> {
        Ok(sel.try_bind(self).map_err(|e| FileIoError(fname.as_ref().to_path_buf(),FileFormatError::Bind(e)))?
        .save(fname.as_ref().to_str().unwrap())?)
    }

    /// Computes the Solvet Accessible Surface Area (SASA).
    pub fn sasa(&self, sel: &Sel) -> Result<SasaResults, MeasureError> {
        let b = sel.try_bind(self)?;
        Ok(molar_powersasa::compute_sasa(
            b.len(),
            0.14,
            |i| unsafe {
                let ind = b.get_index_unchecked(i);
                b.coords_ptr().add(ind) as *mut f32
            },
            |i: usize| b.get_particle(i).unwrap().atom.vdw(),
        ))
    }

    /// Append selection derived from self
    pub fn append_self_sel(&mut self, sel: &Sel) -> Result<Sel,SelectionError> {
        let old_last = self.len()-1;
        let ind = sel.0.len() - 1;
        let last_ind = unsafe { *sel.0.get_unchecked(ind) };
        if last_ind >= self.top.len() {
            return Err(BindError::Mut(last_ind, self.top.len()))?;
        }
        let pos: Vec<_> = sel.0.iter().map(|i| &self.st.coords[*i]).cloned().collect();
        let atoms: Vec<_> = sel.0.iter().map(|i| &self.top.atoms[*i]).cloned().collect();
        self.st.add_coords(pos.into_iter());
        self.top.add_atoms(atoms.into_iter());
        Ok(Sel::from_iter(old_last+1..self.len()))
    }

    pub fn append_atoms_pos<'a>(&mut self, atoms: impl AtomIterator<'a>, coords: impl PosIterator<'a>) -> Result<Sel,SelectionError> {
        let old_last = self.len()-1;
        self.st.add_coords(coords.cloned());
        self.top.add_atoms(atoms.cloned());
        check_topology_state_sizes(&self.top, &self.st)?;
        Ok(Sel::from_iter(old_last+1..self.len()))
    }

    pub fn append(&mut self, data: &(impl AtomIterProvider + PosIterProvider)) -> Result<Sel,SelectionError> {
        let old_last = self.len()-1;
        self.st.add_coords(data.iter_pos().cloned());
        self.top.add_atoms(data.iter_atoms().cloned());
        check_topology_state_sizes(&self.top, &self.st)?;
        Ok(Sel::from_iter(old_last+1..self.len()))
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
                    let added = self.append_self_sel(&all)?;
                    let shift =
                        m.column(0) * x as f32 + m.column(1) * y as f32 + m.column(2) * z as f32;
                    added.bind_mut(self).translate(&shift);
                }
            }
        }
        // Scale the box
        self.get_box_mut().unwrap().scale_vectors([nbox[0] as f32, nbox[1] as f32, nbox[2] as f32])?;

        // Re-assign resindex
        self.top.assign_resindex();
        Ok(())
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
/// Read only subsystem
/// Implements only read-only analysis traits
//================================================
pub struct SubSystem<'a> {
    sys: &'a System,
    index: &'a [usize],
}

impl SubSystem<'_> {
    /// Create new sub-selection based on provided definition.
    pub fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        Ok(Sel(def.into_sel_index(
            &self.sys.top,
            &self.sys.st,
            Some(self.index),
        )?))
    }
}

impl LenProvider for SubSystem<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SubSystem<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl AtomPosAnalysis for SubSystem<'_> {
    fn atoms_ptr(&self) -> *const Atom {
        self.sys.top.atoms.as_ptr()
    }

    fn coords_ptr(&self) -> *const Pos {
        self.sys.st.coords.as_ptr()
    }
}

impl NonAtomPosAnalysis for SubSystem<'_> {
    fn top_ptr(&self) -> *const Topology {
        &self.sys.top
    }

    fn st_ptr(&self) -> *const State {
        &self.sys.st
    }
}

impl TopologyWrite for SubSystem<'_> {}
impl StateWrite for SubSystem<'_> {}
impl TopologyStateWrite for SubSystem<'_> {}
//================================================
/// Read-write subsystem having access to all fields of Topology and State
//================================================
pub struct SubSystemMut<'a> {
    sys: &'a mut System,
    index: &'a [usize],
}

impl SubSystemMut<'_> {
    /// Create new sub-selection based on provided definition.
    pub fn select(&self, def: impl SelectionDef) -> Result<Sel, SelectionError> {
        Ok(Sel(def.into_sel_index(
            &self.sys.top,
            &self.sys.st,
            Some(self.index),
        )?))
    }
}

impl LenProvider for SubSystemMut<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SubSystemMut<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl AtomPosAnalysis for SubSystemMut<'_> {
    fn atoms_ptr(&self) -> *const Atom {
        self.sys.top.atoms.as_ptr()
    }

    fn coords_ptr(&self) -> *const Pos {
        self.sys.st.coords.as_ptr()
    }
}

impl AtomPosAnalysisMut for SubSystemMut<'_> {
    fn atoms_ptr_mut(&mut self) -> *mut Atom {
        self.sys.top.atoms.as_mut_ptr()
    }

    fn coords_ptr_mut(&mut self) -> *mut Pos {
        self.sys.st.coords.as_mut_ptr()
    }
}

impl NonAtomPosAnalysis for SubSystemMut<'_> {
    fn st_ptr(&self) -> *const State {
        &self.sys.st
    }

    fn top_ptr(&self) -> *const Topology {
        &self.sys.top
    }
}


impl NonAtomPosAnalysisMut for SubSystemMut<'_> {
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
pub struct SubSystemParMut<'a> {
    pos_ptr: *mut Pos,
    atom_ptr: *mut Atom,
    index: &'a [usize],
}
unsafe impl Sync for SubSystemParMut<'_> {}
unsafe impl Send for SubSystemParMut<'_> {}

impl LenProvider for SubSystemParMut<'_> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

impl IndexProvider for SubSystemParMut<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        *self.index.get_unchecked(i)
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        self.index.iter().cloned()
    }
}

impl AtomPosAnalysis for SubSystemParMut<'_> {
    fn atoms_ptr(&self) -> *const Atom {
        self.atom_ptr
    }

    fn coords_ptr(&self) -> *const Pos {
        self.pos_ptr
    }
}

impl AtomPosAnalysisMut for SubSystemParMut<'_> {
    fn atoms_ptr_mut(&mut self) -> *mut Atom {
        self.atom_ptr
    }

    fn coords_ptr_mut(&mut self) -> *mut Pos {
        self.pos_ptr
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

/// Bound parallel split
pub struct ParSplitBound<'a>(Vec<SubSystemParMut<'a>>);

impl ParSplit {
    /// Binds split to System for read-write access
    pub fn bind_mut<'a>(
        &'a self,
        sys: &'a mut System,
    ) -> Result<ParSplitBound<'a>, SelectionError> {
        // The cost calling is just one comparison
        if self.max_index < sys.top.len() {
            // Bind each selection in mutable parallel way
            let bound_sels: Vec<_> = self
                .selections
                .iter()
                .map(|sel| SubSystemParMut {
                    index: &sel.0,
                    pos_ptr: sys.st.coords.as_mut_ptr(),
                    atom_ptr: sys.top.atoms.as_mut_ptr(),
                })
                .collect();
            Ok(ParSplitBound(bound_sels))
        } else {
            Err(SelectionError::OutOfBounds(0, 0))
        }
    }
}

impl<'a> ParSplitBound<'a> {
    /// Returns parallel iterator over stored parallel selections.
    pub fn iter_par(&'a self) -> rayon::slice::Iter<'a, SubSystemParMut<'a>> {
        self.0.par_iter()
    }

    /// Returns parallel mutable iterator over stored parallel selections.
    pub fn iter_par_mut(&'a mut self) -> rayon::slice::IterMut<'a, SubSystemParMut<'a>> {
        self.0.par_iter_mut()
    }

    /// Returns serial iterator over stored parallel selections.
    pub fn iter(&'a self) -> impl Iterator<Item = &'a SubSystemParMut<'a>> {
        self.0.iter()
    }

    /// Returns serial mutable iterator over stored parallel selections.
    pub fn iter_mut(&'a mut self) -> impl Iterator<Item = &'a mut SubSystemParMut<'a>> {
        self.0.iter_mut()
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

        for at in sel1.try_bind(&sys)?.iter_atoms() {
            println!("{} {}", at.name, at.resname);
        }

        let mut lock = sel1.try_bind_mut(&mut sys)?;
        for at in lock.iter_atoms_mut() {
            at.bfactor += 1.0;
        }

        let lock2 = sel1.try_bind(&sys)?;
        for at in lock2.iter_atoms() {
            println!("{} ", at.bfactor);
        }

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

        par.bind_mut(&mut sys)?.iter_par_mut().try_for_each(|sel| {
            println!("{} {}", sel.len(), sel.first_atom().resname);
            if sel.first_atom().resname == "ALA" {
                sel.first_pos_mut().coords.add_scalar_mut(1.0);
                println!("{}", sel.first_pos());
            }
            Ok::<_, SelectionError>(())
        })?;

        // Add serial selection
        let ca = sys.select("name CA")?;
        let ca_b = ca.try_bind(&sys)?;
        let cb = sys.select("name CB")?;
        let cb_b = cb.try_bind(&sys)?;
        println!("#ca: {} {}", ca_b.len(), ca_b.center_of_mass()?);

        for a in cb_b.iter_atoms() {
            println!("{}", a.name);
        }

        //Iter serial views
        let b = par.bind_mut(&mut sys)?;
        let serials: Vec<_> = b.iter().collect();
        println!("serial #5: {}", serials[5].first_atom().resname);

        Ok(())
    }
}
