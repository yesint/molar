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
pub struct Sel(SVec);

impl Sel {
    pub fn from_vec(index: Vec<usize>) -> Result<Self, SelectionError> {
        if index.is_empty() {
            Err(SelectionError::EmptySlice)
        } else {
            Ok(Self(SVec::from_unsorted(index)))
        }
    }

    /// Binds Sel to System for read-only access
    pub fn bind<'a>(&'a self, sys: &'a System) -> Result<SubSystem<'a>, BindError> {
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

    /// Binds Sel to System for read-write access
    pub fn bind_mut<'a>(&'a self, sys: &'a mut System) -> Result<SubSystemMut<'a>, BindError> {
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

    /// Creates a string in Gromacs index format representing self.
    pub fn as_gromacs_ndx_str(&self, name: impl AsRef<str>) -> String {
        use itertools::Itertools;
        let name = name.as_ref();
        let mut s = format!("[ {} ]\n", name);
        for chunk in &self.0.iter().chunks(15) {
            let line: String = chunk.map(|i| (i + 1).to_string()).join(" ");
            s.push_str(&line);
            s.push('\n');
        }
        s
    }
}

impl LenProvider for Sel {
    fn len(&self) -> usize {
        self.0.len()
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

    // /// Binds Sel to System for read-only access
    // pub fn with<'a>(&'a self, sel: &'a Sel) -> Result<SubSystem<'a>, SelectionError> {
    //     let last_ind = sel.0.len() - 1;
    //     if unsafe { *sel.0.get_unchecked(last_ind) } < self.top.len() {
    //         Ok(SubSystem {
    //             sys: self,
    //             index: &sel.0,
    //         })
    //     } else {
    //         Err(SelectionError::OutOfBounds(0, 0))
    //     }
    // }

    // /// Binds Sel to System for read-write access
    // pub fn with_mut<'a>(&'a mut self, sel: &'a Sel) -> Result<SubSystemMut<'a>, SelectionError> {
    //     // The cost calling is just one comparison
    //     let last_ind = sel.0.len() - 1;
    //     if unsafe { *sel.0.get_unchecked(last_ind) } < self.top.len() {
    //         Ok(SubSystemMut {
    //             sys: self,
    //             index: &sel.0,
    //         })
    //     } else {
    //         Err(SelectionError::OutOfBounds(0, 0))
    //     }
    // }

    // Measuring

    pub fn min_max(&self, sel: &Sel) -> Result<(Pos, Pos), MeasureError> {
        Ok(sel.bind(self)?.min_max())
    }

    pub fn center_of_geometry(&self, sel: &Sel) -> Result<Pos,MeasureError> {
        Ok(sel.bind(self)?.center_of_geometry())
    }

    pub fn center_of_geometry_pbc(&self, sel: &Sel) -> Result<Pos,MeasureError> {
        Ok(sel.bind(self)?.center_of_geometry_pbc()?)
    }

    pub fn center_of_geometry_pbc_dims(&self, sel: &Sel, dims: PbcDims) -> Result<Pos,MeasureError> {
        Ok(sel.bind(self)?.center_of_geometry_pbc_dims(dims)?)
    }

    pub fn rmsd(&self, sel1: &Sel, sel2: &Sel) -> Result<f32, MeasureError> {
        let b1 = sel1.bind(self)?;
        let b2 = sel2.bind(self)?;
        Ok(MeasurePos::rmsd(&b1, &b2)?)
    }

    pub fn rmsd_mw(&self, sel1: &Sel, sel2: &Sel) -> Result<f32, MeasureError> {
        let b1 = sel1.bind(self)?;
        let b2 = sel2.bind(self)?;
        Ok(rmsd_mw(&b1, &b2)?)
    }

    pub fn center_of_mass(&self, sel: &Sel) -> Result<Pos,MeasureError> {
        Ok(sel.bind(self)?.center_of_mass()?)
    }

    pub fn center_of_mass_pbc(&self, sel: &Sel) -> Result<Pos,MeasureError> {
        Ok(sel.bind(self)?.center_of_mass_pbc()?)
    }

    pub fn center_of_mass_pbc_dims(&self, sel: &Sel, dims: PbcDims) -> Result<Pos,MeasureError> {
        Ok(sel.bind(self)?.center_of_mass_pbc_dims(dims)?)
    }

    pub fn gyration(&self, sel: &Sel) -> Result<f32,MeasureError> {
        Ok(sel.bind(self)?.gyration()?)
    }

    pub fn gyration_pbc(&self, sel: &Sel) -> Result<f32,MeasureError> {
        Ok(sel.bind(self)?.gyration_pbc()?)
    }

    pub fn inertia(&self, sel: &Sel) -> Result<(Vector3f, Matrix3f), MeasureError> {
        Ok(sel.bind(self)?.inertia()?)
    }

    pub fn inertia_pbc(&self, sel: &Sel) -> Result<(Vector3f, Matrix3f), MeasureError> {
        Ok(sel.bind(self)?.inertia_pbc()?)
    }

    pub fn principal_transform(&self, sel: &Sel) -> Result<IsometryMatrix3<f32>, MeasureError> {
        Ok(sel.bind(self)?.principal_transform()?)
    }

    pub fn principal_transform_pbc(&self, sel: &Sel) -> Result<IsometryMatrix3<f32>, MeasureError> {
        Ok(sel.bind(self)?.principal_transform_pbc()?)
    }

    pub fn fit_transform(&self, sel1: &Sel, sel2: &Sel) -> Result<IsometryMatrix3<f32>, MeasureError> {
        let b1 = sel1.bind(self)?;
        let b2 = sel2.bind(self)?;
        Ok(fit_transform(&b1, &b2)?)
    }

    pub fn fit_transform_at_origin(&self, sel1: &Sel, sel2: &Sel) -> Result<IsometryMatrix3<f32>, MeasureError> {
        let b1 = sel1.bind(self)?;
        let b2 = sel2.bind(self)?;
        Ok(fit_transform_at_origin(&b1, &b2)?)
    }

    pub fn lipid_tail_order(
        &self,
        sel: &Sel,
        order_type: OrderType,
        normals: &Vec<Vector3f>,
        bond_orders: &Vec<u8>,
    ) -> Result<nalgebra::DVector<f32>, MeasureError> {
        Ok(sel.bind(self)?.lipid_tail_order(order_type,normals,bond_orders)?)
    }

    // Modifying
    pub fn translate<S>(&mut self, sel: &Sel, shift: &nalgebra::Matrix<f32, Const<3>, Const<1>, S>) -> Result<(),MeasureError>
    where
        S: nalgebra::storage::Storage<f32, Const<3>, Const<1>>,
    {
        Ok(sel.bind_mut(self)?.translate(shift))
    }

    pub fn rotate(&mut self, sel: &Sel, ax: &Unit<Vector3f>, ang: f32) -> Result<(),MeasureError> {
        Ok(sel.bind_mut(self)?.rotate(ax,ang))
    }

    pub fn apply_transform(&mut self, sel: &Sel, tr: &nalgebra::IsometryMatrix3<f32>) -> Result<(),MeasureError> {
        Ok(sel.bind_mut(self)?.apply_transform(tr))
    }

    pub fn unwrap_simple_dim(&mut self, sel: &Sel, dims: PbcDims) -> Result<(), MeasureError> {
        Ok(sel.bind_mut(self)?.unwrap_simple_dim(dims)?)
    }

    pub fn unwrap_simple(&mut self, sel: &Sel) -> Result<(), MeasureError> {
        Ok(sel.bind_mut(self)?.unwrap_simple()?)
    }

    pub fn unwrap_connectivity(&mut self, sel: &Sel, cutoff: f32) -> Result<(), MeasureError> {
        Ok(sel.bind_mut(self)?.unwrap_connectivity(cutoff)?)
    }

    pub fn unwrap_connectivity_dim(&mut self, sel: &Sel, cutoff: f32, dims: PbcDims) -> Result<(), MeasureError> {
        Ok(sel.bind_mut(self)?.unwrap_connectivity_dim(cutoff,dims)?)
    }

    // Saving
    pub fn save_sel(&self, sel: &Sel, fname: impl AsRef<Path>) -> Result<(),FileIoError> {
        Ok(sel.bind(self).map_err(|e| FileIoError(fname.as_ref().to_path_buf(),FileFormatError::Bind(e)))?
        .save(fname.as_ref().to_str().unwrap())?)
    }

    /// Computes the Solvet Accessible Surface Area (SASA).
    pub fn sasa(&self, sel: &Sel) -> Result<SasaResults, MeasureError> {
        let b = sel.bind(self)?;
        Ok(molar_powersasa::compute_sasa(
            b.len(),
            0.14,
            |i| unsafe {
                let ind = b.get_index_unchecked(i);
                b.pos_ptr().add(ind) as *mut f32
            },
            |i: usize| b.get_particle(i).unwrap().atom.vdw(),
        ))
    }

    /// Append selection derived from self
    pub fn append_sel(&mut self, sel: &Sel) -> Result<(),SelectionError> {
        let ind = sel.0.len() - 1;
        let last_ind = unsafe { *sel.0.get_unchecked(ind) };
        if last_ind >= self.top.len() {
            return Err(BindError::Mut(last_ind, self.top.len()))?;
        }
        let pos: Vec<_> = sel.0.iter().map(|i| &self.st.coords[*i]).cloned().collect();
        let atoms: Vec<_> = sel.0.iter().map(|i| &self.top.atoms[*i]).cloned().collect();
        self.st.add_coords(pos.into_iter());
        self.top.add_atoms(atoms.into_iter());
        Ok(())
    }

    pub fn append<'a>(&mut self, coords: impl PosIterator<'a>, atoms: impl AtomIterator<'a>) -> Result<(),SelectionError> {
        self.st.add_coords(coords.cloned());
        self.top.add_atoms(atoms.cloned());
        check_topology_state_sizes(&self.top, &self.st)?;
        Ok(())
    }

    pub fn remove_coords(
        &mut self,
        removed: impl Iterator<Item = usize> + Clone,
    ) -> Result<(), BuilderError> {
        self.st.remove_coords(removed.clone())?;
        self.top.remove_atoms(removed)?;
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
    fn atom_ptr(&self) -> *const Atom {
        self.top.atoms.as_ptr()
    }

    fn pos_ptr(&self) -> *const Pos {
        self.st.coords.as_ptr()
    }
}

impl AtomPosAnalysisMut for System {
    fn atom_mut_ptr(&mut self) -> *mut Atom {
        self.top.atoms.as_mut_ptr()
    }

    fn pos_mut_ptr(&mut self) -> *mut Pos {
        self.st.coords.as_mut_ptr()
    }
}

impl NonAtomPosAnalysis for System {
    fn top_ref(&self) -> &Topology {
        &self.top
    }

    fn st_ref(&self) -> &State {
        &self.st
    }
}

impl NonAtomPosAnalysisMut for System {
    fn top_ref_mut(&mut self) -> &mut Topology {
        &mut self.top
    }

    fn st_ref_mut(&mut self) -> &mut State {
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
    fn atom_ptr(&self) -> *const Atom {
        self.sys.top.atoms.as_ptr()
    }

    fn pos_ptr(&self) -> *const Pos {
        self.sys.st.coords.as_ptr()
    }
}

impl NonAtomPosAnalysis for SubSystem<'_> {
    fn top_ref(&self) -> &Topology {
        &self.sys.top
    }

    fn st_ref(&self) -> &State {
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
    fn atom_ptr(&self) -> *const Atom {
        self.sys.top.atoms.as_ptr()
    }

    fn pos_ptr(&self) -> *const Pos {
        self.sys.st.coords.as_ptr()
    }
}

impl AtomPosAnalysisMut for SubSystemMut<'_> {
    fn atom_mut_ptr(&mut self) -> *mut Atom {
        self.sys.top.atoms.as_mut_ptr()
    }

    fn pos_mut_ptr(&mut self) -> *mut Pos {
        self.sys.st.coords.as_mut_ptr()
    }
}

impl NonAtomPosAnalysis for SubSystemMut<'_> {
    fn st_ref(&self) -> &State {
        &self.sys.st
    }

    fn top_ref(&self) -> &Topology {
        &self.sys.top
    }
}


impl NonAtomPosAnalysisMut for SubSystemMut<'_> {
    fn st_ref_mut(&mut self) -> &mut State {
        &mut self.sys.st
    }

    fn top_ref_mut(&mut self) -> &mut Topology {
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
    fn atom_ptr(&self) -> *const Atom {
        self.atom_ptr
    }

    fn pos_ptr(&self) -> *const Pos {
        self.pos_ptr
    }
}

impl AtomPosAnalysisMut for SubSystemParMut<'_> {
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
    selections: Vec<Sel>,
    max_index: usize,
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

//================================================
/// Umbrella trait for implementing read-only analysis traits involving atoms and positions
//================================================
pub(crate) trait AtomPosAnalysis: LenProvider + IndexProvider + Sized {
    // Raw pointer to the beginning of array of atoms
    fn atom_ptr(&self) -> *const Atom;
    // Raw pointer to the beginning of array of coords
    fn pos_ptr(&self) -> *const Pos;

    /// Creates a parallel split based on provided closure.
    /// A closure takes a [Particle] and returns a distinct value for each piece.
    /// New selection is created whenever new return value differes from the previous one.
    fn split_par<F, R>(&self, func: F) -> Result<ParSplit, SelectionError>
    where
        F: Fn(Particle) -> Option<R>,
        R: Default + PartialOrd,
        Self: Sized,
    {
        let selections: Vec<Sel> = self.split_iter(func).collect();

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
    fn split_iter<RT, F>(&self, func: F) -> impl Iterator<Item = Sel>
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

    fn split_resindex_iter(&self) -> impl Iterator<Item = Sel> {
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
                self.pos_ptr().add(ind) as *mut f32
            },
            |i: usize| self.get_particle(i).unwrap().atom.vdw(),
        )
    }
}

//================================================
/// Umbrella trait for implementing read-only analysis traits NOT involving atoms and positions
//================================================
trait NonAtomPosAnalysis: LenProvider + IndexProvider + Sized {
    fn top_ref(&self) -> &Topology;
    fn st_ref(&self) -> &State;

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
trait AtomPosAnalysisMut: LenProvider + IndexProvider + Sized {
    // Raw pointer to the beginning of array of atoms
    fn atom_mut_ptr(&mut self) -> *mut Atom;
    // Raw pointer to the beginning of array of coords
    fn pos_mut_ptr(&mut self) -> *mut Pos;
}

//================================================
/// Umbrella trait for implementing read-write analysis traits NOT involving atoms and positions
//================================================
trait NonAtomPosAnalysisMut: Sized {
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

//██████  Measure traits

impl<T: AtomPosAnalysis> MeasurePos for T {}
impl<T: AtomPosAnalysis + NonAtomPosAnalysis> MeasurePeriodic for T {}
impl<T: AtomPosAnalysis> MeasureMasses for T {}
impl<T: AtomPosAnalysis> MeasureRandomAccess for T {}

//██████  Modify traits

impl<T: AtomPosAnalysisMut> ModifyPos for T {}
impl<T: AtomPosAnalysisMut + NonAtomPosAnalysisMut + NonAtomPosAnalysis> ModifyPeriodic for T {}
impl<T: AtomPosAnalysisMut + AtomPosAnalysis + NonAtomPosAnalysisMut + NonAtomPosAnalysis> ModifyRandomAccess for T {}

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
        let sel1 = Sel::from_vec(vec![1, 2, 6, 7])?;

        for at in sel1.bind(&sys)?.iter_atoms() {
            println!("{} {}", at.name, at.resname);
        }

        let mut lock = sel1.bind_mut(&mut sys)?;
        for at in lock.iter_atoms_mut() {
            at.bfactor += 1.0;
        }

        let lock2 = sel1.bind(&sys)?;
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
        let ca_b = ca.bind(&sys)?;
        let cb = sys.select("name CB")?;
        let cb_b = cb.bind(&sys)?;
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
