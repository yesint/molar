use crate::{
    core::selection::sel::private::{AllowsSelecting, SelectionPrivate},
    prelude::*,
};
use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator};
use std::{marker::PhantomData, path::Path};
use triomphe::Arc;

// Private module containing internal methods implemented by all selection types
mod private {
    use crate::core::{SVec, State, Topology};
    use triomphe::Arc;

    // Trait for things that allow creating selections from them
    //(System, Sel, SelPar)
    pub trait AllowsSelecting {
        fn get_topology_arc(&self) -> &Arc<Topology>;
        fn set_topology_arc(&mut self, top: Arc<Topology>);
        fn get_state_arc(&self) -> &Arc<State>;
        fn set_state_arc(&mut self, st: Arc<State>);
        fn index_slice(&self) -> Option<&[usize]>;
    }

    // Internal trait for selections
    pub trait SelectionPrivate: AllowsSelecting {
        fn index_arc(&self) -> &Arc<SVec>;

        // Calling new is forbidden for users, thus it is in the private trait
        fn new_sel(topology: Arc<Topology>, state: Arc<State>, index: Arc<SVec>) -> Self
        where
            Self: Sized;
    }
}

//----------------------------------------------------------
/// Trait for things that can provide refs to Topology and State (Systems and Selections)
/// Used for implementing analysis traits.
pub trait HasTopState: LenProvider + IndexProvider {
    fn get_topology(&self) -> &Topology;
    fn get_state(&self) -> &State;
}

/// Trait for things that supports selecting atoms (Systems and Selections)
pub trait Selectable: private::AllowsSelecting + LenProvider + IndexProvider {
    type DerivedSel: private::SelectionPrivate + Selectable;

    /// Sets new [Topology] and returns a shared pointer to the old one
    /// This is "shallow" operation that only affects current object and doesn't influence
    /// other objects that refer to the same [Topology].
    fn set_topology(
        &mut self,
        topology: impl Into<Arc<Topology>>,
    ) -> Result<Arc<Topology>, SelectionError> {
        let topology: Arc<Topology> = topology.into();
        if !self.get_topology_arc().interchangeable(&topology) {
            return Err(SelectionError::IncompatibleState);
        }

        let ret = Arc::clone(self.get_topology_arc());
        self.set_topology_arc(topology);
        Ok(ret)
    }

    /// Sets new [State] and returns a shared pointer to the old one
    /// This is "shallow" operation that only affects current object and doesn't influence
    /// other objects that refer to the same [State].
    fn set_state(&mut self, state: impl Into<Arc<State>>) -> Result<Arc<State>, SelectionError> {
        let state: Arc<State> = state.into();
        if !self.get_state_arc().interchangeable(&state) {
            return Err(SelectionError::IncompatibleState);
        }
        
        let ret = Arc::clone(self.get_state_arc());
        self.set_state_arc(state);
        Ok(ret)
    }

    /// Sets new [State] grabbed from other [Selectable] object
    /// This is "shallow" operation that only affects current object and doesn't influence
    /// other objects that refer to the same [State].
    fn set_state_from(&mut self, other: &impl Selectable) -> Result<Arc<State>, SelectionError> {
        let cloned_arc = Arc::clone(other.get_state_arc());
        Ok(self.set_state(cloned_arc)?)
    }

    /// Create new selection based on provided definition.
    /// If self is already a selection then sub-selection is performed.
    fn select(&self, def: impl SelectionDef) -> Result<Self::DerivedSel, SelectionError>
    where
        Self: Sized,
    {
        use private::SelectionPrivate;

        let ind = def.into_sel_index(
            self.get_topology_arc(),
            self.get_state_arc(),
            self.index_slice(),
        )?;
        Ok(Self::DerivedSel::new_sel(
            Arc::clone(&self.get_topology_arc()),
            Arc::clone(&self.get_state_arc()),
            Arc::new(ind),
        ))
    }

    /// Creates a parallel split based on provided closure.
    /// A closure takes a [Particle] and returns a distinct value for each piece.
    /// New selection is created whenever new return value differes from the previous one.
    fn split_par<F, R>(&self, func: F) -> Result<ParSplit, SelectionError>
    where
        F: Fn(Particle) -> Option<R>,
        R: Default + PartialOrd,
        Self: Sized,
    {
        let selections: Vec<SelPar> = split_iter_as(self, func).collect();

        if selections.is_empty() {
            return Err(SelectionError::EmptySplit);
        }

        Ok(ParSplit {
            selections,
            _phantom: Default::default(),
        })
    }

    /// Split based on provided closure and return iterator over pieces.
    /// A closure takes a [Particle] and returns a distinct value for each piece.
    /// New selection is created whenever new return value differes from the previous one.
    fn split_iter<RT, F>(&self, func: F) -> impl Iterator<Item = Self::DerivedSel>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(Particle) -> Option<RT>,
        Self: Sized,
        Self::DerivedSel: Selectable,
    {
        split_iter_as(self, func)
    }

    fn split_resindex_iter(&self) -> impl Iterator<Item = Self::DerivedSel>
    where
        Self: Sized,
    {
        split_iter_as(self, |p| Some(p.atom.resindex))
    }

    /// Splits by molecule and returns an iterator over them. 
    /// If molecule is only partially contained in self then only this part is returned (molecules are clipped).
    /// If there are no molecules in [Topology] return an empty iterator.
    fn split_mol_iter(&self) -> impl Iterator<Item = Self::DerivedSel>
    where
        Self: Sized,
    {
        // Iterate over molecules and find those inside selection
        let (first, last) = match self.index_slice() {
            Some(ind) => (*ind.first().unwrap(), *ind.last().unwrap()),
            None => (0, self.len()),
        };

        let mut molid = 0;

        let next_fn = move || {
            if self.get_topology_arc().num_molecules() == 0 {
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
            .map(|r| self.select(r).unwrap());

            molid += 1;
            res
        };

        std::iter::from_fn(next_fn)
    }

    /// Creates a string in Gromacs index format representing self.
    fn as_gromacs_ndx_str(&self, name: impl AsRef<str>) -> String {
        use itertools::Itertools;
        let name = name.as_ref();
        let mut s = format!("[ {} ]\n", name);
        for chunk in &self.iter_index().chunks(15) {
            let line: String = chunk.map(|i| (i + 1).to_string()).join(" ");
            s.push_str(&line);
            s.push('\n');
        }
        s
    }
}

/// Trait for selections
pub trait Selection: private::SelectionPrivate + Selectable {
    /// Get selection index as a slice
    fn get_index(&self) -> &[usize] {
        self.index_slice().unwrap()
    }

    /// Returns first selection index
    fn get_first_index(&self) -> usize {
        *self.index_slice().unwrap().first().unwrap()
    }

    /// Returns last selection index
    fn get_last_index(&self) -> usize {
        *self.index_slice().unwrap().last().unwrap()
    }

    /// Creates new view of the current selection. 
    /// All views share the same index and thus are very cheap to create (no allocations needed).
    fn new_view(&self) -> Self::DerivedSel {
        Self::DerivedSel::new_sel(
            Arc::clone(self.get_topology_arc()),
            Arc::clone(self.get_state_arc()),
            Arc::clone(self.index_arc()),
        )
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
    fn whole_attr<T>(&self, attr_fn: fn(&Atom) -> &T) -> Self::DerivedSel
    where
        T: Eq + std::hash::Hash + Copy,
        Self: Sized,
    {
        // Collect all properties from the inner
        let mut properties = std::collections::HashSet::<T>::new();
        for at in self.iter_atoms() {
            properties.insert(*attr_fn(at));
        }

        let mut ind = vec![];
        // Loop over all atoms with the same property
        for (i, at) in self.get_topology().iter_atoms().enumerate() {
            let cur_prop = attr_fn(at);
            if properties.contains(cur_prop) {
                ind.push(i);
            }
        }

        Self::DerivedSel::new_sel(
            Arc::clone(&self.get_topology_arc()),
            Arc::clone(&self.get_state_arc()),
            Arc::new(unsafe { SVec::from_sorted(ind) }),
        )
    }

    /// Selects whole residiues present in the current selection (in terms of resindex)
    fn whole_residues(&self) -> Self::DerivedSel
    where
        Self: Sized,
    {
        self.whole_attr(|at| &at.resindex)
    }

    /// Selects whole chains present in the current selection
    fn whole_chains(&self) -> Self::DerivedSel
    where
        Self: Sized,
    {
        self.whole_attr(|at| &at.chain)
    }

    /// Computes the Solvet Accessible Surface Area (SASA).
    fn sasa(&self) -> SasaResults
    where
        Self: Sized,
    {
        molar_powersasa::compute_sasa(
            self.len(),
            0.14,
            |i| unsafe {
                let ind = *self.index_arc().get_unchecked(i);
                self.get_state_arc()
                    .get_pos_mut_unchecked(ind)
                    .coords
                    .as_mut_ptr()
            },
            |i: usize| self.get_particle(i).unwrap().atom.vdw(),
        )
    }
}

/// Trait for mutable selections
pub trait MutableSelectable: HasTopState {
    /// Sets same name to all selected atoms
    fn set_same_name(&self, val: &str)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.name = val.into();
        }
    }

    /// Sets same resname to all selected atoms
    fn set_same_resname(&self, val: &str)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.resname = val.into();
        }
    }

    /// Sets same resid to all selected atoms
    fn set_same_resid(&self, val: i32)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.resid = val;
        }
    }

    /// Sets same chain to all selected atoms
    fn set_same_chain(&self, val: char)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.chain = val;
        }
    }

    /// Sets same mass to all selected atoms
    fn set_same_mass(&self, val: f32)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.mass = val;
        }
    }

    /// Sets same B-factor to all selected atoms
    fn set_same_bfactor(&self, val: f32)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.bfactor = val;
        }
    }
}

impl<T: Selectable> HasTopState for T {
    fn get_topology(&self) -> &Topology {
        self.get_topology_arc()
    }

    fn get_state(&self) -> &State {
        self.get_state_arc()
    }
}

// Internal splitting function
fn split_iter_as<RT, F, S, SO>(sel: &S, func: F) -> impl Iterator<Item = SO> + use<'_, RT, F, S, SO>
where
    RT: Default + std::cmp::PartialEq,
    F: Fn(Particle) -> Option<RT>,
    S: Selectable,
    SO: private::SelectionPrivate,
{
    let mut cur_val = RT::default();
    let mut cur = 0usize;

    let next_fn = move || {
        let mut index = Vec::<usize>::new();

        while cur < sel.len() {
            let p = unsafe { sel.get_particle_unchecked(cur) };
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
                    return Some(SO::new_sel(
                        Arc::clone(sel.get_topology_arc()),
                        Arc::clone(sel.get_state_arc()),
                        Arc::new(unsafe { SVec::from_sorted(index) }),
                    ));
                }
            }
            // Next particle
            cur += 1;
        }

        // Return any remaining index as last selection
        if !index.is_empty() {
            return Some(SO::new_sel(
                Arc::clone(sel.get_topology_arc()),
                Arc::clone(sel.get_state_arc()),
                Arc::new(unsafe { SVec::from_sorted(index) }),
            ));
        }

        // If we are here stop iterating
        None
    };

    std::iter::from_fn(next_fn)
}

//-----------------------------------------

macro_rules! impl_selection {
    ($t:ty, $d:ident) => {
        impl Selectable for $t {
            type DerivedSel = $d;
        }

        impl private::AllowsSelecting for $t {
            fn index_slice(&self) -> Option<&[usize]> {
                Some(&self.index_storage)
            }

            fn get_state_arc(&self) -> &Arc<State> {
                &self.state
            }

            fn set_state_arc(&mut self, st: Arc<State>) {
                self.state = st;
            }

            fn get_topology_arc(&self) -> &Arc<Topology> {
                &self.topology
            }

            fn set_topology_arc(&mut self, top: Arc<Topology>) {
                self.topology = top;
            }
        }

        impl LenProvider for $t {
            fn len(&self) -> usize {
                self.index_storage.len()
            }
        }

        impl IndexProvider for $t {
            fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
                self.index_storage.iter().cloned()
            }

            unsafe fn get_index_unchecked(&self, i: usize) -> usize {
                *self.index_storage.get_unchecked(i)
            }
        }
    };
}

//-----------------------------------------
/// Serial selection. 
/// This a primary selection type that should be used by default. 
/// Most of [Sel] functionality is provided by implemented traits.
pub struct Sel {
    topology: Arc<Topology>,
    state: Arc<State>,
    index_storage: Arc<SVec>,
    _phantom: PhantomData<*const ()>,
}

impl MutableSelectable for Sel {}

impl private::SelectionPrivate for Sel {
    fn index_arc(&self) -> &Arc<SVec> {
        &self.index_storage
    }

    fn new_sel(topology: Arc<Topology>, state: Arc<State>, index: Arc<SVec>) -> Self
    where
        Self: Sized,
    {
        Self {
            topology,
            state,
            index_storage: index,
            _phantom: Default::default(),
        }
    }
}

impl Selection for Sel {}

impl_selection!(Sel, Sel);

impl Sel {
    pub fn to_par_immut(self) -> SelParImmut {
        SelParImmut {
            topology: self.topology,
            state: self.state,
            index_storage: self.index_storage,
        }
    }
}

//-------------------------------------------------------
/// Immutable parallel selection
pub struct SelParImmut {
    topology: Arc<Topology>,
    state: Arc<State>,
    index_storage: Arc<SVec>,
}

impl private::SelectionPrivate for SelParImmut {
    fn index_arc(&self) -> &Arc<SVec> {
        &self.index_storage
    }

    fn new_sel(topology: Arc<Topology>, state: Arc<State>, index: Arc<SVec>) -> Self
    where
        Self: Sized,
    {
        Self {
            topology,
            state,
            index_storage: index,
        }
    }
}

impl Selection for SelParImmut {}

impl_selection!(SelParImmut, SelParImmut);

impl SelParImmut {
    pub unsafe fn as_mut(&self) -> &Sel {
        std::mem::transmute(self)
    }
}

//-------------------------------------------------------

/// System is a container for matching [Topology] and [State].
/// Most of [System] functionality is provided by implemented traits.
pub struct System {
    topology: Arc<Topology>,
    state: Arc<State>,
    _phantom: PhantomData<*const ()>,
}

impl private::AllowsSelecting for System {
    fn index_slice(&self) -> Option<&[usize]> {
        None
    }

    fn get_state_arc(&self) -> &Arc<State> {
        &self.state
    }

    fn set_state_arc(&mut self, st: Arc<State>) {
        self.state = st;
    }

    fn get_topology_arc(&self) -> &Arc<Topology> {
        &self.topology
    }

    fn set_topology_arc(&mut self, top: Arc<Topology>) {
        self.topology = top;
    }
}

impl Selectable for System {
    type DerivedSel = Sel;
}

impl MutableSelectable for System {}

impl LenProvider for System {
    fn len(&self) -> usize {
        self.state.len()
    }
}

impl IndexProvider for System {
    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        0..self.len()
    }

    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        i
    }
}

impl System {
    pub fn new(
        top: impl Into<Arc<Topology>>,
        st: impl Into<Arc<State>>,
    ) -> Result<Self, SelectionError> {
        let top: Arc<Topology> = top.into();
        let st: Arc<State> = st.into();
        check_topology_state_sizes(&top, &st)?;
        Ok(Self {
            topology: top,
            state: st,
            _phantom: Default::default(),
        })
    }

    pub fn from_file(fname: impl AsRef<Path>) -> Result<Self, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        Ok(Self::new(top, st)?)
    }

    pub fn new_empty() -> Self {
        Self {
            topology: Arc::new(Topology::default()),
            state: Arc::new(State::default()),
            _phantom: Default::default(),
        }
    }

    /// Release and return [Topology] and [State].
    /// Fails if any selection is still pointing to them.
    pub fn release(self) -> Result<(Topology, State), SelectionError> {
        if self.topology.is_unique() && self.state.is_unique() {
            Ok((
                Arc::unwrap_or_clone(self.topology),
                Arc::unwrap_or_clone(self.state),
            ))
        } else {
            Err(SelectionError::Release)
        }
    }

    pub fn select_all(&self) -> Result<Sel, SelectionError> {
        self.select(0..self.len())
    }

    pub fn append(&self, data: &(impl PosIterProvider + AtomIterProvider)) -> Sel {
        let first_added_index = self.len();
        self.topology
            .get_storage_mut()
            .add_atoms(data.iter_atoms().cloned());
        self.state
            .get_storage_mut()
            .add_coords(data.iter_pos().cloned());
        let last_added_index = self.len();
        self.select(first_added_index..last_added_index).unwrap()
    }

    pub fn append_atoms_pos(
        &self,
        atoms: impl Iterator<Item = Atom>,
        coords: impl Iterator<Item = Pos>,
    ) -> Sel {
        let first_added_index = self.len();
        self.topology.get_storage_mut().add_atoms(atoms);
        self.state.get_storage_mut().add_coords(coords);
        let last_added_index = self.len();
        self.select(first_added_index..last_added_index).unwrap()
    }

    // This method only works if this system has no selections
    pub fn remove(&self, to_remove: &impl IndexProvider) -> Result<(), SelectionError> {
        if self.topology.is_unique() && self.state.is_unique() {
            self.topology
                .get_storage_mut()
                .remove_atoms(to_remove.iter_index())?;
            self.state
                .get_storage_mut()
                .remove_coords(to_remove.iter_index())?;
            Ok(())
        } else {
            Err(SelectionError::Release)
        }
    }

    /// Sets periodic box replacing current one
    pub fn set_box(&self, new_box: Option<PeriodicBox>) {
        self.state.get_storage_mut().pbox = new_box;
    }

    /// Replace periodic box with the box from other object
    pub fn set_box_from(&self, box_provider: &impl BoxProvider) {
        self.state.get_storage_mut().pbox = box_provider.get_box().cloned();
    }

    pub fn set_time(&self, t: f32) {
        self.state.set_time(t);
    }

    pub fn multiply_periodically(&self, nbox: [usize; 3]) -> Result<(), SelectionError> {
        if self.get_box().is_none() {
            return Err(PeriodicBoxError::NoPbc)?;
        }
        let b = self.get_box_mut().unwrap();
        let m = b.get_matrix();
        let all = self.select_all()?;
        for x in 0..=nbox[0] {
            for y in 0..=nbox[1] {
                for z in 0..=nbox[2] {
                    if x == 0 && y == 0 && z == 0 {
                        continue;
                    }
                    let added = self.append(&all);
                    let shift =
                        m.column(0) * x as f32 + m.column(1) * y as f32 + m.column(2) * z as f32;
                    added.translate(&shift);
                }
            }
        }
        // Scale the box
        b.scale_vectors([nbox[0] as f32, nbox[1] as f32, nbox[2] as f32])?;
        // Re-assign resindex
        self.topology.assign_resindex();
        Ok(())
    }
}

//---------------------------------------------------

pub struct SelPar {
    topology: Arc<Topology>,
    state: Arc<State>,
    index_storage: Arc<SVec>,
}

impl private::SelectionPrivate for SelPar {
    fn index_arc(&self) -> &Arc<SVec> {
        &self.index_storage
    }

    fn new_sel(topology: Arc<Topology>, state: Arc<State>, index: Arc<SVec>) -> Self
    where
        Self: Sized,
    {
        Self {
            topology,
            state,
            index_storage: index,
        }
    }
}

impl Selection for SelPar {}

impl MutableSelectable for SelPar {}

impl_selection!(SelPar, Sel);

//--------------------------------------------------------------------

macro_rules! impl_logical_ops {
    ($t:ident) => {
        impl<S: Selection> std::ops::BitOr<&S> for &$t {
            type Output = <$t as Selectable>::DerivedSel;
            
            /// Creates new selection which is a logical OR (union) of operands
            fn bitor(self, rhs: &S) -> Self::Output {
                let ind = union_sorted(&self.index_storage, rhs.index_arc());
                Self::Output::new_sel(
                    Arc::clone(self.get_topology_arc()),
                    Arc::clone(self.get_state_arc()),
                    Arc::new(ind),
                )
            }
        }

        impl<S: Selection> std::ops::BitAnd<&S> for &$t {
            type Output = <$t as Selectable>::DerivedSel;
            
            /// Creates new selection which is a logical AND (intersection) of operands
            fn bitand(self, rhs: &S) -> Self::Output {
                let ind = intersection_sorted(&self.index_storage, rhs.index_arc());
                if ind.is_empty() {
                    panic!("'and' of selections is empty");
                }
                Self::Output::new_sel(
                    Arc::clone(self.get_topology_arc()),
                    Arc::clone(self.get_state_arc()),
                    Arc::new(ind),
                )
            }
        }

        impl<S: Selection> std::ops::Sub<&S> for &$t {
            type Output = <$t as Selectable>::DerivedSel;

            /// Creates new selection with atoms from rhs removed from self (logical difference).
            fn sub(self, rhs: &S) -> Self::Output {
                let ind = difference_sorted(&self.index_storage, rhs.index_arc());
                Self::Output::new_sel(
                    Arc::clone(self.get_topology_arc()),
                    Arc::clone(self.get_state_arc()),
                    Arc::new(ind),
                )
            }
        }

        impl<S: Selection> std::ops::Add<&S> for &$t {
            type Output = <$t as Selectable>::DerivedSel;

            /// Creates new selection which is a logical OR (union) of operands. The same as `self | rhs`.
            fn add(self, rhs: &S) -> Self::Output {
                let ind = union_sorted(&self.index_storage, rhs.index_arc());
                Self::Output::new_sel(
                    Arc::clone(self.get_topology_arc()),
                    Arc::clone(self.get_state_arc()),
                    Arc::new(ind),
                )
            }
        }

        impl std::ops::Not for &$t {
            type Output = <$t as Selectable>::DerivedSel;

            /// Inverts selection i.e. selects all atoms that were not selected.
            /// ## Panics
            /// Panics if inverted selection is empty (happens if you invert "all" selection).
            fn not(self) -> Self::Output {
                let ind = difference_sorted(
                    unsafe { &SVec::from_sorted((0..self.topology.len()).collect()) },
                    self.index_arc(),
                );
                if ind.is_empty() {
                    panic!("negated selection is empty");
                }
                Self::Output::new_sel(
                    Arc::clone(self.get_topology_arc()),
                    Arc::clone(self.get_state_arc()),
                    Arc::new(ind),
                )
            }
        }
    };
}

impl_logical_ops!(Sel);
impl_logical_ops!(SelPar);
impl_logical_ops!(SelParImmut);

//--------------------------------------------------------------------
pub struct ParSplit {
    selections: Vec<SelPar>,
    _phantom: PhantomData<*const ()>,
}

impl ParSplit {
    /// Creates a split from iterator over serial selections
    pub fn from_serial_selections<'a>(
        sels: impl IntoIterator<Item = &'a Sel>,
    ) -> Result<Self, SelectionError> {
        // For oerlap check
        let mut used = vec![];
        let mut selections = vec![];
        let mut top_arc: Option<&Arc<Topology>> = None;
        let mut st_arc: Option<&Arc<State>> = None;

        for sel in sels {
            // On first selection resize used
            if used.is_empty() {
                used.resize(sel.get_topology().len(), false);
                top_arc = Some(sel.get_topology_arc());
                st_arc = Some(sel.get_state_arc());
            }

            let top_arc = top_arc.unwrap();
            let st_arc = st_arc.unwrap();

            if !(Arc::ptr_eq(top_arc, sel.get_topology_arc())
                && Arc::ptr_eq(st_arc, sel.get_state_arc()))
            {
                return Err(SelectionError::ParSplitDifferentSystems);
            }

            for i in sel.iter_index() {
                if used[i] {
                    return Err(SelectionError::ParSplitOverlap);
                } else {
                    used[i] = true;
                }
            }

            selections.push(SelPar::new_sel(
                Arc::clone(top_arc),
                Arc::clone(st_arc),
                Arc::clone(sel.index_arc()),
            ));
        }
        Ok(Self {
            selections,
            _phantom: Default::default(),
        })
    }

    /// Returns parallel iterator over stored parallel selections.
    pub fn par_iter(&mut self) -> rayon::slice::Iter<'_, SelPar> {
        self.selections.par_iter()
    }

    /// Returns parallel mutable iterator over stored parallel selections.
    pub fn par_iter_mut(&mut self) -> rayon::slice::IterMut<'_, SelPar> {
        self.selections.par_iter_mut()
    }

    pub fn iter(&self) -> impl Iterator<Item = Sel> + '_ {
        self.selections.iter().map(|sel| sel.new_view())
    }
}

//═══════════════════════════════════════════════════════════
//  Blanket trait implementations for Selections and System
//═══════════════════════════════════════════════════════════

//██████  IO traits

impl<T: HasTopState> WritableToFile for T {}

impl<T: HasTopState> TopologyIoProvider for T {}
impl<T: HasTopState> StateIoProvider for T {}

impl<T: HasTopState> TimeProvider for T {
    fn get_time(&self) -> f32 {
        self.get_state().get_time()
    }
}

//██████  Immutable analysis traits

impl<T: HasTopState> BoxProvider for T {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.get_state().get_box()
    }
}

impl<T: HasTopState> PosIterProvider for T {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        unsafe {
            self.iter_index()
                .map(|i| self.get_state().get_pos_unchecked(i))
        }
    }
}

impl<T: HasTopState> AtomIterProvider for T {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        unsafe {
            self.iter_index()
                .map(|i| self.get_topology().get_atom_unchecked(i))
        }
    }
}

impl<T: HasTopState + MutableSelectable> AtomIterMutProvider for T {
    fn iter_atoms_mut(&self) -> impl AtomMutIterator<'_> {
        unsafe {
            self.iter_index()
                .map(|i| self.get_topology().get_atom_mut_unchecked(i))
        }
    }
}

impl<T: HasTopState> RandomPosProvider for T {
    unsafe fn get_pos_unchecked(&self, i: usize) -> &Pos {
        let ind = self.get_index_unchecked(i);
        self.get_state().get_pos_unchecked(ind)
    }
}

impl<T: HasTopState> RandomAtomProvider for T {
    unsafe fn get_atom_unchecked(&self, i: usize) -> &Atom {
        let ind = self.get_index_unchecked(i);
        self.get_topology().get_atom_unchecked(ind)
    }
}

impl<T: HasTopState + MutableSelectable> RandomAtomMutProvider for T {
    unsafe fn get_atom_mut_unchecked(&self, i: usize) -> &mut Atom {
        let ind = self.get_index_unchecked(i);
        self.get_topology().get_atom_mut_unchecked(ind)
    }
}

impl<T: HasTopState> MoleculesProvider for T {
    fn num_molecules(&self) -> usize {
        self.get_topology().num_molecules()
    }

    fn iter_molecules(&self) -> impl Iterator<Item = &[usize; 2]> {
        self.get_topology().iter_molecules()
    }

    unsafe fn get_molecule_unchecked(&self, i: usize) -> &[usize; 2] {
        self.get_topology().get_molecule_unchecked(i)
    }
}

impl<T: HasTopState> BondsProvider for T {
    fn num_bonds(&self) -> usize {
        self.get_topology().num_bonds()
    }

    fn iter_bonds(&self) -> impl Iterator<Item = &[usize; 2]> {
        self.get_topology().iter_bonds()
    }

    unsafe fn get_bond_unchecked(&self, i: usize) -> &[usize; 2] {
        self.get_topology().get_bond_unchecked(i)
    }
}

impl<T: HasTopState> RandomParticleProvider for T {
    unsafe fn get_particle_unchecked(&self, i: usize) -> Particle<'_> {
        let ind = self.get_index_unchecked(i);
        Particle {
            id: ind,
            atom: self.get_topology().get_atom_unchecked(ind),
            pos: self.get_state().get_pos_unchecked(ind),
        }
    }
}

//██████  Measure traits

impl<T: HasTopState> MeasurePos for T {}
impl<T: HasTopState> MeasurePeriodic for T {}
impl<T: HasTopState> MeasureMasses for T {}
impl<T: HasTopState> MeasureRandomAccess for T {}

//██████  Mutable analysis traits

impl BoxMutProvider for System {
    fn get_box_mut(&self) -> Option<&mut PeriodicBox> {
        self.get_state().get_box_mut()
    }
}

impl<T: HasTopState + MutableSelectable> PosIterMutProvider for T {
    fn iter_pos_mut(&self) -> impl PosMutIterator<'_> {
        unsafe {
            self.iter_index()
                .map(|i| self.get_state().get_pos_mut_unchecked(i))
        }
    }
}

impl<T: HasTopState + MutableSelectable> RandomPosMutProvider for T {
    unsafe fn get_pos_mut_unchecked(&self, i: usize) -> &mut Pos {
        let ind = self.get_index_unchecked(i);
        self.get_state().get_pos_mut_unchecked(ind)
    }
}

impl<T: HasTopState + MutableSelectable> RandomParticleMutProvider for T {
    unsafe fn get_particle_mut_unchecked(&self, i: usize) -> ParticleMut {
        let ind = self.get_index_unchecked(i);
        ParticleMut {
            id: ind,
            atom: self.get_topology().get_atom_mut_unchecked(ind),
            pos: self.get_state().get_pos_mut_unchecked(ind),
        }
    }
}

//██████  Modify traits

impl<T: HasTopState + MutableSelectable> ModifyPos for T {}
impl<T: HasTopState + MutableSelectable> ModifyPeriodic for T {}
impl<T: HasTopState + MutableSelectable> ModifyRandomAccess for T {}

//========================================================

#[allow(unused_imports)]
mod tests {
    use super::*;
    use crate::prelude::*;

    #[test]
    fn test1() -> anyhow::Result<()> {
        let (top, st) = FileHandler::open("tests/albumin.pdb")?.read()?;
        let sys = System::new(top, st)?;

        let mut par = sys.split_par(|p| {
            if p.atom.resname != "SOL" {
                Some(p.atom.resindex)
            } else {
                None
            }
        })?;

        par.par_iter().try_for_each(|sel| {
            println!("{} {}", sel.len(), sel.first_atom().resname);
            if sel.first_atom().resname == "ALA" {
                let subsel = sel.select("name CA")?;
                let pos = subsel.first_pos();
                println!("{}", pos);
            }
            Ok::<_, SelectionError>(())
        })?;

        // Add serial selection
        let ca = sys.select("name CA")?;
        println!("#ca: {}", ca.len());

        //Iter serial views
        let serials: Vec<Sel> = par.iter().collect();
        println!("serial #5: {}", serials[5].first_atom().resname);

        Ok(())
    }

    #[test]
    fn test1_set_state() -> anyhow::Result<()> {
        let (top, st1) = FileHandler::open("tests/albumin.pdb")?.read()?;
        let st2 = FileHandler::open("tests/albumin.pdb")?
            .read_state()?
            .unwrap();
        st2.set_time(100.0);

        let sys = System::new(top, st1)?;
        let mut sel1 = sys.select("name CA")?;
        let sel2 = sys.select("name CB")?;

        println!(
            "Before set_state(): sys: {} sel1: {} sel2: {}",
            sys.get_time(),
            sel1.get_time(),
            sel2.get_time()
        );
        sel1.set_state(Arc::new(st2))?;
        println!(
            "After set_state(): sys: {} sel1: {} sel2: {}",
            sys.get_time(),
            sel1.get_time(),
            sel2.get_time()
        );
        Ok(())
    }
}
