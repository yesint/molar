use crate::prelude::*;
use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator};
use std::marker::PhantomData;
use sync_unsafe_cell::SyncUnsafeCell;
use triomphe::Arc;

struct SelIndex(SyncUnsafeCell<SVec>);

impl SelIndex {
    fn get(&self) -> &SVec {
        unsafe { &*self.0.get() }
    }

    fn get_mut(&self) -> &mut SVec {
        unsafe { &mut *self.0.get() }
    }

    fn from_sorted_vec(v: Vec<usize>) -> Self {
        unsafe { Self(SyncUnsafeCell::new(SVec::from_sorted(v))) }
    }

    fn from_svec(v: SVec) -> Self {
        Self(SyncUnsafeCell::new(v))
    }
}

// Private module containing internal methods implemented by all selection types
mod private {
    use crate::core::{
        selection::sel2::SelIndex, IndexProvider, LenProvider, SVec, State, Topology,
    };
    use triomphe::Arc;

    // Trait for things that can provide refs to Topology and State
    // (for implementing analysis traits, both Systems and Selections)
    pub trait HasTopState: LenProvider + IndexProvider {
        fn get_topology(&self) -> &Topology;
        fn get_state(&self) -> &State;
    }

    // Trait for things that allow creating selections from them
    //(System, Sel, SelPar)
    pub trait AllowsSelecting {
        fn topology_arc(&self) -> &Arc<Topology>;
        fn state_arc(&self) -> &Arc<State>;
        fn index_slice(&self) -> Option<&[usize]>;
    }

    // Internal trait for selections
    pub trait SelectionPrivate: AllowsSelecting {
        fn index_mut(&self) -> &mut SVec;

        #[allow(private_interfaces)]
        fn new_sel(topology: Arc<Topology>, state: Arc<State>, index: Arc<SelIndex>) -> Self
        where
            Self: Sized;
    }
}

//----------------------------------------------------------

/// Trait for things that supports selecting atoms (Systems and Selections)
pub trait Selectable: private::AllowsSelecting + private::HasTopState {
    type DerivedSel: private::SelectionPrivate;

    // Sub-selection
    fn select(&self, def: impl SelectionDef) -> Result<Self::DerivedSel, SelectionError>
    where
        Self: Sized,
    {
        use private::SelectionPrivate;

        let ind = def.into_sel_index(self.topology_arc(), self.state_arc(), self.index_slice())?;
        Ok(Self::DerivedSel::new_sel(
            Arc::clone(&self.topology_arc()),
            Arc::clone(&self.state_arc()),
            Arc::new(SelIndex::from_svec(ind)),
        ))
    }

    // Parallel split
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

    // Serial split
    fn split_iter<RT, F, S, SO>(&self, func: F) -> impl Iterator<Item = Self::DerivedSel>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(Particle) -> Option<RT>,
        Self: Sized,
        Self::DerivedSel: Selectable,
    {
        split_iter_as(self, func)
    }

    fn set_topology(
        &mut self,
        topology: impl Into<Arc<Topology>>,
    ) -> Result<Arc<Topology>, SelectionError> {
        let topology: Arc<Topology> = topology.into();
        if !self.topology_arc().interchangeable(&topology) {
            return Err(SelectionError::IncompatibleState);
        }
        let p = self.topology_arc() as *const Arc<Topology> as *mut Arc<Topology>;
        Ok(unsafe { std::ptr::replace(p, topology) })
    }

    /// Sets new [State] in [Sel].
    /// This is "shallow" operation, which doesn't affect other
    /// [Sel]s that point to the same [State].
    fn set_state(&self, state: impl Into<Arc<State>>) -> Result<Arc<State>, SelectionError> {
        let state: Arc<State> = state.into();
        if !self.state_arc().interchangeable(&state) {
            return Err(SelectionError::IncompatibleState);
        }
        let p = self.state_arc() as *const Arc<State> as *mut Arc<State>;
        Ok(unsafe { std::ptr::replace(p, state) })
    }

    //==============================================
    // Setting same properties for selections
    //==============================================

    fn set_same_name(&self, val: &str)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.name = val.into();
        }
    }

    fn set_same_resname(&self, val: &str)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.resname = val.into();
        }
    }

    fn set_same_resid(&self, val: i32)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.resid = val;
        }
    }

    fn set_same_chain(&self, val: char)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.chain = val;
        }
    }

    fn set_same_mass(&self, val: f32)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.mass = val;
        }
    }

    fn set_same_bfactor(&self, val: f32)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.bfactor = val;
        }
    }

    fn set_same_occupancy(&self, val: f32)
    where
        Self: Sized,
    {
        for a in self.iter_atoms_mut() {
            a.occupancy = val;
        }
    }
}

// Public trait for selections
pub trait Selection: private::SelectionPrivate {
    fn new_view<S: private::SelectionPrivate>(&self) -> S;

    //==============================================
    // Logic on selections that modify existing ones
    //==============================================

    fn add(&self, other: impl SelectionDef) -> Result<(), SelectionError> {
        self.index_mut().extend(
            other
                .into_sel_index(self.topology_arc(), self.state_arc(), None)?
                .iter()
                .cloned(),
        );
        Ok(())
    }

    // pub fn remove_global(&mut self, other: impl SelectionDef) -> Result<(), SelectionError> {
    //     let lhs = rustc_hash::FxHashSet::<usize>::from_iter(self.iter_index());
    //     let rhs = rustc_hash::FxHashSet::<usize>::from_iter(
    //         other
    //             .into_sel_index(&self.topology, &self.state, None)?
    //             .iter()
    //             .cloned(),
    //     );
    //     let v: Vec<usize> = lhs.difference(&rhs).cloned().collect();
    //     self.index_storage = SortedSet::from_unsorted(v);
    //     Ok(())
    // }

    // pub fn remove_local(&mut self, other: impl IntoIterator<Item = usize>) {
    //     let lhs = rustc_hash::FxHashSet::<usize>::from_iter(0..self.len());
    //     let rhs = rustc_hash::FxHashSet::<usize>::from_iter(other);
    //     let v: Vec<usize> = lhs
    //         .difference(&rhs)
    //         .map(|i| self.index_storage[*i])
    //         .collect();
    //     self.index_storage = SortedSet::from_unsorted(v);
    // }

    // pub fn invert(&mut self) {
    //     let lhs = rustc_hash::FxHashSet::<usize>::from_iter(0..self.len());
    //     let rhs = rustc_hash::FxHashSet::<usize>::from_iter(self.iter_index());
    //     let v: Vec<usize> = lhs.difference(&rhs).cloned().collect();
    //     self.index_storage = SortedSet::from_unsorted(v);
    // }
}

impl<T: Selectable> private::HasTopState for T {
    fn get_topology(&self) -> &Topology {
        self.topology_arc()
    }

    fn get_state(&self) -> &State {
        self.state_arc()
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
                        Arc::clone(sel.topology_arc()),
                        Arc::clone(sel.state_arc()),
                        Arc::new(SelIndex::from_sorted_vec(index)),
                    ));
                }
            }
            // Next particle
            cur += 1;
        }

        // Return any remaining index as last selection
        if !index.is_empty() {
            return Some(SO::new_sel(
                Arc::clone(sel.topology_arc()),
                Arc::clone(sel.state_arc()),
                Arc::new(SelIndex::from_sorted_vec(index)),
            ));
        }

        // If we are here stop iterating
        None
    };

    std::iter::from_fn(next_fn)
}

/// Marker trait for serial selections
pub trait SerialSelection: private::SelectionPrivate {}

//-----------------------------------------

macro_rules! impl_selection {
    ($t:ty, $d:ident) => {
        impl Selectable for $t {
            type DerivedSel = $d;
        }

        impl private::AllowsSelecting for $t {
            fn index_slice(&self) -> Option<&[usize]> {
                Some(self.index_storage.get())
            }

            fn state_arc(&self) -> &Arc<State> {
                &self.state
            }

            fn topology_arc(&self) -> &Arc<Topology> {
                &self.topology
            }
        }

        impl LenProvider for $t {
            fn len(&self) -> usize {
                self.index_storage.get().len()
            }
        }

        impl IndexProvider for $t {
            fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
                self.index_storage.get().iter().cloned()
            }

            unsafe fn get_index_unchecked(&self, i: usize) -> usize {
                *self.index_storage.get().get_unchecked(i)
            }
        }
    };
}

//-----------------------------------------
// Primary serial selections
pub struct SelSerial {
    topology: Arc<Topology>,
    state: Arc<State>,
    index_storage: Arc<SelIndex>,
    _phantom: PhantomData<*const ()>,
}

impl private::SelectionPrivate for SelSerial {
    fn index_mut(&self) -> &mut SVec {
        self.index_storage.get_mut()
    }

    #[allow(private_interfaces)]
    fn new_sel(topology: Arc<Topology>, state: Arc<State>, index: Arc<SelIndex>) -> Self
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

impl Selection for SelSerial {
    fn new_view<S: private::SelectionPrivate>(&self) -> S {
        S::new_sel(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            Arc::clone(&self.index_storage),
        )
    }
}

impl_selection!(SelSerial, SelSerial);

impl SerialSelection for SelSerial {}

//-------------------------------------------------------

pub struct System {
    topology: Arc<Topology>,
    state: Arc<State>,
    _phantom: PhantomData<*const ()>,
}

impl private::AllowsSelecting for System {
    fn index_slice(&self) -> Option<&[usize]> {
        None
    }

    fn state_arc(&self) -> &Arc<State> {
        &self.state
    }

    fn topology_arc(&self) -> &Arc<Topology> {
        &self.topology
    }
}

impl Selectable for System {
    type DerivedSel = SelSerial;
}

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
    pub fn new(top: Topology, st: State) -> Result<Self, SelectionError> {
        check_topology_state_sizes(&top, &st)?;
        Ok(Self {
            topology: Arc::new(top),
            state: Arc::new(st),
            _phantom: Default::default(),
        })
    }

    pub fn select_all(&self) -> Result<SelSerial, SelectionError> {
        self.select(0..self.len())
    }

    pub fn append(&self, data: &(impl PosIterProvider + AtomIterProvider)) -> SelSerial {
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

    pub fn append_atoms(
        &self,
        atoms: impl Iterator<Item = Atom>,
        coords: impl Iterator<Item = Pos>,
    ) -> SelSerial {
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
        Ok(())
    }
}

//---------------------------------------------------

pub struct SelPar {
    topology: Arc<Topology>,
    state: Arc<State>,
    index_storage: Arc<SelIndex>,
}

impl private::SelectionPrivate for SelPar {
    fn index_mut(&self) -> &mut SVec {
        self.index_storage.get_mut()
    }

    #[allow(private_interfaces)]
    fn new_sel(topology: Arc<Topology>, state: Arc<State>, index: Arc<SelIndex>) -> Self
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

impl Selection for SelPar {
    fn new_view<S: private::SelectionPrivate>(&self) -> S {
        S::new_sel(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            Arc::clone(&self.index_storage),
        )
    }
}

impl_selection!(SelPar, SelSerial);

//--------------------------------------------------------------------
pub struct ParSplit {
    selections: Vec<SelPar>,
    _phantom: PhantomData<*const ()>,
}

impl ParSplit {
    /// Returns parallel iterator over stored parallel selections.
    pub fn par_iter(&mut self) -> rayon::slice::Iter<'_, SelPar> {
        self.selections.par_iter()
    }

    /// Returns parallel mutable iterator over stored parallel selections.
    pub fn par_iter_mut(&mut self) -> rayon::slice::IterMut<'_, SelPar> {
        self.selections.par_iter_mut()
    }

    pub fn iter<S: SerialSelection>(&self) -> impl Iterator<Item = S> + '_ {
        self.selections.iter().map(|sel| sel.new_view())
    }
}

//═══════════════════════════════════════════════════════════
//  Blanket trait implementations for Selections and System
//═══════════════════════════════════════════════════════════

//██████  IO traits

impl<T: private::HasTopState> WritableToFile for T {}

impl<T: private::HasTopState> TopologyIoProvider for T {}
impl<T: private::HasTopState> StateIoProvider for T {}

impl<T: private::HasTopState> TimeProvider for T {
    fn get_time(&self) -> f32 {
        self.get_state().get_time()
    }
}

//██████  Immutable analysis traits

impl<T: private::HasTopState> BoxProvider for T {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.get_state().get_box()
    }
}

impl<T: private::HasTopState> PosIterProvider for T {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        unsafe {
            self.iter_index()
                .map(|i| self.get_state().get_pos_unchecked(i))
        }
    }
}

impl<T: private::HasTopState> AtomIterProvider for T {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        unsafe {
            self.iter_index()
                .map(|i| self.get_topology().get_atom_unchecked(i))
        }
    }
}

impl<T: private::HasTopState> AtomIterMutProvider for T {
    fn iter_atoms_mut(&self) -> impl AtomMutIterator<'_> {
        unsafe {
            self.iter_index()
                .map(|i| self.get_topology().get_atom_mut_unchecked(i))
        }
    }
}

impl<T: private::HasTopState> MassIterProvider for T {
    fn iter_masses(&self) -> impl Iterator<Item = f32> {
        unsafe {
            self.iter_index()
                .map(|i| self.get_topology().get_atom_unchecked(i).mass)
        }
    }
}

impl<T: private::HasTopState> RandomPosProvider for T {
    unsafe fn get_pos_unchecked(&self, i: usize) -> &Pos {
        let ind = self.get_index_unchecked(i);
        self.get_state().get_pos_unchecked(ind)
    }
}

impl<T: private::HasTopState> RandomAtomProvider for T {
    unsafe fn get_atom_unchecked(&self, i: usize) -> &Atom {
        let ind = self.get_index_unchecked(i);
        self.get_topology().get_atom_unchecked(ind)
    }
}

impl<T: private::HasTopState> RandomAtomMutProvider for T {
    unsafe fn get_atom_mut_unchecked(&self, i: usize) -> &mut Atom {
        let ind = self.get_index_unchecked(i);
        self.get_topology().get_atom_mut_unchecked(ind)
    }
}

impl<T: private::HasTopState> MoleculesProvider for T {
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

impl<T: private::HasTopState> BondsProvider for T {
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

impl<T: private::HasTopState> RandomParticleProvider for T {
    unsafe fn get_particle_unchecked(&self, i: usize) -> Particle<'_> {
        let ind = self.get_index_unchecked(i);
        Particle {
            id: ind,
            atom: self.get_topology().get_atom_unchecked(ind),
            pos: self.get_state().get_pos_unchecked(ind),
        }
    }
}

impl<T: private::HasTopState> ParticleIterProvider for T {
    fn iter_particle(&self) -> impl Iterator<Item = Particle<'_>> {
        unsafe{ self.iter_index().map(|i| self.get_particle_unchecked(i)) }
    }
}

//██████  Measure traits

impl<T: private::HasTopState> MeasurePos for T {}
impl<T: private::HasTopState> MeasurePeriodic for T {}
impl<T: private::HasTopState> MeasureMasses for T {}
impl<T: private::HasTopState> MeasureRandomAccess for T {}

//██████  Mutable analysis traits

impl BoxMutProvider for System {
    fn get_box_mut(&self) -> Option<&mut PeriodicBox> {
        use private::HasTopState;
        self.get_state().get_box_mut()
    }
}

impl<T: private::HasTopState> PosIterMutProvider for T {
    fn iter_pos_mut(&self) -> impl PosMutIterator<'_> {
        unsafe {
            self.iter_index()
                .map(|i| self.get_state().get_pos_mut_unchecked(i))
        }
    }
}

impl<T: private::HasTopState> RandomPosMutProvider for T {
    unsafe fn get_pos_mut_unchecked(&self, i: usize) -> &mut Pos {
        let ind = self.get_index_unchecked(i);
        self.get_state().get_pos_mut_unchecked(ind)
    }
}

impl<T: private::HasTopState> RandomParticleMutProvider for T {
    unsafe fn get_particle_mut_unchecked(&self, i: usize) -> ParticleMut {
        let ind = self.get_index_unchecked(i);
        ParticleMut {
            id: ind,
            atom: self.get_topology().get_atom_mut_unchecked(ind),
            pos: self.get_state().get_pos_mut_unchecked(ind),
        }
    }
}

impl<T: private::HasTopState> ParticleIterMutProvider for T {
    fn iter_particle_mut(&self) -> impl Iterator<Item = ParticleMut<'_>> {
        unsafe{ self.iter_index().map(|i| self.get_particle_mut_unchecked(i)) }
    }
}

//██████  Modify traits

impl<T: private::HasTopState> ModifyPos for T {}
impl<T: private::HasTopState> ModifyPeriodic for T {}
impl<T: private::HasTopState> ModifyRandomAccess for T {}

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
        let serials: Vec<SelSerial> = par.iter().collect();
        println!("serial #5: {}", serials[5].first_atom().resname);

        Ok(())
    }
}
