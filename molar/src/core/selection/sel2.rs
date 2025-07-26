use crate::{core::selection::sel2::private::SelectionPrivate, prelude::*};
use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator};
use std::marker::PhantomData;
use triomphe::Arc;

// Private module containing internal methods implemented by all selection types
mod private {
    use crate::core::{IndexProvider, LenProvider, SVec, State, Topology};
    use triomphe::Arc;

    pub trait TopSt: LenProvider + IndexProvider {
        fn get_topology(&self) -> &Topology;
        fn get_state(&self) -> &State;
    }

    pub trait SystemPrivate {
        fn topology_arc(&self) -> &Arc<Topology>;
        fn state_arc(&self) -> &Arc<State>;
        fn index_slice(&self) -> Option<&[usize]>;
    }

    pub trait SelectionPrivate: SystemPrivate {
        fn from_arcs(topology: Arc<Topology>, state: Arc<State>, index: Arc<SVec>) -> Self
        where
            Self: Sized;
    }
}

//----------------------------------------------------------

// Public trait for selections
pub trait Selection: private::SystemPrivate + LenProvider + IndexProvider {
    type DerivedSel: private::SelectionPrivate;

    // Sub-selection
    fn select(&self, def: impl SelectionDef) -> Result<Self::DerivedSel, SelectionError>
    where
        Self: Sized,
    {
        use private::SelectionPrivate;

        let ind = def.into_sel_index(self.topology_arc(), self.state_arc(), self.index_slice())?;
        Ok(Self::DerivedSel::from_arcs(
            Arc::clone(&self.topology_arc()),
            Arc::clone(&self.state_arc()),
            Arc::new(ind),
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
        Self::DerivedSel: Selection,
    {
        split_iter_as(self, func)
    }
}

impl<T: Selection> private::TopSt for T {
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
    S: Selection,
    SO: SelectionPrivate,
{
    // Check validity of index in case if called from Builder selection
    // It will check bounds internally and panic immediately if not valid
    // For other selection types ir does nothing
    //let _ = sel.index_arc();

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
                    return Some(SO::from_arcs(
                        Arc::clone(sel.topology_arc()),
                        Arc::clone(sel.state_arc()),
                        Arc::new(unsafe { SVec::from_sorted(index) }),
                    ));
                }
            }
            // Next particle
            cur += 1;
        }

        // Return any remaining index as last selection
        if !index.is_empty() {
            return Some(SO::from_arcs(
                Arc::clone(sel.topology_arc()),
                Arc::clone(sel.state_arc()),
                Arc::new(unsafe { SVec::from_sorted(index) }),
            ));
        }

        // If we are here stop iterating
        None
    };

    std::iter::from_fn(next_fn)
}

// Marker trait for serial selections
pub trait SerialSelection: private::SelectionPrivate {}

//-----------------------------------------

macro_rules! impl_selection {
    ($t:ty, $d:ident) => {
        impl Selection for $t {
            type DerivedSel = $d;
        }

        impl private::SystemPrivate for $t {
            fn index_slice(&self) -> Option<&[usize]> {
                Some(&self.index_storage)
            }

            fn state_arc(&self) -> &Arc<State> {
                &self.state
            }

            fn topology_arc(&self) -> &Arc<Topology> {
                &self.topology
            }
        }

        impl private::SelectionPrivate for $t {
            fn from_arcs(topology: Arc<Topology>, state: Arc<State>, index: Arc<SVec>) -> Self
            where
                Self: Sized,
            {
                Self {
                    topology,
                    state,
                    index_storage: index,
                    ..Default::default()
                }
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

        impl $t {
            fn new_view<S: SerialSelection>(&self) -> S {
                use private::SystemPrivate;
                S::from_arcs(
                    Arc::clone(self.topology_arc()),
                    Arc::clone(self.state_arc()),
                    Arc::clone(&self.index_storage),
                )
            }
        }
    };
}

//-----------------------------------------
// Primary serial selections
#[derive(Default)]
pub struct SelSerial {
    topology: Arc<Topology>,
    state: Arc<State>,
    index_storage: Arc<SVec>,
    _phantom: PhantomData<*const ()>,
}

impl_selection!(SelSerial, SelSerial);

impl SerialSelection for SelSerial {}

impl SelSerial {
    pub fn new_all(top: Topology, st: State) -> Result<Self, SelectionError> {
        check_topology_state_sizes(&top, &st)?;
        let n = top.len();
        Ok(Self {
            topology: Arc::new(top),
            state: Arc::new(st),
            index_storage: Arc::new(unsafe { SVec::from_sorted((0..n).collect::<Vec<_>>()) }),
            _phantom: Default::default(),
        })
    }
}

//--------------------------------------------

#[derive(Default)]
pub struct SelBuilder {
    topology: Arc<Topology>,
    state: Arc<State>,
    index_storage: Arc<SVec>,
    _phantom: PhantomData<*const ()>,
}

//TODO: We don't use macro because custom index_arc() is needed
impl_selection!(SelBuilder, SelBuilder);
//-------------------------------------------------------

pub struct Builder {
    topology: Arc<Topology>,
    state: Arc<State>,
    _phantom: PhantomData<*const ()>,
}

impl private::SystemPrivate for Builder {
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

impl Selection for Builder {
    type DerivedSel = SelBuilder;
}

impl Builder {
    pub fn new(top: Topology, st: State) -> Result<Self, SelectionError> {
        check_topology_state_sizes(&top, &st)?;
        Ok(Self {
            topology: Arc::new(top),
            state: Arc::new(st),
            _phantom: Default::default(),
        })
    }

    pub fn select(&self, def: impl SelectionDef) -> Result<SelBuilder, SelectionError> {
        let ind = def.into_sel_index(&self.topology, &self.state, None)?;
        Ok(SelBuilder {
            topology: Arc::clone(&self.topology),
            state: Arc::clone(&self.state),
            index_storage: Arc::new(ind),
            _phantom: Default::default(),
        })
    }

    pub fn select_all(&self) -> Result<SelBuilder, SelectionError> {
        let ind = (0..self.len()).into_sel_index(&self.topology, &self.state, None)?;
        Ok(SelBuilder {
            topology: Arc::clone(&self.topology),
            state: Arc::clone(&self.state),
            index_storage: Arc::new(ind),
            _phantom: Default::default(),
        })
    }

    pub fn append(&self, data: &(impl PosIterProvider + AtomIterProvider)) -> SelBuilder {
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
    ) -> SelBuilder {
        let first_added_index = self.len();
        self.topology.get_storage_mut().add_atoms(atoms);
        self.state.get_storage_mut().add_coords(coords);
        let last_added_index = self.len();
        self.select(first_added_index..last_added_index).unwrap()
    }

    pub fn remove(&self, to_remove: &impl IndexProvider) -> Result<(), SelectionError> {
        self.topology
            .get_storage_mut()
            .remove_atoms(to_remove.iter_index())?;
        self.state
            .get_storage_mut()
            .remove_coords(to_remove.iter_index())?;
        Ok(())
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

#[derive(Default)]
pub struct SelPar {
    topology: Arc<Topology>,
    state: Arc<State>,
    index_storage: Arc<SVec>,
}

impl_selection!(SelPar, SelSerial);

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

//========================================================
//  Blanket trait implementations for Selection
//========================================================

//══════════════════════════════════════════════
//███  IO traits
//══════════════════════════════════════════════

impl IndexProvider for Builder {
    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        0..self.len()
    }

    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        i
    }
}

impl<T: private::TopSt> WritableToFile for T {}

impl<T: private::TopSt> TopologyIoProvider for T {}

impl<T: private::TopSt> StateIoProvider for T {}

impl<T: private::TopSt> TimeProvider for T {
    fn get_time(&self) -> f32 {
        self.get_state().get_time()
    }
}

//══════════════════════════════════════════════
//███  Immutable analysis traits
//══════════════════════════════════════════════

impl<T: private::TopSt> BoxProvider for T {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.get_state().get_box()
    }
}

impl<T: private::TopSt> MeasurePeriodic for T {}

impl<T: private::TopSt> PosIterProvider for T {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        unsafe {
            self.iter_index()
                .map(|i| self.get_state().get_pos_unchecked(i))
        }
    }
}

impl<T: private::TopSt> MeasurePos for T {}

impl<T: private::TopSt> AtomIterProvider for T {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        unsafe {
            self.iter_index()
                .map(|i| self.get_topology().get_atom_unchecked(i))
        }
    }
}

impl<T: private::TopSt> AtomIterMutProvider for T {
    fn iter_atoms_mut(&self) -> impl AtomMutIterator<'_> {
        unsafe {
            self.iter_index()
                .map(|i| self.get_topology().get_atom_mut_unchecked(i))
        }
    }
}

impl<T: private::TopSt> MassIterProvider for T {
    fn iter_masses(&self) -> impl Iterator<Item = f32> {
        unsafe {
            self.iter_index()
                .map(|i| self.get_topology().get_atom_unchecked(i).mass)
        }
    }
}

impl<T: private::TopSt> MeasureMasses for T {}

impl LenProvider for Builder {
    fn len(&self) -> usize {
        use private::TopSt;
        self.get_state().len()
    }
}

impl<T: private::TopSt> RandomPosProvider for T {
    unsafe fn get_pos_unchecked(&self, i: usize) -> &Pos {
        let ind = self.get_index_unchecked(i);
        self.get_state().get_pos_unchecked(ind)
    }
}

impl<T: private::TopSt> RandomAtomProvider for T {
    unsafe fn get_atom_unchecked(&self, i: usize) -> &Atom {
        let ind = self.get_index_unchecked(i);
        self.get_topology().get_atom_unchecked(ind)
    }
}

impl<T: private::TopSt> RandomAtomMutProvider for T {
    unsafe fn get_atom_mut_unchecked(&self, i: usize) -> &mut Atom {
        let ind = self.get_index_unchecked(i);
        self.get_topology().get_atom_mut_unchecked(ind)
    }
}

impl<T: private::TopSt> MeasureRandomAccess for T {}

impl<T: private::TopSt> MoleculesProvider for T {
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

impl<T: private::TopSt> BondsProvider for T {
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

impl<T: private::TopSt> RandomParticleProvider for T {
    unsafe fn get_particle_unchecked(&self, i: usize) -> Particle<'_> {
        let ind = self.get_index_unchecked(i);
        Particle {
            id: ind,
            atom: self.get_topology().get_atom_unchecked(ind),
            pos: self.get_state().get_pos_unchecked(ind),
        }
    }
}

//═══════════════════════════════════════════════════════════
//███  Mutable analysis traits (only for mutable selections)
//═══════════════════════════════════════════════════════════

impl BoxMutProvider for Builder {
    fn get_box_mut(&self) -> Option<&mut PeriodicBox> {
        use private::TopSt;
        self.get_state().get_box_mut()
    }
}

impl<T: private::TopSt> PosIterMutProvider for T {
    fn iter_pos_mut(&self) -> impl PosMutIterator<'_> {
        unsafe {
            self.iter_index()
                .map(|i| self.get_state().get_pos_mut_unchecked(i))
        }
    }
}

impl<T: private::TopSt> RandomPosMutProvider for T {
    unsafe fn get_pos_mut_unchecked(&self, i: usize) -> &mut Pos {
        let ind = self.get_index_unchecked(i);
        self.get_state().get_pos_mut_unchecked(ind)
    }
}

impl<T: private::TopSt> RandomParticleMutProvider for T {
    unsafe fn get_particle_mut_unchecked(&self, i: usize) -> ParticleMut {
        let ind = self.get_index_unchecked(i);
        ParticleMut {
            id: ind,
            atom: self.get_topology().get_atom_mut_unchecked(ind),
            pos: self.get_state().get_pos_mut_unchecked(ind),
        }
    }
}

impl<T: private::TopSt> ModifyPos for T {}
impl<T: private::TopSt> ModifyPeriodic for T {}
impl<T: private::TopSt> ModifyRandomAccess for T {}

//========================================================

#[allow(unused_imports)]
mod tests {
    use super::*;
    use crate::prelude::*;

    #[test]
    fn test1() -> anyhow::Result<()> {
        let (top, st) = FileHandler::open("tests/albumin.pdb")?.read()?;
        let sel = SelSerial::new_all(top, st)?;

        let mut par = sel.split_par(|p| {
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
        let ca = sel.select("name CA")?;
        println!("#ca: {}", ca.len());

        //Iter serial views
        let serials: Vec<SelSerial> = par.iter().collect();
        println!("serial #5: {}", serials[5].first_atom().resname);

        Ok(())
    }
}
