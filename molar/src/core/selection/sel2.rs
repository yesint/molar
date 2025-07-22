use crate::prelude::*;
use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator};
use std::marker::PhantomData;
use triomphe::Arc;

// Private module containing private internal methods implemented by
// all selection types
mod private {
    use crate::core::{SVec, State, Topology};
    use triomphe::Arc;

    pub trait SelectionBase {
        fn get_index(&self) -> &Arc<SVec>;
        fn get_topology(&self) -> &Arc<Topology>;
        fn get_state(&self) -> &Arc<State>;
        fn new_internal(topology: Arc<Topology>, state: Arc<State>, index: Arc<SVec>) -> Self
        where
            Self: Sized;
    }
}

// Blanket implementations for selections

impl<T: private::SelectionBase> LenProvider for T {
    fn len(&self) -> usize {
        self.get_index().len()
    }
}

impl<T: private::SelectionBase> IndexProvider for T {
    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize> {
        self.get_index().iter().cloned()
    }
}

impl<T: private::SelectionBase> ParticleIterProvider for T {
    fn iter_particle(&self) -> impl ExactSizeIterator<Item = Particle<'_>> {
        self.iter_index().map(|i| unsafe {
            Particle {
                atom: self.get_topology().nth_atom_unchecked(i),
                pos: self.get_state().nth_pos_unchecked(i),
                id: i,
            }
        })
    }
}

impl<T: private::SelectionBase> ParticleIterMutProvider for T {
    fn iter_particle_mut(&self) -> impl ExactSizeIterator<Item = ParticleMut<'_>> {
        self.iter_index().map(|i| unsafe {
            ParticleMut {
                atom: self.get_topology().nth_atom_mut_unchecked(i),
                pos: self.get_state().nth_pos_mut_unchecked(i),
                id: i,
            }
        })
    }
}

impl<T: private::SelectionBase> RandomParticleProvider for T {
    unsafe fn nth_particle_unchecked(&self, i: usize) -> Particle<'_> {
        Particle {
            atom: self.get_topology().nth_atom_unchecked(i),
            pos: self.get_state().nth_pos_unchecked(i),
            id: i,
        }
    }
}

impl<T: private::SelectionBase> RandomParticleMutProvider for T {
    unsafe fn nth_particle_mut_unchecked(&self, i: usize) -> ParticleMut {
        ParticleMut {
            atom: self.get_topology().nth_atom_mut_unchecked(i),
            pos: self.get_state().nth_pos_mut_unchecked(i),
            id: i,
        }
    }
}

impl<T: private::SelectionBase> RandomAtomProvider for T {
    fn num_atoms(&self) -> usize {
        self.get_index().len()
    }

    unsafe fn nth_atom_unchecked(&self, i: usize) -> &Atom {
        let ind = self.get_index().get_unchecked(i);
        self.get_topology().nth_atom_unchecked(*ind)
    }
}

impl<T: private::SelectionBase> RandomAtomMutProvider for T {
    unsafe fn nth_atom_mut_unchecked(&self, i: usize) -> &mut Atom {
        let ind = self.get_index().get_unchecked(i);
        self.get_topology().nth_atom_mut_unchecked(*ind)
    }
}

impl<T: private::SelectionBase> PosIterProvider for T {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.get_index().iter().map(|i| unsafe {
            self.get_state().nth_pos_unchecked(*i)
        })
    }
}

impl<T: private::SelectionBase> PosIterMutProvider for T {
    fn iter_pos_mut(&self) -> impl PosMutIterator<'_> {
        self.get_index().iter().map(|i| unsafe {
            self.get_state().nth_pos_mut_unchecked(*i)
        })
    }
}

impl<T: private::SelectionBase> BoxProvider for T {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.get_state().get_box()
    }
}

impl<T: private::SelectionBase> BoxMutProvider for T {
    fn get_box_mut(&self) -> Option<&mut PeriodicBox> {
        self.get_state().get_box_mut()
    }
}

impl<T: private::SelectionBase> AtomIterProvider for T {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        self.get_index().iter().map(|i| unsafe {
            self.get_topology().nth_atom_unchecked(*i)
        })
    }
}

impl<T: private::SelectionBase> RandomPosProvider for T {
    fn num_pos(&self) -> usize {
        self.get_state().len()
    }

    unsafe fn nth_pos_unchecked(&self, i: usize) -> &Pos {
        let ind = self.get_index().get_unchecked(i);
        self.get_state().nth_pos_unchecked(*ind)
    }
}

impl<T: private::SelectionBase> RandomPosMutProvider for T {
    unsafe fn nth_pos_mut_unchecked(&self, i: usize) -> &mut Pos {
        let ind = self.get_index().get_unchecked(i);
        self.get_state().nth_pos_mut_unchecked(*ind)
    }
}

impl<T: private::SelectionBase> ModifyPos for T {}
impl<T: private::SelectionBase> ModifyPeriodic for T {}
impl<T: private::SelectionBase> ModifyRandomAccess for T {}

//----------------------------------------------------------

pub trait Selection: private::SelectionBase {
    type SubSel: private::SelectionBase;

    fn select(&self, def: impl SelectionDef) -> Result<Self::SubSel, SelectionError>
    where
        Self: Sized,
    {
        use private::SelectionBase;

        let ind = def.into_sel_index(
            self.get_topology(),
            self.get_state(),
            self.get_index().as_slice().into(),
        )?;
        Ok(Self::SubSel::new_internal(
            Arc::clone(&self.get_topology()),
            Arc::clone(&self.get_state()),
            Arc::new(ind),
        ))
    }

    fn new_view<S: SerialSelection>(&self) -> S {
        S::new_internal(
            Arc::clone(self.get_topology()),
            Arc::clone(self.get_state()),
            Arc::clone(self.get_index()),
        )
    }

    fn split_parallel<F, R>(&self, func: F) -> Result<ParSplit, SelectionError>
    where
        F: Fn(Particle) -> Option<R>,
        R: Default + PartialOrd,
        Self: Sized,
    {
        let selections: Vec<SelPar> = split_iter(self,func).collect();

        if selections.is_empty() {
            return Err(SelectionError::EmptySplit);
        }

        Ok(ParSplit {
            selections,
            _phantom: Default::default(),
        })
    }
}

fn split_iter<RT, F, S, SO>(sel: &S, func: F) -> impl Iterator<Item = SO> + use<'_, RT, F, S, SO>
where
    RT: Default + std::cmp::PartialEq,
    F: Fn(Particle) -> Option<RT>,
    S: Selection,
    SO: Selection,
{
    let mut cur_val = RT::default();
    let mut cur = 0usize;

    let next_fn = move || {
        let mut index = Vec::<usize>::new();

        while cur < sel.len() {
            let p = unsafe { sel.nth_particle_unchecked(cur) };
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
                    return Some(SO::new_internal(
                        Arc::clone(sel.get_topology()),
                        Arc::clone(sel.get_state()),
                        Arc::new(unsafe{SVec::from_sorted(index)}),
                    ));
                }
            }
            // Next particle
            cur += 1;
        }

        // Return any remaining index as last selection
        if !index.is_empty() {
            return Some(SO::new_internal(
                Arc::clone(sel.get_topology()),
                Arc::clone(sel.get_state()),
                Arc::new(unsafe{SVec::from_sorted(index)}),
            ));
        }

        // If we are here stop iterating
        None
    };

    std::iter::from_fn(next_fn)
}

pub trait SerialSelection: private::SelectionBase {}
//-----------------------------------------

pub struct SelSerial {
    topology: Arc<Topology>,
    state: Arc<State>,
    index_storage: Arc<SVec>,
    _phantom: PhantomData<*const ()>,
}

impl private::SelectionBase for SelSerial {
    fn get_index(&self) -> &Arc<SVec> {
        &self.index_storage
    }

    fn get_state(&self) -> &Arc<State> {
        &self.state
    }

    fn get_topology(&self) -> &Arc<Topology> {
        &self.topology
    }

    fn new_internal(topology: Arc<Topology>, state: Arc<State>, index: Arc<SVec>) -> Self
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
    type SubSel = Self;
}

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
pub struct SelBuilder {
    topology: Arc<Topology>,
    state: Arc<State>,
    index_storage: Arc<SVec>,
    _phantom: PhantomData<*const ()>,
}

impl private::SelectionBase for SelBuilder {
    fn get_index(&self) -> &Arc<SVec> {
        if self.index_storage.len() > self.topology.len()
            || self.index_storage.len() > self.state.len()
        {
            panic!("selection index is invalidated",);
        }
        &self.index_storage
    }

    fn get_state(&self) -> &Arc<State> {
        &self.state
    }

    fn get_topology(&self) -> &Arc<Topology> {
        &self.topology
    }

    fn new_internal(topology: Arc<Topology>, state: Arc<State>, index: Arc<SVec>) -> Self
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

impl Selection for SelBuilder {
    type SubSel = Self;
}

impl SerialSelection for SelBuilder {}

//-------------------------------------------------------
pub struct Builder {
    topology: Arc<Topology>,
    state: Arc<State>,
    _phantom: PhantomData<*const ()>,
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
        let first_added_index = self.num_atoms();
        self.topology
            .get_storage_mut()
            .add_atoms(data.iter_atoms().cloned());
        self.state
            .get_storage_mut()
            .add_coords(data.iter_pos().cloned());
        let last_added_index = self.num_atoms();
        self.select(first_added_index..last_added_index).unwrap()
    }

    pub fn append_atoms(
        &self,
        atoms: impl Iterator<Item = Atom>,
        coords: impl Iterator<Item = Pos>,
    ) -> SelBuilder {
        let first_added_index = self.num_atoms();
        self.topology.get_storage_mut().add_atoms(atoms);
        self.state.get_storage_mut().add_coords(coords);
        let last_added_index = self.num_atoms();
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

macro_rules! trait_def {
    ( $tr:tt, $t:ty ) => {
        impl $tr for $t
    };
    ( $tr:tt ) => {
        impl <T: private::SelectionBase> $tr for T
    };
}

impl LenProvider for Builder {
    fn len(&self) -> usize {
        self.num_atoms()
    }
}

impl IndexProvider for Builder {
    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize> {
        0..self.topology.num_atoms()
    }
}

impl TopologyIoProvider for Builder {}

impl AtomIterProvider for Builder {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        self.topology.iter_atoms()
    }
}

impl RandomAtomProvider for Builder {
    fn num_atoms(&self) -> usize {
        self.topology.num_atoms()
    }

    unsafe fn nth_atom_unchecked(&self, i: usize) -> &Atom {
        self.topology.nth_atom_unchecked(i)
    }
}

impl StateIoProvider for Builder {}

impl TimeProvider for Builder {
    fn get_time(&self) -> f32 {
        self.state.get_time()
    }
}

impl PosIterProvider for Builder {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.state.iter_pos()
    }
}

impl BoxProvider for Builder {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.state.get_box()
    }
}

impl WritableToFile for Builder {}

impl ParticleIterProvider for Builder {
    fn iter_particle(&self) -> impl ExactSizeIterator<Item = Particle<'_>> {
        self.iter_index()
            .map(|i| unsafe { self.nth_particle_unchecked(i) })
    }
}

impl RandomParticleProvider for Builder {
    unsafe fn nth_particle_unchecked(&self, i: usize) -> Particle<'_> {
        Particle {
            id: i,
            atom: unsafe { self.topology.nth_atom_unchecked(i) },
            pos: unsafe { self.state.nth_pos_unchecked(i) },
        }
    }
}

impl RandomPosProvider for Builder {
    fn num_pos(&self) -> usize {
        self.state.num_pos()
    }

    unsafe fn nth_pos_unchecked(&self, i: usize) -> &Pos {
        self.state.nth_pos_unchecked(i)
    }
}

impl MoleculesProvider for Builder {
    fn num_molecules(&self) -> usize {
        self.topology.num_molecules()
    }

    fn iter_molecules(&self) -> impl Iterator<Item = &[usize; 2]> {
        self.topology.iter_molecules()
    }

    unsafe fn nth_molecule_unchecked(&self, i: usize) -> &[usize; 2] {
        self.topology.nth_molecule_unchecked(i)
    }
}

impl BondsProvider for Builder {
    fn num_bonds(&self) -> usize {
        self.topology.num_bonds()
    }

    fn iter_bonds(&self) -> impl Iterator<Item = &[usize; 2]> {
        self.topology.iter_bonds()
    }

    unsafe fn nth_bond_unchecked(&self, i: usize) -> &[usize; 2] {
        self.topology.nth_bond_unchecked(i)
    }
}

impl MassIterProvider for Builder {
    fn iter_masses(&self) -> impl ExactSizeIterator<Item = f32> {
        unsafe {
            self.iter_index()
                .map(|i| self.topology.nth_atom_unchecked(i).mass)
        }
    }
}

impl MeasureMasses for Builder {}
impl MeasureRandomAccess for Builder {}

//===============================================================
// Traits for read-write sources
//===============================================================

impl PosIterMutProvider for Builder {
    fn iter_pos_mut(&self) -> impl PosMutIterator<'_> {
        self.state.iter_pos_mut()
    }
}

impl RandomPosMutProvider for Builder {
    fn nth_pos_mut(&self, i: usize) -> Option<&mut Pos> {
        self.state.nth_pos_mut(i)
    }

    unsafe fn nth_pos_mut_unchecked(&self, i: usize) -> &mut Pos {
        self.state.nth_pos_mut_unchecked(i)
    }
}

impl AtomIterMutProvider for Builder {
    fn iter_atoms_mut(&self) -> impl AtomMutIterator<'_> {
        self.topology.iter_atoms_mut()
    }
}

impl RandomAtomMutProvider for Builder {
    unsafe fn nth_atom_mut_unchecked(&self, i: usize) -> &mut Atom {
        self.topology.nth_atom_mut_unchecked(i)
    }
}

impl ParticleIterMutProvider for Builder {
    fn iter_particle_mut(&self) -> impl ExactSizeIterator<Item = ParticleMut<'_>> {
        self.iter_index()
            .map(|i| unsafe { self.nth_particle_mut_unchecked(i) })
    }
}

impl RandomParticleMutProvider for Builder {
    unsafe fn nth_particle_mut_unchecked(&self, i: usize) -> ParticleMut {
        ParticleMut {
            id: i,
            atom: unsafe { self.topology.nth_atom_mut_unchecked(i) },
            pos: unsafe { self.state.nth_pos_mut_unchecked(i) },
        }
    }
}

impl ModifyPos for Builder {}
impl ModifyPeriodic for Builder {}
impl ModifyRandomAccess for Builder {}

//===============================================================
// For Serial kinds also BoxMut and TimeMut is implemnted
//===============================================================
impl BoxMutProvider for Builder {
    fn get_box_mut(&self) -> Option<&mut PeriodicBox> {
        self.state.get_box_mut()
    }
}

impl TimeMutProvider for Builder {
    fn set_time(&self, t: f32) {
        self.state.set_time(t);   
    }
}

//---------------------------------------------------

pub struct SelPar {
    topology: Arc<Topology>,
    state: Arc<State>,
    index_storage: Arc<SVec>,
}

impl private::SelectionBase for SelPar {
    fn get_index(&self) -> &Arc<SVec> {
        &self.index_storage
    }

    fn get_state(&self) -> &Arc<State> {
        &self.state
    }

    fn get_topology(&self) -> &Arc<Topology> {
        &self.topology
    }

    fn new_internal(topology: Arc<Topology>, state: Arc<State>, index: Arc<SVec>) -> Self
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
    type SubSel = SelSerial;
}

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

    pub fn iter_serial_views<S: SerialSelection>(&self) -> impl Iterator<Item = S> + '_ {
        self.selections.iter().map(|sel| sel.new_view())
    }
}

#[allow(unused_imports)]
mod tests {
    use super::*;
    use crate::prelude::*;

    #[test]
    fn test1() -> anyhow::Result<()> {
        let (top, st) = FileHandler::open("tests/albumin.pdb")?.read()?;
        let sel = SelSerial::new_all(top, st)?;

        let mut par = sel.split_parallel(|p| {
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

        // Iter serial views
        let serials: Vec<SelSerial> = par.iter_serial_views().collect();
        println!("serial #5: {}", serials[5].first_atom().resname);

        Ok(())
    }
}
