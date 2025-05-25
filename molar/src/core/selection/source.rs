use super::utils::*;
use crate::prelude::*;
use sorted_vec::SortedSet;
use std::{marker::PhantomData, path::Path};

//----------------------------------------
// Source of parallel selections
//----------------------------------------

/// Source for selections.
pub struct Source<K> {
    topology: Holder<Topology, K>,
    state: Holder<State, K>,
    _no_send: PhantomData<*const ()>,
    _kind: PhantomData<K>,
}

//=======================
// Constructors
//=======================

impl Source<()> {
    /// Create the [Source] producing mutable selections that may overlap accessible from a single thread.
    pub fn new_serial(
        topology: Holder<Topology, MutableSerial>,
        state: Holder<State, MutableSerial>,
    ) -> Result<Source<MutableSerial>, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        Ok(Source {
            topology: topology.into(),
            state: state.into(),
            _no_send: Default::default(),
            _kind: Default::default(),
        })
    }

    pub fn new_builder(
        topology: Holder<Topology, BuilderSerial>,
        state: Holder<State, BuilderSerial>,
    ) -> Result<Source<BuilderSerial>, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        Ok(Source {
            topology: topology.into(),
            state: state.into(),
            _no_send: Default::default(),
            _kind: Default::default(),
        })
    }

    pub fn empty_builder() -> Source<BuilderSerial> {
        Source {
            topology: Topology::default().into(),
            state: State::default().into(),
            _no_send: Default::default(),
            _kind: Default::default(),
        }
    }

    /// Creates a source of immutable parallel selections that _may_ overlap.
    pub fn new_parallel(
        topology: Holder<Topology, ImmutableParallel>,
        state: Holder<State, ImmutableParallel>,
    ) -> Result<Source<ImmutableParallel>, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        Ok(Source {
            topology: topology.into(),
            state: state.into(),
            _no_send: Default::default(),
            _kind: Default::default(),
        })
    }

    pub fn serial_from_file(
        fname: impl AsRef<Path>,
    ) -> Result<Source<MutableSerial>, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        Ok(Source::new_serial(top.into(), st.into())?)
    }

    pub fn builder_from_file(
        fname: impl AsRef<Path>,
    ) -> Result<Source<BuilderSerial>, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        Ok(Source::new_builder(top.into(), st.into())?)
    }

    pub fn parallel_from_file(
        fname: impl AsRef<Path>,
    ) -> Result<Source<ImmutableParallel>, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        Ok(Source::new_parallel(top.into(), st.into())?)
    }
}

//=======================
// Public API
//=======================

impl<K: UserCreatableKind> Source<K> {
    pub fn new(
        topology: impl Into<Holder<Topology, K>>,
        state: impl Into<Holder<State, K>>,
    ) -> Result<Source<K>, SelectionError> {
        let topology = topology.into();
        let state = state.into();
        check_topology_state_sizes(&topology, &state)?;
        Ok(Source {
            topology: topology,
            state: state,
            _no_send: Default::default(),
            _kind: Default::default(),
        })
    }

    pub fn from_file(fname: &str) -> Result<Source<K>, SelectionError> {
        let (top, st) = FileHandler::open(fname)?.read()?;
        Ok(Source::new(top, st)?)
    }

    /// Release and return [Topology] and [State]. 
    /// Fails if any selection is still pointing to them..
    pub fn release(self) -> Result<(Topology, State), SelectionError> {
        Ok((self.topology.release()?, self.state.release()?))
    }

    pub fn len(&self) -> usize {
        self.topology.num_atoms()
    }

    pub fn select(&self, def: impl SelectionDef) -> Result<Sel<K>, SelectionError> {
        Sel::from_holders_and_index(
            unsafe {self.topology.new_ref_with_kind()},
            unsafe {self.state.new_ref_with_kind()},
            def.into_sel_index(&self.topology, &self.state, None)?,
        )
    }

    pub unsafe fn select_vec_unchecked(&self, vec: Vec<usize>) -> Result<Sel<K>, SelectionError> {
        Sel::from_holders_and_index(
            self.topology.new_ref_with_kind(),
            self.state.new_ref_with_kind(),
            vec.into(),
        )
    }

    /// Creates selection of all
    pub fn select_all(&self) -> Result<Sel<K>, SelectionError> {
        Sel::from_holders_and_index(
            unsafe {self.topology.new_ref_with_kind()},
            unsafe {self.state.new_ref_with_kind()},
            unsafe { SortedSet::from_sorted((0..self.len()).collect()) },
        )
    }

    /// Sets new [State] in this [Source].
    /// This is "shallow" otheration, selections created earlier from this source
    /// are _not_ affected and still view an old state.
    ///
    /// New state should be compatible with the old one (have the same number of atoms). If not, the error is returned.
    ///
    /// Returns [Holder] with an old state, so it could be reused if needed.
    pub fn set_state(
        &mut self,
        state: impl Into<Holder<State, K>>,
    ) -> Result<Holder<State, K>, SelectionError> {
        let state: Holder<State, K> = state.into();
        if !self.state.interchangeable(&state) {
            return Err(SelectionError::IncompatibleState);
        }
        Ok(std::mem::replace(&mut self.state, state))
    }

    /// Sets new [Topology] in this source.
    /// This is "shallow" otheration, selections created earlier from this source
    /// are not affected and still view an old topology.
    /// New topology should be compatible with the old one (have the same number of atoms, bonds, molecules, etc.). If not, the error is returned.
    ///
    /// Returns [Holder] with an old topology, so it could be reused if needed.
    pub fn set_topology(
        &mut self,
        topology: impl Into<Holder<Topology, K>>,
    ) -> Result<Holder<Topology, K>, SelectionError> {
        let topology = topology.into();
        if !self.topology.interchangeable(&topology) {
            return Err(SelectionError::IncompatibleState);
        }
        Ok(std::mem::replace(&mut self.topology, topology))
    }

    pub fn get_topology(&self) -> Holder<Topology, K> {
        unsafe {self.topology.new_ref_with_kind()}
    }

    pub fn get_state(&self) -> Holder<State, K> {
        unsafe {self.state.new_ref_with_kind()}
    }
}

//=======================
// Builder API
//=======================

impl Source<BuilderSerial> {
    pub fn append(&self, data: &(impl PosIterProvider + AtomIterProvider)) -> Sel<BuilderSerial> {
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
    ) -> Sel<BuilderSerial> {
        let first_added_index = self.num_atoms();
        self.topology.get_storage_mut().add_atoms(atoms);
        self.state.get_storage_mut().add_coords(coords);
        let last_added_index = self.num_atoms();
        self.select(first_added_index..last_added_index).unwrap()
    }

    pub fn remove(&self, to_remove: &impl IndexProvider) -> Result<(), SelectionError> {
        // We are checking index validity inside remove methods
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

//=======================
// Trait impls
//=======================

//===============================================================
// Traits for all sources
//===============================================================

impl<K: SelectionKind> LenProvider for Source<K> {
    fn len(&self) -> usize {
        self.num_atoms()
    }
}

impl<K: SelectionKind> IndexProvider for Source<K> {
    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize> {
        0..self.topology.num_atoms()
    }
}

impl<K: SelectionKind> TopologyIoProvider for Source<K> {}

impl<K: SelectionKind> AtomIterProvider for Source<K> {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        self.topology.iter_atoms()
    }
}

impl<K: SelectionKind> RandomAtomProvider for Source<K> {
    fn num_atoms(&self) -> usize {
        self.topology.num_atoms()
    }

    unsafe fn nth_atom_unchecked(&self, i: usize) -> &Atom {
        self.topology.nth_atom_unchecked(i)
    }
}

impl<K: SelectionKind> StateIoProvider for Source<K> {
    fn get_time(&self) -> f32 {
        self.state.get_time()
    }
}

impl<K: SelectionKind> PosIterProvider for Source<K> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.state.iter_pos()
    }
}

impl<K: SelectionKind> BoxProvider for Source<K> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.state.get_box()
    }
}

impl<K: SelectionKind> WritableToFile for Source<K> {}

impl<K: SelectionKind> ParticleIterProvider for Source<K> {
    fn iter_particle(&self) -> impl ExactSizeIterator<Item = Particle<'_>> {
        self.iter_index()
            .map(|i| unsafe { self.nth_particle_unchecked(i) })
    }
}

impl<K: SelectionKind> RandomParticleProvider for Source<K> {
    unsafe fn nth_particle_unchecked(&self, i: usize) -> Particle<'_> {
        Particle {
            id: i,
            atom: unsafe { self.topology.nth_atom_unchecked(i) },
            pos: unsafe { self.state.nth_pos_unchecked(i) },
        }
    }
}

impl<K: SelectionKind> RandomPosProvider for Source<K> {
    fn num_pos(&self) -> usize {
        self.state.num_pos()
    }

    unsafe fn nth_pos_unchecked(&self, i: usize) -> &Pos {
        self.state.nth_pos_unchecked(i)
    }
}

impl<K: SelectionKind> MoleculesProvider for Source<K> {
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

impl<K: SelectionKind> BondsProvider for Source<K> {
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

impl<K: SelectionKind> MassIterProvider for Source<K> {
    fn iter_masses(&self) -> impl ExactSizeIterator<Item = f32> {
        unsafe {
            self.iter_index()
                .map(|i| self.topology.nth_atom_unchecked(i).mass)
        }
    }
}

impl<K: SelectionKind> MeasureMasses for Source<K> {}
impl<K: SelectionKind> MeasureRandomAccess for Source<K> {}

//===============================================================
// Traits for read-write sources
//===============================================================

impl<K: MutableKind> PosIterMutProvider for Source<K> {
    fn iter_pos_mut(&self) -> impl PosMutIterator<'_> {
        self.state.iter_pos_mut()
    }
}

impl<K: MutableKind> RandomPosMutProvider for Source<K> {
    fn nth_pos_mut(&self, i: usize) -> Option<&mut Pos> {
        self.state.nth_pos_mut(i)
    }

    unsafe fn nth_pos_mut_unchecked(&self, i: usize) -> &mut Pos {
        self.state.nth_pos_mut_unchecked(i)
    }
}

impl<K: MutableKind> AtomsIterMutProvider for Source<K> {
    fn iter_atoms_mut(&self) -> impl AtomMutIterator<'_> {
        self.topology.iter_atoms_mut()
    }
}

impl<K: MutableKind> RandomAtomMutProvider for Source<K> {
    unsafe fn nth_atom_mut_unchecked(&self, i: usize) -> &mut Atom {
        self.topology.nth_atom_mut_unchecked(i)
    }
}

impl<K: MutableKind> ParticleIterMutProvider for Source<K> {
    fn iter_particle_mut(&self) -> impl ExactSizeIterator<Item = ParticleMut<'_>> {
        self.iter_index()
            .map(|i| unsafe { self.nth_particle_mut_unchecked(i) })
    }
}

impl<K: MutableKind> RandomParticleMutProvider for Source<K> {
    unsafe fn nth_particle_mut_unchecked(&self, i: usize) -> ParticleMut {
        ParticleMut {
            id: i,
            atom: unsafe { self.topology.nth_atom_mut_unchecked(i) },
            pos: unsafe { self.state.nth_pos_mut_unchecked(i) },
        }
    }
}

impl<K: MutableKind> ModifyPos for Source<K> {}
impl<K: MutableKind> ModifyPeriodic for Source<K> {}
impl<K: MutableKind> ModifyRandomAccess for Source<K> {}

//===============================================================
// For Builder also BoxMut is implemnted
//===============================================================
impl BoxMutProvider for Source<BuilderSerial> {
    fn get_box_mut(&self) -> Option<&mut PeriodicBox> {
        self.state.get_box_mut()
    }
}
