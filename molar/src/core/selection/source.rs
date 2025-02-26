use super::utils::*;
use crate::prelude::*;
use sorted_vec::SortedSet;
use std::marker::PhantomData;

//----------------------------------------
// Source of parallel selections
//----------------------------------------

/// Source for selections, which are processed in parallel.
///
/// Selections are stored inside the source and the user can only access them directly by reference form
/// the same thread where [SourceParaltopologylel] was created (an attempt to send them to other thread won't compile).
/// In order to process selections in parallel user calls `map_par` method to run an arbitrary closure on each stored
/// selection in separate threads.
///
/// It is safe to change the [State] contained inside [SourceParallel] by calling [set_state](SourceParallel::set_state)
/// because it is guaranteed that no other threads are accessing stored selections at the same time.
///
/// # Example 2: immutable overlapping selections, using par_iter()
/// ```
/// # use molar::prelude::*;
/// # use anyhow::Result;
/// # use rayon::prelude::*;
/// # let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
/// let mut src = Source::new_parallel(top.into(), st.into())?;
/// // Overlapping selections
/// let mut sels = vec![];
/// sels.push( src.select_iter(0..10)? );
/// sels.push( src.select_iter(5..15)? );
/// // Process them using par_iter explicitly
/// let res = sels.par_iter().map(|sel| {
///    Ok(sel.center_of_mass()?)
/// }).collect::<Result<Vec<_>>>()?;
/// println!("{:?}",res);
/// #  Ok::<(), anyhow::Error>(())
/// ```
//------------------------------------------------------------------

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

    pub fn parallel_from_file(fname: &str) -> Result<Source<ImmutableParallel>, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        Ok(Source::new_parallel(top.into(), st.into())?)
    }

    pub fn serial_from_file(fname: &str) -> Result<Source<MutableSerial>, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        Ok(Source::new_serial(top.into(), st.into())?)
    }

    pub fn builder_from_file(fname: &str) -> Result<Source<BuilderSerial>, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        Ok(Source::new_builder(top.into(), st.into())?)
    }
}

//=======================
// Public API
//=======================

impl<K: UserCreatableKind> Source<K> {
    /// Release and return [Topology] and [State]. Fails if any selections created from this [Source] are still alive.
    pub fn release(self) -> Result<(Topology, State), SelectionError> {
        Ok((self.topology.release()?, self.state.release()?))
    }

    /// Creates new selection from an iterator of indexes. Indexes are bound checked, sorted and duplicates are removed.
    /// If any index is out of bounds the error is returned.
    pub fn select_iter(&self, iter: impl Iterator<Item = usize>) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_iter(iter, self.topology.num_atoms())?;
        self.new_sel(vec)
    }

    pub fn select_vec(&self, vec: Vec<usize>) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_vec(vec, self.topology.num_atoms())?;
        self.new_sel(vec)
    }

    pub unsafe fn select_vec_unchecked(&self, vec: Vec<usize>) -> Result<Sel<K>, SelectionError> {
        let vec = SortedSet::from_sorted(vec);
        self.new_sel(vec)
    }

    /// Creates selection of all
    pub fn select_all(&self) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_all(self.topology.num_atoms());
        self.new_sel(vec)
    }

    /// Creates new selection from a selection expression string. Selection expression is constructed internally but
    /// can't be reused. Consider using [add_expr](Self::add_expr) if you already have selection expression.
    pub fn select_str(&self, selstr: impl AsRef<str>) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_str(selstr.as_ref(), &self.topology, &self.state)?;
        self.new_sel(vec)
    }

    /// Creates new selection from an existing selection expression.
    pub fn select_expr(&self, expr: &mut SelectionExpr) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_expr(expr, &self.topology, &self.state)?;
        self.new_sel(vec)
    }

    /// Creates new selection from a range of indexes.
    /// If rangeis out of bounds the error is returned.
    pub fn select_range(&self, range: std::ops::Range<usize>) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_range(range, self.topology.num_atoms())?;
        self.new_sel(vec)
    }

    /// Sets new [State] in this source. All selections created from this source will automatically view
    /// the new state.
    ///
    /// New state should be compatible with the old one (have the same number of atoms). If not, the error is returned.
    ///
    /// Returns [Holder] with old state, so it could be reused if needed.
    pub fn set_state(
        &mut self,
        state: State,
    ) -> Result<State, SelectionError> {
        //let state: Holder<State, K>  = Holder::new(state);
        if !self.state.interchangeable(&state) {
            return Err(SelectionError::SetState);
        }
        let p = self.state.arc.as_ptr() as *mut State;
        //let mut p2 = state.arc.as_ptr() as *mut State;
        // We physically spap memory at these locations
        // this is cheap because coordinates are allocated on heap
        // and only pointers to allocations are swapped
        Ok(unsafe{ std::ptr::replace(p, state) })
    }

    /// Sets new [Topology] in this source. All selections created from this Source will automatically view
    /// the new topology.
    ///
    /// New topology should be compatible with the old one (have the same number of atoms, bonds, molecules, etc.). If not, the error is returned.
    ///
    /// Returns unique pointer to the old topology, so it could be reused if needed.
    ///
    /// # Safety
    /// If such change happens when selections from different threads are accessing the data
    /// inconsistent results may be produced, however this should never lead to issues with memory safety.
    pub fn set_topology(
        &mut self,
        topology: Holder<Topology, K>,
    ) -> Result<Holder<Topology, K>, SelectionError> {
        // Check if the states are compatible
        if !self.topology.interchangeable(&topology) {
            return Err(SelectionError::SetTopology);
        }

        Ok(std::mem::replace(&mut self.topology, topology))
    }

    pub fn get_topology(&self) -> Holder<Topology, K> {
        self.topology.clone_with_kind()
    }

    pub fn get_state(&self) -> Holder<State, K> {
        self.state.clone_with_kind()
    }

    fn new_sel(&self, index: SortedSet<usize>) -> Result<Sel<K>, SelectionError> {
        Sel::from_holders_and_index(
            self.topology.clone_with_kind(),
            self.state.clone_with_kind(),
            index,
        )
    }
}

//=======================
// Builder API
//=======================

impl Source<BuilderSerial> {
    pub fn append(&self, data: &(impl PosProvider + AtomProvider)) -> Sel<BuilderSerial> {
        let first_added_index = self.num_atoms();
        self.topology
            .get_storage_mut()
            .add_atoms(data.iter_atoms().cloned());
        self.state
            .get_storage_mut()
            .add_coords(data.iter_pos().cloned());
        let last_added_index = self.num_atoms();
        self.select_range(first_added_index..last_added_index)
            .unwrap()
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
        self.select_range(first_added_index..last_added_index)
            .unwrap()
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

/// Macros for implementing common traits for read-only sources
macro_rules! impl_read_only_source_traits {
    ( $t:ty ) => {
        impl TopologyProvider for $t {
            
        }
        
        impl AtomProvider for $t {
            fn iter_atoms(&self) -> impl AtomIterator<'_> {
                self.topology.iter_atoms()
            }
        }

        impl RandomAtomProvider for $t {
            fn num_atoms(&self) -> usize {
                self.topology.num_atoms()
            }

            unsafe fn nth_atom_unchecked(&self, i: usize) -> &Atom {
                self.topology.nth_atom_unchecked(i)
            }
        }
        
        impl StateProvider for $t {
            fn get_time(&self) -> f32 {
                self.state.get_time()
            }
        }
        
        impl PosProvider for $t {
            fn iter_pos(&self) -> impl PosIterator<'_> {
                self.state.iter_pos()
            }
        }
        
        impl BoxProvider for $t {
            fn get_box(&self) -> Option<&PeriodicBox> {
                self.state.get_box()
            }
        }
        
        impl WritableToFile for $t {}

        impl ParticleProvider for $t {
            fn iter_particle(&self) -> impl ExactSizeIterator<Item = Particle<'_>> {
                self.iter_index().map(|i| Particle {
                    id: i,
                    atom: unsafe{self.topology.nth_atom_unchecked(i)},
                    pos: unsafe{self.state.nth_pos_unchecked(i)},
                })
            }
        }

        impl LenProvider for $t {
            fn len(&self) -> usize {
                self.state.num_coords()
            }
        }
        
        impl RandomPosProvider for $t {
            fn num_coords(&self) -> usize {
                self.state.num_coords()
            }

            unsafe fn nth_pos_unchecked(&self, i: usize) -> &Pos {
                self.state.nth_pos_unchecked(i)
            }
        }

        impl MoleculesProvider for $t {
            fn num_molecules(&self) -> usize {
                self.topology.num_molecules()
            }
        
            fn iter_molecules(&self) -> impl Iterator<Item = &[usize; 2]> {
                self.topology.iter_molecules()
            }
        
            unsafe fn nth_molecule_unchecked(&self, i: usize) -> &[usize;2] {
                self.topology.nth_molecule_unchecked(i)
            }
        }
        
        impl BondsProvider for $t {
            fn num_bonds(&self) -> usize {
                self.topology.num_bonds()
            }
        
            fn iter_bonds(&self) -> impl Iterator<Item = &[usize; 2]> {
                self.topology.iter_bonds()
            }
        
            unsafe fn nth_bond_unchecked(&self, i: usize) -> &[usize;2] {
                self.topology.nth_bond_unchecked(i)
            }
        }
    };
}

/// Macros for implementing common traits for read-write sources 
macro_rules! impl_read_write_source_traits {
    ( $t:ty ) => {
        impl_read_only_source_traits!($t);

        impl PosMutProvider for $t {
            fn iter_pos_mut(&self) -> impl PosMutIterator<'_> {
                self.state.iter_pos_mut()
            }
        }

        impl RandomPosMut for $t {
            fn nth_pos_mut(&self, i: usize) -> Option<&mut Pos> {
                self.state.nth_pos_mut(i)
            }

            unsafe fn nth_pos_mut_unchecked(&self, i: usize) -> &mut Pos {
                self.state.nth_pos_mut_unchecked(i)
            }
        }

        impl AtomsMutProvider for $t {
            fn iter_atoms_mut(&self) -> impl AtomMutIterator<'_> {
                self.topology.iter_atoms_mut()
            }
        }
        
        impl ParticleMutProvider for $t {
            fn iter_particle_mut(&self) -> impl ExactSizeIterator<Item = ParticleMut<'_>> {
                self.iter_index().map(|i| ParticleMut {
                    id: i,
                    atom: unsafe{self.topology.nth_atom_unchecked_mut(i)},
                    pos: unsafe{self.state.nth_pos_unchecked_mut(i)},
                })
            }
        }
    }
}

// All sources provide an index trivially
impl<K: SelectionKind> IndexProvider for Source<K> {
    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize> {
        0..self.topology.num_atoms()
    }
}

// For serial sources all traits are implemented
impl_read_write_source_traits!(Source<MutableSerial>);
impl_read_write_source_traits!(Source<BuilderSerial>);

// For immutable parallel source read-only traits are implemented
impl_read_only_source_traits!(Source<ImmutableParallel>);

// For mutable parallel source no traits are implemented!

// Bor Builder also BoxMut is implemnted
impl BoxMutProvider for Source<BuilderSerial> {
    fn get_box_mut(&self) -> Option<&mut PeriodicBox> {
        self.state.get_box_mut()
    }
}
