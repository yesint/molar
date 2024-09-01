use super::utils::*;
use crate::prelude::*;
use sorted_vec::SortedSet;
use std::marker::PhantomData;
use triomphe::{Arc, UniqueArc};

//----------------------------------------
// Source of parallel selections
//----------------------------------------

/// Source for selections, which are processed in parallel.
///
/// Selections are stored inside the source and the user can only access them directly by reference form
/// the same thread where [SourceParallel] was created (an attempt to send them to other thread won't compile).
/// In order to process selections in parallel user calls `map_par` method to run an arbitrary closure on each stored
/// selection in separate threads.
///
/// It is safe to change the [State] contained inside [SourceParallel] by calling [set_state](SourceParallel::set_state)
/// because it is guaranteed that no other threads are accessing stored selections at the same time.
///
/// # Example 1: mutable non-overlapping selections
/// ```
/// # use molar::prelude::*;
/// # use anyhow::Result;
/// # use rayon::prelude::*;
/// # let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
/// let mut src = Source::new_parallel_mut(top, st)?;
/// // Add a bunch of non-overlapping selections
/// let mut sels = vec![];
/// sels.push( src.select_iter(0..10)? );
/// sels.push( src.select_iter(11..15)? );
/// sels.push( src.select_iter(16..25)? );
/// // Process them
/// let v = Vector3f::new(1.0, 2.0, 3.0);
/// let res = sels.par_iter().map(|sel| {
///    sel.translate(&v);
///    Ok(sel.center_of_mass()?)
/// }).collect::<Result<Vec<_>>>()?;
/// println!("{:?}",res);
/// #  Ok::<(), anyhow::Error>(())
/// ```
/// # Example 2: immutable overlapping selections, using par_iter()
/// ```
/// # use molar::prelude::*;
/// # use anyhow::Result;
/// # use rayon::prelude::*;
/// # let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
/// let mut src = Source::new_parallel(top, st)?;
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
/// # Example 3: using subselections during parallel processing
/// ```
/// # use molar::prelude::*;
/// # use anyhow::Result;
/// # use rayon::prelude::*;
/// # let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
/// let mut src = Source::new_parallel(top, st)?;
/// // Add a bunch of non-overlapping selections for residues
/// let mut sels = vec![];
/// for r in 545..550 {
///     sels.push( src.select_str(format!("resid {r}"))? );
/// }
/// // Process them in parallel. For each residue
/// // get CA and N atoms and compute a distance between them
/// let dist = sels.par_iter().map(|res| {
///    //res.unwrap_simple()?;
///    let ca = res.subsel_from_str("name CA")?;
///    let n = res.subsel_from_str("name N")?;
///    Ok::<_,anyhow::Error>(ca.first_particle().pos-n.first_particle().pos)
/// }).collect::<Result<Vec<_>>>()?;
/// println!("{:?}",dist);
/// #  Ok::<(), anyhow::Error>(())
/// ```

/// # It's Ok to move parallel selection to the other thread
///```
/// # use molar::prelude::*;
/// # use std::thread;
/// # let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
/// // This won't compile
/// let mut src = Source::new_parallel_mut(top, st)?;
/// let sel = src.select_str("not resname TIP3 POT CLA")?;
/// thread::spawn( move || sel.translate(&Vector3f::new(10.0, 10.0, 10.0)));
/// #  Ok::<(), anyhow::Error>(())
/// ```

pub struct Source<K> {
    topology: Arc<Topology>,
    state: Arc<State>,
    used: rustc_hash::FxHashSet<usize>,
    _no_send: PhantomData<*const ()>,
    _kind: PhantomData<K>,
}

impl Source<()> {
    /// Create the [Source] producing mutable selections that may overlap accessible from a single thread.
    pub fn new_serial(
        topology: UniqueArc<Topology>,
        state: UniqueArc<State>,
    ) -> Result<Source<MutableSerial>, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        Ok(Source {
            topology: topology.shareable(),
            state: state.shareable(),
            used: Default::default(),
            _no_send: Default::default(),
            _kind: Default::default()
        })
    }

    pub fn new_builder(
        topology: UniqueArc<Topology>,
        state: UniqueArc<State>,
    ) -> Result<Source<BuilderSerial>, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        Ok(Source {
            topology: topology.shareable(),
            state: state.shareable(),
            used: Default::default(),
            _no_send: Default::default(),
            _kind: Default::default()
        })
    }

    pub fn empty_builder() -> Source<BuilderSerial> {
        Source {
            topology: Arc::new(Topology::default()),
            state: Arc::new(State::default()),
            used: Default::default(),
            _no_send: Default::default(),
            _kind: Default::default()
        }
    }

    /// Creates a source of mutable parallel selections that _can't_ overlap.
    pub fn new_parallel_mut(
        topology: UniqueArc<Topology>,
        state: UniqueArc<State>,
    ) -> Result<Source<MutableParallel>, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        Ok(Source {
            topology: topology.shareable(),
            state: state.shareable(),
            used: Default::default(),
            _no_send: Default::default(),
            _kind: Default::default()
        })
    }

    /// Creates a source of immutable parallel selections that _may_ overlap.
    pub fn new_parallel(
        topology: UniqueArc<Topology>,
        state: UniqueArc<State>,
    ) -> Result<Source<ImmutableParallel>, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        Ok(Source {
            topology: topology.shareable(),
            state: state.shareable(),
            used: Default::default(),
            _no_send: Default::default(),
            _kind: Default::default()
        })
    }

    pub fn parallel_from_file(fname: &str) -> Result<Source<ImmutableParallel>, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        Ok(Source::new_parallel(top, st)?)
    }

    pub fn serial_from_file(fname: &str) -> Result<Source<MutableSerial>, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        Ok(Source::new_serial(top, st)?)
    }

    pub fn builder_from_file(fname: &str) -> Result<Source<BuilderSerial>, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        Ok(Source::new_builder(top, st)?)
    }

    pub fn parallel_mut_from_file(fname: &str) -> Result<Source<MutableParallel>, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        Ok(Source::new_parallel_mut(top, st)?)
    }
}

impl<K: SelectionKind> Source<K> {
    /// Release and return [Topology] and [State]. Fails if any selections created from this [Source] are still alive.
    pub fn release(self) -> Result<(UniqueArc<Topology>, UniqueArc<State>), SelectionError> {
        Ok((
            Arc::try_unique(self.topology).map_err(|_| SelectionError::Release)?, 
            Arc::try_unique(self.state).map_err(|_| SelectionError::Release)?,
        ))
    }

    /// Adds new selection from an iterator of indexes. Indexes are bound checked, sorted and duplicates are removed.
    /// If any index is out of bounds the error is returned.
    pub fn select_iter(
        &mut self,
        iter: impl Iterator<Item = usize>,
    ) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_iter(iter, self.topology.num_atoms())?;
        K::check_overlap(&vec, &mut self.used)?;
        Ok(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ))
    }

    pub fn select_vec(
        &mut self,
        vec: &Vec<usize>,
    ) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_vec(vec, self.topology.num_atoms())?;
        K::check_overlap(&vec, &mut self.used)?;
        Ok(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ))
    }

    pub unsafe fn select_vec_unchecked(
        &mut self,
        vec: Vec<usize>,
    ) -> Result<Sel<K>, SelectionError> {
        let vec = SortedSet::from_sorted(vec);
        K::check_overlap(&vec, &mut self.used)?;
        Ok(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ))
    }

    /// Adds selection of all
    pub fn select_all(&mut self) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_all(self.topology.num_atoms());
        K::check_overlap(&vec, &mut self.used)?;
        Ok(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ))
    }

    /// Adds new selection from a selection expression string. Selection expression is constructed internally but
    /// can't be reused. Consider using [add_expr](Self::add_expr) if you already have selection expression.
    pub fn select_str(&mut self, selstr: impl AsRef<str>) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_str(selstr.as_ref(), &self.topology, &self.state)?;
        K::check_overlap(&vec, &mut self.used)?;
        Ok(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ))
    }

    /// Adds new selection from an existing selection expression.
    pub fn select_expr(&mut self, expr: &SelectionExpr) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_expr(expr, &self.topology, &self.state)?;
        K::check_overlap(&vec, &mut self.used)?;
        Ok(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ))
    }

    /// Adds new selection from a range of indexes.
    /// If rangeis out of bounds the error is returned.
    pub fn select_range(&mut self, range: &std::ops::Range<usize>) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_range(range, self.topology.num_atoms())?;
        K::check_overlap(&vec, &mut self.used)?;
        Ok(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ))
    }
}


// Specific methods for serial selections 
impl<K: SerialSel> Source<K> {
    /// Sets new [State] in this source. All selections created from this source will automatically view
    /// the new state.
    ///
    /// New state should be compatible with the old one (have the same number of atoms). If not, the error is returned.
    ///
    /// Returns unique pointer to the old state, so it could be reused if needed.
    ///
    pub fn set_state(&mut self, state: State) -> Result<State, SelectionError> {
        // Check if the states are compatible
        if !self.state.interchangeable(&state) {
            return Err(SelectionError::SetState);
        }
        unsafe {
            std::ptr::swap(self.state.get_storage_mut(), state.get_storage_mut());
        };

        Ok(state)
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
    pub fn set_topology(&mut self, topology: Topology) -> Result<Topology, SelectionError> {
        // Check if the states are compatible
        if !self.topology.interchangeable(&topology) {
            return Err(SelectionError::SetTopology);
        }

        unsafe {
            std::ptr::swap(
                self.topology.get_storage_mut(),
                topology.get_storage_mut(),
            );
        };

        Ok(topology)
    }

    pub fn get_shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    pub fn get_shared_state(&self) -> Arc<State> {
        Arc::clone(&self.state)
    }

    pub fn new_from_shared(
        topology: &Arc<Topology>,
        state: &Arc<State>,
    ) -> Result<Source<MutableSerial>, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        Ok(Source {
            topology: Arc::clone(topology),
            state: Arc::clone(state),
            used: Default::default(),
            _no_send: Default::default(),
            _kind: Default::default(),
        })
    }

    pub fn set_shared_topology(
        &mut self,
        topology: Arc<Topology>,
    ) -> Result<Arc<Topology>, SelectionError> {
        if !self.topology.interchangeable(&topology) {
            return Err(SelectionError::SetTopology);
        }
        Ok(std::mem::replace(&mut self.topology, topology))
    }

    pub fn set_shared_state(&mut self, state: Arc<State>) -> Result<Arc<State>, SelectionError> {
        if !self.state.interchangeable(&state) {
            return Err(SelectionError::SetState);
        }
        Ok(std::mem::replace(&mut self.state, state))
    }
}

// Specific methods of builder source
impl Source<BuilderSerial> {
    pub fn append(&mut self, data: &(impl PosProvider + AtomsProvider)) -> Sel<BuilderSerial> {
        let first_added_index = self.num_atoms();
        self.topology.get_storage_mut().add_atoms(data.iter_atoms());
        self.state.get_storage_mut().add_coords(data.iter_pos());
        let last_added_index = self.num_atoms();
        self.select_range(&(first_added_index..last_added_index)).unwrap()
    }

    pub fn remove(&mut self, to_remove: &impl IndexProvider) -> Result<(), SelectionError> {
        // We are checking index validity inside remove methods
        self.topology
            .get_storage_mut()
            .remove_atoms(to_remove.iter_index())?;
        self.state
            .get_storage_mut()
            .remove_coords(to_remove.iter_index())?;
        Ok(())
    }
}

//--------------------------------------------------
// For serial selections all traits are implemented
//--------------------------------------------------
impl TopologyProvider for Source<MutableSerial> {
    fn num_atoms(&self) -> usize {
        self.state.num_coords()
    }
}

impl AtomsProvider for Source<MutableSerial> {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        self.topology.iter_atoms()
    }
}

impl StateProvider for Source<MutableSerial> {
    fn get_time(&self) -> f32 {
        self.state.get_time()
    }

    fn num_coords(&self) -> usize {
        self.state.num_coords()
    }
}

impl PosProvider for Source<MutableSerial> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.state.iter_pos()
    }
}

impl BoxProvider for Source<MutableSerial> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.state.get_box()
    }
}

impl WritableToFile for Source<MutableSerial> {}

//--------------------------------------------------
// For builder selections all traits are implemented
//--------------------------------------------------
impl TopologyProvider for Source<BuilderSerial> {
    fn num_atoms(&self) -> usize {
        self.state.num_coords()
    }
}

impl AtomsProvider for Source<BuilderSerial> {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        self.topology.iter_atoms()
    }
}

impl StateProvider for Source<BuilderSerial> {
    fn get_time(&self) -> f32 {
        self.state.get_time()
    }

    fn num_coords(&self) -> usize {
        self.state.num_coords()
    }
}

impl PosProvider for Source<BuilderSerial> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.state.iter_pos()
    }
}

impl BoxProvider for Source<BuilderSerial> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.state.get_box()
    }
}

impl WritableToFile for Source<BuilderSerial> {}


#[cfg(test)]
mod tests {
    use crate::prelude::*;
    use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};
    #[test]
    fn par_iter2() -> anyhow::Result<()> {
        let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
        let mut src = Source::new_parallel_mut(top, st)?;
        let mut sels = vec![];
        sels.push( src.select_str("resid 545")? );
        sels.push( src.select_str("resid 546")? );
        sels.push( src.select_str("resid 547")? );
        sels.par_iter_mut()
            .for_each(|sel| sel.translate(&Vector3f::new(10.0, 10.0, 10.0)));
        Ok::<(), anyhow::Error>(())
    }
}
