use std::{marker::PhantomData, ptr};
use triomphe::{Arc, UniqueArc};

use super::utils::*;
use crate::prelude::*;

/// Source of serial selections.
///
/// [Source] produces mutable selections, which may overlap and are can only be
/// used from the same thread where they are created.
///
/// [Source] takes ownership of [Topology] and [State] so that after creating a [Source] they are no longer accessible outside it.
/// This guarantees correct access without data integrity issues.

pub struct Source<K> {
    topology: triomphe::Arc<Topology>,
    state: triomphe::Arc<State>,
    _marker: PhantomData<K>,
}

impl Source<()> {
    /// Create the [Source] producing mutable selections that may overlap accessible from a single thread.
    pub fn new(
        topology: UniqueArc<Topology>,
        state: UniqueArc<State>,
    ) -> Result<Source<MutableSerial>, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        Ok(Source {
            topology: topology.shareable(),
            state: state.shareable(),
            _marker: Default::default(),
        })
    }

    pub fn from_system(system: System) -> Result<Source<MutableSerial>, SelectionError> {
        check_topology_state_sizes(&system.topology, &system.state)?;
        Ok(Source {
            topology: Arc::new(system.topology),
            state: Arc::new(system.state),
            _marker: Default::default(),
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
            _marker: Default::default(),
        })
    }

    pub fn empty_builder() -> Result<Source<BuilderSerial>, SelectionError> {
        Ok(Source {
            topology: Arc::new(Topology::default()),
            state: Arc::new(State::default()),
            _marker: Default::default(),
        })
    }

    pub fn from_file(fname: &str) -> Result<Source<MutableSerial>, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        Ok(Source::new(top, st)?)
    }

    pub fn from_file_builder(fname: &str) -> Result<Source<BuilderSerial>, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        Ok(Source::new_builder(top, st)?)
    }

    // Constructor for internal usage
    pub(crate) fn new_internal(
        topology: Arc<Topology>,
        state: Arc<State>,
    ) -> Result<Source<MutableSerial>, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        Ok(Source {
            topology,
            state,
            _marker: Default::default(),
        })
    }
}

impl<K: SerialSel> Source<K> {
    /// Release and return [Topology] and [State]. Fails if any selections created from this [Source] are still alive.
    pub fn release(self) -> Result<(Topology, State), SelectionError> {
        let top = Arc::try_unwrap(self.topology).or_else(|_| Err(SelectionError::Release))?;
        let st = Arc::try_unwrap(self.state).or_else(|_| Err(SelectionError::Release))?;
        Ok((top, st))
    }

    //---------------------------------
    // Creating selections
    //---------------------------------

    /// Creates new selection from an iterator of indexes. Indexes are bound checked, sorted and duplicates are removed.
    /// If any index is out of bounds the error is returned.
    pub fn select_from_iter(
        &self,
        iter: impl Iterator<Item = usize>,
    ) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_iter(iter, self.topology.num_atoms())?;
        Ok(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ))
    }

    pub unsafe fn select_from_vec_unchecked(
        &self,
        vec: Vec<usize>,    
    ) -> Sel<K> {
        Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            sorted_vec::SortedSet::from_sorted(vec),
        )
    }

    pub fn select_from_vec(
        &self,
        vec: Vec<usize>,    
    ) -> Sel<K> {
        Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            sorted_vec::SortedSet::from_unsorted(vec),
        )
    }

    /// Selects all
    pub fn select_all(&self) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_all(self.topology.num_atoms());
        Ok(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ))
    }

    /// Creates new selection from a selection expression string. Selection expression is constructed internally but
    /// can't be reused. Consider using [select_expr](Self::select_expr) if you already have selection expression.
    pub fn select_str(&self, selstr: &str) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_str(selstr, &self.topology, &self.state)?;
        Ok(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ))
    }

    /// Creates new selection from an existing selection expression.
    pub fn select_expr(&self, expr: &SelectionExpr) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_expr(expr, &self.topology, &self.state)?;
        Ok(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ))
    }

    /// Creates new selection from a range of indexes.
    /// If rangeis out of bounds the error is returned.
    pub fn select_range(&self, range: &std::ops::Range<usize>) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_range(range, self.topology.num_atoms())?;
        Ok(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ))
    }

    /// Sets new [State] in this source. All selections created from this source will automatically view
    /// the new state.
    ///
    /// New state should be compatible with the old one (have the same number of atoms). If not, the error is returned.
    ///
    /// Returns unique pointer to the old state, so it could be reused if needed.
    ///
    pub fn set_state(
        &mut self,
        state: UniqueArc<State>,
    ) -> Result<UniqueArc<State>, SelectionError> {
        // Check if the states are compatible
        if !self.state.interchangeable(&state) {
            return Err(SelectionError::SetState);
        }
        unsafe {
            ptr::swap(self.state.get_storage_mut(), state.get_storage_mut());
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
    pub fn set_topology(
        &mut self,
        topology: UniqueArc<Topology>,
    ) -> Result<UniqueArc<Topology>, SelectionError> {
        // Check if the states are compatible
        if !self.topology.interchangeable(&topology) {
            return Err(SelectionError::SetTopology);
        }
        unsafe {
            ptr::swap(self.topology.get_storage_mut(), topology.get_storage_mut());
        };
        Ok(topology)
    }
}

// Specific methods of serial source
impl Source<MutableSerial> {
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
            _marker: Default::default(),
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

//------------------
// IO traits impls
//-----------------
impl<K: SelectionKind> TopologyProvider for Source<K> {
    fn num_atoms(&self) -> usize {
        self.state.num_coords()
    }
}

impl<K: SelectionKind> AtomsProvider for Source<K> {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        self.topology.iter_atoms()
    }
}

impl<K: SelectionKind> StateProvider for Source<K> {
    fn get_time(&self) -> f32 {
        self.state.get_time()
    }

    fn num_coords(&self) -> usize {
        self.state.num_coords()
    }
}

impl<K: SelectionKind> PosProvider for Source<K> {
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
