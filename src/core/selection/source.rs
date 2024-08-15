use std::ptr;
use crate::prelude::*;
use super::utils::*;

/// Source of serial selections.
/// 
/// [Source] produces mutable selections, which may overlap and are can only be
/// used from the same thread where they are created.
/// 
/// [Source] takes ownership of [Topology] and [State] so that after creating a [Source] they are no longer accessible outside it.
/// This guarantees correct access without data integrity issues.

pub struct Source<K> {
    topology: Topology<K>,
    state: State<K>,
}

impl Source<()> {
    /// Create the [Source] producing mutable selections that may overlap accessible from a single thread.
    pub fn new(
        topology: impl TryInto<Topology<MutableSerial>, Error=SelectionError>,
        state: impl TryInto<State<MutableSerial>, Error=SelectionError>,
    ) -> Result<Source<MutableSerial>, SelectionError> {
        let top = topology.try_into()?;
        let st = state.try_into()?;
        check_topology_state_sizes(&top, &st)?;
        Ok(Source {
            topology: top,
            state: st,
        })
    }

    pub fn new_builder(
        topology: impl TryInto<Topology<BuilderSerial>, Error=SelectionError>,
        state: impl TryInto<State<BuilderSerial>, Error=SelectionError>,
    ) -> Result<Source<BuilderSerial>, SelectionError> {
        let top = topology.try_into()?;
        let st = state.try_into()?;
        check_topology_state_sizes(&top, &st)?;
        Ok(Source {
            topology: top,
            state: st,
        })
    }

    pub fn from_file(fname: &str) -> Result<Source<MutableSerial>, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top,st) = fh.read()?;
        Ok(Source::new(top,st)?)
    }

    pub fn from_file_builder(fname: &str) -> Result<Source<BuilderSerial>, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top,st) = fh.read()?;
        Ok(Source::new_builder(top,st)?)
    }

    /*
    // Constructor for internal usage
    pub(crate) fn new_from_arc_system(system: triomphe::Arc<System>) -> Result<Source<MutableSerial>, SelectionError> {
        check_topology_state_sizes(&system.topology, &system.state)?;
        Ok(Source {
            system,
            _marker: Default::default(),
        })
    }
    */
}

impl<K: SerialSel> Source<K> {

    pub fn get_topology(self) -> Topology<K> {
        self.topology.new_arc()
    }

    pub fn get_state(self) -> State<K> {
        self.state.new_arc()
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
        Ok(Sel::new(
            self.topology.new_arc(),
            self.state.new_arc(),
            vec,            
        ))
    }

    /// Selects all
    pub fn select_all(&self) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_all(self.topology.num_atoms());
        Ok(Sel::new(
            self.topology.new_arc(),
            self.state.new_arc(),
            vec,
        ))
    }

    /// Creates new selection from a selection expression string. Selection expression is constructed internally but
    /// can't be reused. Consider using [select_expr](Self::select_expr) if you already have selection expression.
    pub fn select_str(&self, selstr: &str) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_str(selstr, &self.topology, &self.state)?;
        Ok(Sel::new(
            self.topology.new_arc(),
            self.state.new_arc(),
            vec,            
        ))
    }

    /// Creates new selection from an existing selection expression.
    pub fn select_expr(&self, expr: &SelectionExpr) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_expr(expr, &self.topology, &self.state)?;
        Ok(Sel::new(
            self.topology.new_arc(),
            self.state.new_arc(),
            vec,
        ))
    }

    /// Creates new selection from a range of indexes.
    /// If rangeis out of bounds the error is returned.
    pub fn select_range(&self, range: &std::ops::Range<usize>) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_range(range, self.topology.num_atoms())?;
        Ok(Sel::new(
            self.topology.new_arc(),
            self.state.new_arc(),
            vec,
        ))
    }

    /// Sets new [State] in this source. All selections created from this source will automatically view
    /// the new state. 
    /// 
    /// New state should be compatible with the old one (have the same number of atoms). If not, the error is returned.
    /// 
    /// Returns old state, so it could be reused if needed.
    /// 
    pub fn set_state(&mut self, state: State<K>) -> Result<State<K>, SelectionError> {
        // Check if the states are compatible
        if !self.state.interchangeable(&state) {
            return Err(SelectionError::SetState);
        }
        unsafe {
            ptr::swap(
                self.state.get_storage_mut(), 
                state.get_storage_mut()
            );
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
    pub fn set_topology(&mut self, topology: Topology<K>) -> Result<Topology<K>, SelectionError> {
        // Check if the states are compatible
        if !self.topology.interchangeable(&topology) {
            return Err(SelectionError::SetTopology);
        }

        unsafe {
            ptr::swap(
                self.topology.get_storage_mut(), 
                topology.get_storage_mut()
            );
        };

        Ok(topology)
    }
}

// Specific methods of builder source
impl Source<BuilderSerial> {    
    pub fn append(&mut self, data: &(impl PosProvider + AtomsProvider)) {
        self.topology.get_storage_mut().add_atoms(data.iter_atoms());
        self.state.get_storage_mut().add_coords(data.iter_pos());
    }

    pub fn remove(&mut self, to_remove: &impl IndexProvider) -> Result<(), SelectionError> {
        // We are checking index validity inside remove methods
        self.topology.get_storage_mut().remove_atoms(to_remove.iter_index())?;
        self.state.get_storage_mut().remove_coords(to_remove.iter_index())?;
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
