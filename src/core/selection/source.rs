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
/// # Example 1: mutable non-overlapping selections
/// ```
/// # use molar::prelude::*;
/// # use anyhow::Result;
/// # use rayon::prelude::*;
/// # let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
/// let mut src = Source::new_parallel_mut(top.into(), st.into())?;
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
/// # Example 3: using subselections during parallel processing
/// ```
/// # use molar::prelude::*;
/// # use anyhow::Result;
/// # use rayon::prelude::*;
/// use itertools::Itertools;
/// # let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
/// let mut src = Source::new_parallel_mut(top.into(), st.into())?;
/// // Add a bunch of non-overlapping selections for residues
/// let mut sels = vec![];
/// for r in 545..550 {
///     sels.push( src.select_str(format!("resid {r}"))? );
/// }
/// // Process them in parallel. For each residue
/// // get CA and N atoms and compute a distance between them
/// let dist = sels.into_par_iter().map(|res| {
///     // Remove jumps over periodic boundary
///     res.unwrap_simple()?;
///     // Consume residue selection and convert it
///     // into two selections for Ca and N atoms
///     let (ca,n) = res.into_fragments(|p| {
///         match p.atom.name.as_str() {
///             "CA" => Some(0),
///             "N" => Some(1),
///             _ => None,
///         }
///     }).collect_tuple().unwrap();
///     Ok::<_,anyhow::Error>(ca.first_particle().pos-n.first_particle().pos)
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
/// let mut src = Source::new_parallel_mut(top.into(), st.into())?;
/// let sel = src.select_str("not resname TIP3 POT CLA")?;
/// thread::spawn( move || sel.translate(&Vector3f::new(10.0, 10.0, 10.0)));
/// #  Ok::<(), anyhow::Error>(())
/// ```

//------------------------------------------------------------------

pub struct Source<K: SelectionKind> {
    topology: Holder<Topology,K>,
    state: Holder<State,K>,
    _no_send: PhantomData<*const ()>,
    _kind: PhantomData<K>,
}

impl Source<MutableSerial> {
    /// Create the [Source] producing mutable selections that may overlap accessible from a single thread.
    pub fn new_serial(
        topology: Holder<Topology,MutableSerial>,
        state: Holder<State,MutableSerial>,
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
        topology: Holder<Topology,BuilderSerial>,
        state: Holder<State,BuilderSerial>,
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

    /// Creates a source of mutable parallel selections that _can't_ overlap.
    pub fn new_parallel_mut(
        topology: Holder<Topology,MutableParallel>,
        state: Holder<State,MutableParallel>,
    ) -> Result<Source<MutableParallel>, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        Ok(Source {
            topology: topology.into(),
            state: state.into(),            
            _no_send: Default::default(),
            _kind: Default::default(),
        })
    }

    /// Creates a source of immutable parallel selections that _may_ overlap.
    pub fn new_parallel(
        topology: Holder<Topology,ImmutableParallel>,
        state: Holder<State,ImmutableParallel>,
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

    pub fn parallel_mut_from_file(fname: &str) -> Result<Source<MutableParallel>, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        Ok(Source::new_parallel_mut(top.into(), st.into())?)
    }
}

impl<K: SelectionKind> Source<K> {
    /// Release and return [Topology] and [State]. Fails if any selections created from this [Source] are still alive.
    pub fn release(self) -> Result<(Topology, State), SelectionError> {
        Ok((
            self.topology.release()?,
            self.state.release()?,
        ))
    }

    fn new_sel(&self, index: SortedSet<usize>) -> Result<Sel<K>, SelectionError> {
        Sel::from_holders_and_index(
            self.topology.clone_with_index(&index)?,
            self.state.clone_with_index(&index)?,
            index,
        )
    }

    /// Creates new selection from an iterator of indexes. Indexes are bound checked, sorted and duplicates are removed.
    /// If any index is out of bounds the error is returned.
    pub fn select_iter(
        &mut self,
        iter: impl Iterator<Item = usize>,
    ) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_iter(iter, self.topology.num_atoms())?;        
        self.new_sel(vec)
    }

    pub fn select_vec(&mut self, vec: Vec<usize>) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_vec(vec, self.topology.num_atoms())?;
        self.new_sel(vec)
    }

    pub unsafe fn select_vec_unchecked(
        &mut self,
        vec: Vec<usize>,
    ) -> Result<Sel<K>, SelectionError> {
        let vec = SortedSet::from_sorted(vec);
        self.new_sel(vec)
    }

    /// Creates selection of all
    pub fn select_all(&mut self) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_all(self.topology.num_atoms());
        self.new_sel(vec)
    }

    /// Creates new selection from a selection expression string. Selection expression is constructed internally but
    /// can't be reused. Consider using [add_expr](Self::add_expr) if you already have selection expression.
    pub fn select_str(&mut self, selstr: impl AsRef<str>) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_str(selstr.as_ref(), &self.topology, &self.state)?;
        self.new_sel(vec)
    }

    /// Creates new selection from an existing selection expression.
    pub fn select_expr(&mut self, expr: &SelectionExpr) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_expr(expr, &self.topology, &self.state)?;
        self.new_sel(vec)
    }

    /// Creates new selection from a range of indexes.
    /// If rangeis out of bounds the error is returned.
    pub fn select_range(
        &mut self,
        range: &std::ops::Range<usize>,
    ) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_range(range, self.topology.num_atoms())?;
        self.new_sel(vec)
    }
}

// Specific methods for all sources except MutableParallel
impl<K> Source<K> 
where
    K: SelectionKind<UsedIndexesType = ()>,
{
    /// Sets new [State] in this source. All selections created from this source will automatically view
    /// the new state.
    ///
    /// New state should be compatible with the old one (have the same number of atoms). If not, the error is returned.
    ///
    /// Returns unique pointer to the old state, so it could be reused if needed.
    ///
    pub fn set_state(&mut self, state: Holder<State,K>) -> Result<Holder<State,K>, SelectionError> {
        // Check if the states are compatible
        if !self.state.interchangeable(&state) {
            return Err(SelectionError::SetState);
        }
        Ok(std::mem::replace(&mut self.state, state))
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
    pub fn set_topology(&mut self, topology: Holder<Topology,K>) -> Result<Holder<Topology,K>, SelectionError> {
        // Check if the states are compatible
        if !self.topology.interchangeable(&topology) {
            return Err(SelectionError::SetTopology);
        }

        Ok(std::mem::replace(&mut self.topology, topology))
    }

    pub fn get_topology(&self) -> Holder<Topology, K> {
        self.topology.clone()
    }

    pub fn get_state(&self) -> Holder<State, K> {
        self.state.clone()
    }
}

// Specific methods of builder source
impl Source<BuilderSerial> {
    pub fn append(&mut self, data: &(impl PosProvider + AtomsProvider)) -> Sel<BuilderSerial> {
        let first_added_index = self.num_atoms();
        self.topology
            .get_storage_mut()
            .add_atoms(data.iter_atoms().cloned());
        self.state
            .get_storage_mut()
            .add_coords(data.iter_pos().cloned());
        let last_added_index = self.num_atoms();
        self.select_range(&(first_added_index..last_added_index))
            .unwrap()
    }

    pub fn append_atoms(
        &mut self,
        atoms: impl Iterator<Item = Atom>,
        coords: impl Iterator<Item = Pos>,
    ) -> Sel<BuilderSerial> {
        let first_added_index = self.num_atoms();
        self.topology.get_storage_mut().add_atoms(atoms);
        self.state.get_storage_mut().add_coords(coords);
        let last_added_index = self.num_atoms();
        self.select_range(&(first_added_index..last_added_index))
            .unwrap()
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

    pub fn set_box(&mut self, new_box: Option<PeriodicBox>) {
        self.state.get_storage_mut().pbox = new_box;
    }
}

//--------------------------------------------------
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

//--------------------------------------------------
#[cfg(test)]
mod tests {
    use crate::prelude::*;
    use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

    #[test]
    fn par_iter2() -> anyhow::Result<()> {
        let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
        let mut src = Source::new_parallel_mut(top.into(), st.into())?;
        let mut sels = vec![];
        sels.push(src.select_str("resid 545")?);
        sels.push(src.select_str("resid 546")?);
        sels.push(src.select_str("resid 547")?);
        sels.par_iter_mut()
            .for_each(|sel| sel.translate(&Vector3f::new(10.0, 10.0, 10.0)));
        Ok::<(), anyhow::Error>(())
    }
}
