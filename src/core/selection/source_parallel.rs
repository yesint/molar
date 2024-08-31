use super::utils::*;
use crate::prelude::*;
use rayon::prelude::*;
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
/// # let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
/// let mut src = SourceParallel::new_mut(top, st)?;
/// // Add a bunch of non-overlapping selections
/// src.add_from_iter(0..10)?;
/// src.add_from_iter(11..15)?;
/// src.add_from_iter(16..25)?;
/// // Process them
/// let v = Vector3f::new(1.0, 2.0, 3.0);
/// let res: Vec<_> = src.collect_par(|sel| {
///    sel.translate(&v);
///    Ok::<_,anyhow::Error>(sel.center_of_mass()?)
/// })?;
/// println!("{:?}",res);
/// #  Ok::<(), anyhow::Error>(())
/// ```
/// # Example 2: immutable overlapping selections, using par_iter()
/// ```
/// # use molar::prelude::*;
/// # use anyhow::Result;
/// # let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
/// let mut src = SourceParallel::new(top, st)?;
/// // Overlapping selections
/// src.add_from_iter(0..10)?;
/// src.add_from_iter(5..15)?;
/// // Process them using par_iter explicitly
/// let res = src.par_iter().map(|sel| {
///    Ok(sel.center_of_mass()?)
/// }).collect::<Result<Vec<_>>>()?;
/// println!("{:?}",res);
/// #  Ok::<(), anyhow::Error>(())
/// ```
/// # Example 3: using subselections during parallel processing
/// ```
/// # use molar::prelude::*;
/// # use anyhow::Result;
/// # let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
/// let mut src = SourceParallel::new_mut(top, st)?;
/// // Add a bunch of non-overlapping selections for residues
/// for r in 545..550 {
///     src.add_str(format!("resid {r}"))?;
/// }
/// // Process them in parallel. For each residue
/// // get CA and N atoms and compute a distance between them
/// let dist: Vec<_> = src.collect_par(|res| {
///    res.unwrap_simple()?;
///    let ca = res.subsel_from_str("name CA")?;
///    let n = res.subsel_from_str("name N")?;
///    Ok::<_,anyhow::Error>(ca.first_particle().pos-n.first_particle().pos)
/// })?;
/// println!("{:?}",dist);
/// #  Ok::<(), anyhow::Error>(())
/// ```

/// # Safety guarantees
/// An attampt to move a selection to other thread fails to compile:
///```compile_fail
/// # use molar::prelude::*;
/// # use std::thread;
/// # let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
/// // This won't compile
/// let mut src = SourceParallel::new_mut(top, st)?;
/// let sel = src.add_str("not resname TIP3 POT CLA")?;
/// thread::spawn( move || sel.translate(&Vector3f::new(10.0, 10.0, 10.0)));
/// #  Ok::<(), anyhow::Error>(())
/// ```
///
/// It is also impossible to leak the [Sel] from an iterator:
///```compile_fail
/// # use molar::prelude::*;
/// # use std::thread;
/// # let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
/// // This won't compile
/// let mut src = SourceParallel::new_mut(top, st)?;
/// src.add_str("resid 545")?;
/// src.add_str("resid 546")?;
/// src.add_str("resid 547")?;
/// let sel = src.iter().next().unwrap();
/// thread::spawn( move || sel.translate(&Vector3f::new(10.0, 10.0, 10.0)));
/// #  Ok::<(), anyhow::Error>(())
/// ```
/// However, it's possible to use parallel iterator in usual way:
///```
/// # use molar::prelude::*;
/// # let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
/// let mut src = SourceParallel::new_mut(top, st)?;
/// src.add_str("resid 545")?;
/// src.add_str("resid 546")?;
/// src.add_str("resid 547")?;
/// src.par_iter().for_each(|sel| sel.translate(&Vector3f::new(10.0, 10.0, 10.0)));
/// #  Ok::<(), anyhow::Error>(())
/// ```

pub struct SourceParallel<K> {
    topology: Arc<Topology>,
    state: Arc<State>,
    selections: Vec<Sel<K>>,
    used: rustc_hash::FxHashSet<usize>,
    _marker: PhantomData<*const ()>,
}

impl SourceParallel<()> {
    /// Creates a source of mutable parallel selections that _can't_ overlap.
    pub fn new_mut(
        topology: UniqueArc<Topology>,
        state: UniqueArc<State>,
    ) -> Result<SourceParallel<MutableParallel>, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        Ok(SourceParallel {
            topology: topology.shareable(),
            state: state.shareable(),
            selections: Default::default(),
            used: Default::default(),
            _marker: Default::default(),
        })
    }

    /// Creates a source of immutable parallel selections that _may_ overlap.
    pub fn new(
        topology: UniqueArc<Topology>,
        state: UniqueArc<State>,
    ) -> Result<SourceParallel<ImmutableParallel>, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        Ok(SourceParallel {
            topology: topology.shareable(),
            state: state.shareable(),
            selections: Default::default(),
            used: Default::default(),
            _marker: Default::default(),
        })
    }

    pub fn from_file(fname: &str) -> Result<SourceParallel<ImmutableParallel>, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        Ok(SourceParallel::new(top, st)?)
    }

    pub fn from_file_mut(fname: &str) -> Result<SourceParallel<MutableParallel>, SelectionError> {
        let mut fh = FileHandler::open(fname)?;
        let (top, st) = fh.read()?;
        Ok(SourceParallel::new_mut(top, st)?)
    }
}

impl<K: ParallelSel> SourceParallel<K> {
    /// Release and return [Topology] and [State]. Fails if any selections created from this [Source] are still alive.
    pub fn release(self) -> Result<(UniqueArc<Topology>, UniqueArc<State>), SelectionError> {
        Ok((
            Arc::try_unique(self.topology).map_err(|_| SelectionError::Release)?, 
            Arc::try_unique(self.state).map_err(|_| SelectionError::Release)?,
        ))
    }

    /// Executes provided closure on each stored selection in parallel and
    /// collect the closures' outputs into a container.
    ///
    /// Closure return `Result<RT,_>`, where `RT` is anything
    /// sharable between the threads.
    pub fn collect_par<RT, C, F, E>(&self, func: F) -> Result<C, E>
    where
        RT: Send + Sync,
        E: Send + Sync,
        F: Fn(&Sel<K>) -> Result<RT, E> + Send + Sync,
        C: rayon::iter::FromParallelIterator<RT>,
    {
        self.selections.par_iter().map(func).collect()
    }

    /// Executes provide closure on each stored selection serially
    /// and collect the closures' outputs into a container.
    pub fn collect<RT, C, F, E>(&self, func: F) -> Result<C, E>
    where
        F: Fn(&Sel<K>) -> Result<RT, E> + Send + Sync,
        C: FromIterator<RT>,
    {
        self.selections.iter().map(func).collect()
    }

    /// Returns serial iterator over stored selections
    pub fn iter(&self) -> impl Iterator<Item = &Sel<K>> {
        self.selections.iter()
    }

    /// Returns parallel iterator over stored selections
    pub fn par_iter(&self) -> impl ParallelIterator<Item = &Sel<K>> {
        self.selections.par_iter()
    }

    /// Adds an existing serial selection. Passed selection is consumed
    /// and its index is re-evaluated against the new [SourceParallel].
    pub fn add_existing(&mut self, sel: Sel<MutableSerial>) -> Result<&Sel<K>, SelectionError> {
        self.add_from_iter(sel.iter_index())
    }

    /// Adds new selection from an iterator of indexes. Indexes are bound checked, sorted and duplicates are removed.
    /// If any index is out of bounds the error is returned.
    pub fn add_from_iter(
        &mut self,
        iter: impl Iterator<Item = usize>,
    ) -> Result<&Sel<K>, SelectionError> {
        let vec = index_from_iter(iter, self.topology.num_atoms())?;
        K::check_overlap(&vec, &mut self.used)?;
        self.selections.push(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ));
        Ok(&self.selections.last().unwrap())
    }

    /// Adds selection of all
    pub fn select_all(&mut self) -> Result<&Sel<K>, SelectionError> {
        let vec = index_from_all(self.topology.num_atoms());
        K::check_overlap(&vec, &mut self.used)?;
        self.selections.push(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ));
        Ok(&self.selections.last().unwrap())
    }

    /// Adds new selection from a selection expression string. Selection expression is constructed internally but
    /// can't be reused. Consider using [add_expr](Self::add_expr) if you already have selection expression.
    pub fn add_str(&mut self, selstr: impl AsRef<str>) -> Result<&Sel<K>, SelectionError> {
        let vec = index_from_str(selstr.as_ref(), &self.topology, &self.state)?;
        K::check_overlap(&vec, &mut self.used)?;
        self.selections.push(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ));
        Ok(&self.selections.last().unwrap())
    }

    pub fn get_sel(&self, i: usize) -> Result<&Sel<K>, SelectionError> {
        self.selections
            .get(i)
            .ok_or_else(|| SelectionError::OutOfBounds(i, self.state.num_coords()))
    }

    /// Adds new selection from an existing selection expression.
    pub fn add_expr(&mut self, expr: &SelectionExpr) -> Result<&Sel<K>, SelectionError> {
        let vec = index_from_expr(expr, &self.topology, &self.state)?;
        K::check_overlap(&vec, &mut self.used)?;
        self.selections.push(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ));
        Ok(&self.selections.last().unwrap())
    }

    /// Adds new selection from a range of indexes.
    /// If rangeis out of bounds the error is returned.
    pub fn add_range(&mut self, range: &std::ops::Range<usize>) -> Result<&Sel<K>, SelectionError> {
        let vec = index_from_range(range, self.topology.num_atoms())?;
        K::check_overlap(&vec, &mut self.used)?;
        self.selections.push(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ));
        Ok(&self.selections.last().unwrap())
    }

    /// Converts `Self` into a serial [Source]. All stored parallel selections
    /// are converted into serial mutable selections and returned as a vector.
    pub fn into_source_with_sels(self) -> (Source<MutableSerial>, Vec<Sel<MutableSerial>>) {
        let sels: Vec<Sel<MutableSerial>> = self
            .selections
            .into_iter()
            .map(|sel| Sel::from_parallel(sel))
            .collect();
        // This should never fail
        let src = Source::new_internal(self.topology, self.state).unwrap();
        (src, sels)
    }

    /// Converts `Self` into a serial [Source]. All stored parallel selections
    /// are dropped.
    pub fn into_source(mut self) -> Source<MutableSerial> {
        // Drop all selections
        self.selections.clear();
        let (top, st) = self.release().unwrap();
        // This should never fail
        Source::new(top, st).unwrap()
    }

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
}

//-----------------------------------------
// IO traits impls
// Only for immutable parallel selections!
//-----------------------------------------
impl TopologyProvider for SourceParallel<ImmutableParallel> {
    fn num_atoms(&self) -> usize {
        self.state.num_coords()
    }
}

impl AtomsProvider for SourceParallel<ImmutableParallel> {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        self.topology.iter_atoms()
    }
}

impl StateProvider for SourceParallel<ImmutableParallel> {
    fn get_time(&self) -> f32 {
        self.state.get_time()
    }

    fn num_coords(&self) -> usize {
        self.state.num_coords()
    }
}

impl PosProvider for SourceParallel<ImmutableParallel> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.state.iter_pos()
    }
}

impl BoxProvider for SourceParallel<ImmutableParallel> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.state.get_box()
    }
}

impl WritableToFile for SourceParallel<ImmutableParallel> {}

#[cfg(test)]
mod tests {
    use crate::prelude::*;
    use rayon::iter::ParallelIterator;
    #[test]
    fn par_iter2() -> anyhow::Result<()> {
        let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
        let mut src = SourceParallel::new_mut(top, st)?;
        src.add_str("resid 545")?;
        src.add_str("resid 546")?;
        src.add_str("resid 547")?;
        src.par_iter()
            .for_each(|sel| sel.translate(&Vector3f::new(10.0, 10.0, 10.0)));
        Ok::<(), anyhow::Error>(())
    }
}
