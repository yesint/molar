use std::marker::PhantomData;

use crate::prelude::*;
use anyhow::{bail, Context, Result, anyhow};
use rayon::prelude::*;
use super::utils::*;

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
/// let res: Vec<_> = src.map_par(|sel| {
///    sel.translate(&v);
///    Ok(sel.center_of_mass()?)
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
/// let dist: Vec<_> = src.map_par(|res| {
///    res.unwrap_simple()?;
///    let ca = res.subsel_from_str("name CA")?;
///    let n = res.subsel_from_str("name N")?;
///    Ok(ca.first_particle().pos-n.first_particle().pos)
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
    system: triomphe::Arc<System>,
    selections: Vec<Sel<K>>,
    used: rustc_hash::FxHashSet<usize>,
    _marker: PhantomData<*const ()>,
}

impl SourceParallel<()> {    
    /// Creates a source of mutable parallel selections that _can't_ overlap.
    pub fn new_mut (
        topology: Topology,
        state: State,
    ) -> Result<SourceParallel<MutableParallel>> {
        check_topology_state_sizes(&topology, &state)?;
        Ok(SourceParallel {
            system: triomphe::Arc::new(System{topology,state}),
            selections: Default::default(),
            used: Default::default(),
            _marker: Default::default(),
        })
    }

    /// Creates a source of immutable parallel selections that _may_ overlap.
    pub fn new (
        topology: Topology,
        state: State,
    ) -> Result<SourceParallel<ImmutableParallel>> {
        check_topology_state_sizes(&topology, &state)?;
        Ok(SourceParallel {
            system: triomphe::Arc::new(System{topology,state}),
            selections: Default::default(),
            used: Default::default(),
            _marker: Default::default(),
        })
    }
}

impl<K: ParallelSel> SourceParallel<K> {
    /// Release and return [Topology] and [State]. Fails if any selections created from this [Source] are still alive.
    pub fn release(self) -> anyhow::Result<(Topology, State)> {
        let sys = triomphe::Arc::try_unwrap(self.system)
                .or_else(|_| bail!("Can't release Source: multiple references are active!"))?;
        Ok((sys.topology,sys.state))
    }

    /// Executes provided closure on each stored selection in parallel and
    /// collect the closures' outputs into a container.
    /// 
    /// Closure return `Result<RT>`, where `RT` is anything 
    /// sharable between the threads. 
    pub fn map_par<RT,C,F>(&self, func: F) -> Result<C> 
    where
        RT: Send + Sync,
        F: Fn(&Sel<K>)->Result<RT> + Send + Sync,
        C: rayon::iter::FromParallelIterator<RT>,
    {
        self.selections.par_iter().map(func).collect()
    }

    /// Executes provide closure on each stored selection serially
    /// and collect the closures' outputs into a container.
    pub fn map<RT,C,F>(&self, func: F) -> Result<C> 
    where
        RT: Send + Sync,
        F: Fn(&Sel<K>)->Result<RT> + Send + Sync,
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
    pub fn add_existing(&mut self, sel: Sel<MutableSerial>) -> anyhow::Result<usize> {
        self.add_from_iter(sel.iter_index())
    }

    /// Adds new selection from an iterator of indexes. Indexes are bound checked, sorted and duplicates are removed.
    /// If any index is out of bounds the error is returned.
    pub fn add_from_iter(
        &mut self,
        iter: impl Iterator<Item = usize>,
    ) -> anyhow::Result<usize> {
        let vec = index_from_iter(iter, self.system.topology.num_atoms())?;        
        K::check_overlap(&vec, &mut self.used)?;
        self.selections.push(Sel::new(
            triomphe::Arc::clone(&self.system),
            vec,            
        ));
        Ok(self.selections.len()-1)
    }

    /// Adds selection of all
    pub fn select_all(&mut self) -> anyhow::Result<&Sel<K>> {
        let vec = index_from_all(self.system.topology.num_atoms());        
        K::check_overlap(&vec, &mut self.used)?;
        self.selections.push(Sel::new(
            triomphe::Arc::clone(&self.system),
            vec,            
        ));
        Ok(&self.selections.last().unwrap())
    }

    /// Adds new selection from a selection expression string. Selection expression is constructed internally but
    /// can't be reused. Consider using [add_expr](Self::add_expr) if you already have selection expression.
    pub fn add_str(&mut self, selstr: impl AsRef<str>) -> anyhow::Result<&Sel<K>> {
        let vec = index_from_str(selstr.as_ref(), &self.system.topology, &self.system.state)?;        
        K::check_overlap(&vec, &mut self.used)
            .with_context(|| format!("Adding str selection '{}'",selstr.as_ref()))?;
        self.selections.push(Sel::new(
            triomphe::Arc::clone(&self.system),
            vec,            
        ));
        Ok(&self.selections.last().unwrap())
    }

    pub fn get_sel(&self, i: usize) -> anyhow::Result<&Sel<K>> {
        self.selections.get(i).ok_or_else(|| anyhow!("Invalid selection index"))
    }

    /// Adds new selection from an existing selection expression.
    pub fn add_expr(&mut self, expr: &SelectionExpr) -> anyhow::Result<&Sel<K>> {
        let vec = index_from_expr(expr, &self.system.topology, &self.system.state)?;
        K::check_overlap(&vec, &mut self.used)?;
        self.selections.push(Sel::new(
            triomphe::Arc::clone(&self.system),
            vec,            
        ));
        Ok(&self.selections.last().unwrap())
    }

    /// Adds new selection from a range of indexes.
    /// If rangeis out of bounds the error is returned.
    pub fn add_range(&mut self, range: &std::ops::Range<usize>) -> anyhow::Result<&Sel<K>> {
        let vec = index_from_range(range, self.system.topology.num_atoms())?;
        K::check_overlap(&vec, &mut self.used)?;
        self.selections.push(Sel::new(
            triomphe::Arc::clone(&self.system),
            vec,            
        ));
        Ok(&self.selections.last().unwrap())
    }

    /// Converts `Self` into a serial [Source]. All stored parallel selections 
    /// are converted into serial mutable selections and returned as a vector. 
    pub fn into_serial_with_sels(self) -> (Source<MutableSerial>,Vec<Sel<MutableSerial>>) {                
        let sels: Vec<Sel<MutableSerial>> = self.selections.into_iter().map(|sel|             
            Sel::from_parallel(sel)
        ).collect();
        // This should never fail
        let src = Source::new_from_system(self.system).unwrap(); 
        (src,sels)
    }

    /// Converts `Self` into a serial [Source]. All stored parallel selections 
    /// are dropped. 
    pub fn into_serial(mut self) -> Source<MutableSerial> {
        // Drop all selections
        self.selections.clear();
        let (top,st) = self.release().unwrap();
        // This should never fail
        Source::new(top,st).unwrap()
    }

    /// Sets new [State] in this source. All selections created from this source will automatically view
    /// the new state. 
    /// 
    /// New state should be compatible with the old one (have the same number of atoms). If not, the error is returned.
    /// 
    /// Returns unique pointer to the old state, so it could be reused if needed.
    /// 
    pub fn set_state(&mut self, state: State) -> Result<State> {
        // Check if the states are compatible
        if !self.system.state.interchangeable(&state) {
            bail!("Can't set state: states are incompatible!")
        }
        unsafe {
            std::ptr::swap(
                self.system.state.get_storage_mut(), 
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
    pub fn set_topology(&mut self, topology: Topology) -> Result<Topology> {
        // Check if the states are compatible
        if !self.system.topology.interchangeable(&topology) {
            bail!("Can't set topology: topologies are incompatible!")
        }

        unsafe {
            std::ptr::swap(
                self.system.topology.get_storage_mut(), 
                topology.get_storage_mut()
            );
        };

        Ok(topology)
    }
    
}

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
        src.par_iter().for_each(|sel| sel.translate(&Vector3f::new(10.0, 10.0, 10.0)));
        Ok::<(), anyhow::Error>(())
    }
}