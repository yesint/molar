use std::{marker::PhantomData, ptr};
use crate::prelude::*;
use sorted_vec::SortedSet;
use std::{
    marker::PhantomData,
    ops::Deref,
    sync::Mutex,
};

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
/// }).collect::<Result<Vec<_>>>()?;topology
/// println!("{:?}",res);
/// #  Ok::<(), anyhow::Error>(())
/// ```
/// # Example 2: immutable overlapping selections, using par_iter()
/// ```
/// # use molar::prelude::*;
/// # use anyhow::Result;
/// # use rayon::prelude::*;
/// # let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
/// let mut src = Source::new_parallel(top, st)?;topology
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
/// use itertools::Itertools;topology
/// # let (top, st) = FileHandler::open("tests/protein.pdb")?.read()?;
/// let mut src = Source::new_parallel_mut(top, st)?;
/// // Add a bunch of non-overlapping selections for residues
/// let mut sels = vec![];
/// for r in 545..550 {
///     sels.push( src.select_str(format!("resid {r}"))? );
/// }
/// // Process them in parallel. For each residue
/// // get CA and N atoms and compute a distance between them
/// let dist = sels.into_par_iter().map(|res| {
///     // Remove jusmp over periodic boundary
///     res.unwrap_simple()?;
///     // Consume residue selection and convert it
///     // into two selections for Ca and N atoms
///     let (ca,n) = res.into_split_contig(|p| {
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
/// let mut src = Source::new_parallel_mut(top, st)?;
/// let sel = src.select_str("not resname TIP3 POT CLA")?;
/// thread::spawn( move || sel.translate(&Vector3f::new(10.0, 10.0, 10.0)));
/// #  Ok::<(), anyhow::Error>(())
/// ```

//-------------------------------------------------------------------------

/// Smart pointer wrapper for sharing [Topology] and [State] between serial selections.
/// Acts like Rc parameterized by selection kind, so the user can't accidentally mix incompatible
/// selections pointing to the same data.
/// Can't be sent to other threads.
/// Normally this type should not be used directly by the user.
pub(crate) struct BaseHolder<T, K, U> {
    pub(crate) arc: triomphe::Arc<T>,
    used: U,
    _kind: PhantomData<K>,
}

pub type UniqueHolder<T> = BaseHolder<T,(),()>;
pub type Holder<T,K: SelectionKind> = BaseHolder<T,K,K::UsedIndexType>;

impl<T> UniqueHolder<T> {
    // This is equivalent to construction of normal smart pointer
    // and creates an "untyped" holder
    pub fn new(data: T) -> Self {
        Self {
            arc: triomphe::Arc::new(data),
            used: (),
            _kind: Default::default(),
        }
    }
}

// Conversions from "untyped" holder.
// This should not fail because there is no legal way to make
// non-unique untyped holder.
// Just in case we check this and pinic if this is violated.
macro_rules! impl_from_untyped_holder {
    ( $t:ty ) => {
        impl<T> From<UniqueHolder<T>> for Holder<T, $t> {
            fn from(value: UniqueHolder<T>) -> Self {
                if !value.arc.is_unique() {
                    panic!("Untyped Holder is not uniquesly owned!")
                }
                Self {
                    arc: value.arc,
                    used: Default::default(),
                    _kind: Default::default(),
                }
            }
        }
    };
}

impl_from_untyped_holder!(MutableSerial);
impl_from_untyped_holder!(BuilderSerial);
impl_from_untyped_holder!(MutableParallel);
impl_from_untyped_holder!(ImmutableParallel);

// All holders except MutableParallel allow cloning
impl<T,K> Clone for BaseHolder<T,K,()> {
    fn clone(&self) -> Self {
        Self {
            arc: self.arc.clone(),
            used: (),
            _kind: Default::default(),
        }
    }
}

pub(crate) type UsedHashMap = triomphe::Arc<Mutex<rustc_hash::FxHashSet<usize>>>;

impl<T,K: SelectionKind> Holder<T, K> {
    pub(crate) fn clone_with_index(
        &self,
        ind: &impl IndexProvider,
    ) -> Result<Self, SelectionError> {
        K::try_add_used(ind, &self.used)?;
        Ok(Self {
            arc: self.arc.clone(),
            used: self.used.clone(),
            _kind: Default::default(),
        })
    }    
}

// All holders are dereferenced as usual smart pointers
impl<T, K, U> Deref for BaseHolder<T, K, U> {
    type Target = T;
    fn deref(&self) -> &Self::Target {
        &self.arc
    }
}

// All holders can release a UniqueHolder if
// they are uniquesly owned
impl<T, K: SelectionKind> Holder<T, K> {
    pub fn release(self) -> Result<UniqueHolder<T>, SelectionError> {
        if self.arc.is_unique() {
            Ok(UniqueHolder {
                arc: self.arc,
                _kind: Default::default(),
                used: (),
            })
        } else {
            Err(SelectionError::Release)
        }
    }
}

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
        topology: UniqueHolder<Topology>,
        state: UniqueHolder<State>,
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
        topology: UniqueHolder<Topology>,
        state: UniqueHolder<State>,
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
            topology: BaseHolder::new(Topology::default()).into(),
            state: BaseHolder::new(State::default()).into(),            
            _no_send: Default::default(),
            _kind: Default::default(),
        }
    }

    /// Creates a source of mutable parallel selections that _can't_ overlap.
    pub fn new_parallel_mut(
        topology: UniqueHolder<Topology>,
        state: UniqueHolder<State>,
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
        topology: UniqueHolder<Topology>,
        state: UniqueHolder<State>,
    ) -> Result<Source<ImmutableParallel>, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        Ok(Source {
            topology: topology.into(),
            state: state.into(),            
            _no_send: Default::default(),
            _kind: Default::default(),
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

    // Constructor for internal usage
    pub(crate) fn new_from_arc_system(system: triomphe::Arc<System>) -> Result<Source<MutableSerial>, SelectionError> {
        check_topology_state_sizes(&system.topology, &system.state)?;
        Ok(Source {
            system,
            _marker: Default::default(),
        })
    }
}

impl<K: SerialSel> Source<K> {
    /// Release and return [Topology] and [State]. Fails if any selections created from this [Source] are still alive.
    pub fn release(self) -> Result<(UniqueHolder<Topology>, UniqueHolder<State>), SelectionError> {
        Ok((
            self.topology.release()?,
            self.state.release()?,
        ))
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
        Ok(Sel{
            topology: &self.topology.clone_with_index(&vec)?,
            Arc::clone(&self.state),
            vec,
    })
    }

    pub fn select_vec(&mut self, vec: Vec<usize>) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_vec(vec, self.topology.num_atoms())?;
        K::try_add_used(&vec, &mut self.used)?;
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
        K::try_add_used(&vec, &mut self.used)?;
        Ok(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ))
    }

    /// Creates selection of all
    pub fn select_all(&mut self) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_all(self.topology.num_atoms());
        K::try_add_used(&vec, &mut self.used)?;
        Ok(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ))
    }

    /// Creates new selection from a selection expression string. Selection expression is constructed internally but
    /// can't be reused. Consider using [add_expr](Self::add_expr) if you already have selection expression.
    pub fn select_str(&mut self, selstr: impl AsRef<str>) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_str(selstr.as_ref(), &self.topology, &self.state)?;
        K::try_add_used(&vec, &mut self.used)?;
        Ok(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ))
    }

    /// Creates new selection from an existing selection expression.
    pub fn select_expr(&mut self, expr: &SelectionExpr) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_expr(expr, &self.topology, &self.state)?;
        K::try_add_used(&vec, &mut self.used)?;
        Ok(Sel::new_internal(
            Arc::clone(&self.topology),
            Arc::clone(&self.state),
            vec,
        ))
    }

    /// Creates new selection from a range of indexes.
    /// If rangeis out of bounds the error is returned.
    pub fn select_range(
        &mut self,
        range: &std::ops::Range<usize>,
    ) -> Result<Sel<K>, SelectionError> {
        let vec = index_from_range(range, self.topology.num_atoms())?;
        K::try_add_used(&vec, &mut self.used)?;
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
    pub fn set_state(&mut self, state: State) -> Result<State, SelectionError> {
        // Check if the states are compatible
        if !self.system.state.interchangeable(&state) {
            return Err(SelectionError::SetState);
        }
        unsafe {
            ptr::swap(
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
    pub fn set_topology(&mut self, topology: Topology) -> Result<Topology, SelectionError> {
        // Check if the states are compatible
        if !self.system.topology.interchangeable(&topology) {
            return Err(SelectionError::SetTopology);
        }

        unsafe {
            ptr::swap(
                self.system.topology.get_storage_mut(), 
                topology.get_storage_mut()
            );
        };

        Ok(topology)
    }

    pub fn get_shared_topology(&self) -> BaseHolder<Topology, K> {
        BaseHolder::from_arc(Arc::clone(&self.topology))
    }

    pub fn get_shared_state(&self) -> BaseHolder<State, K> {
        BaseHolder::from_arc(Arc::clone(&self.state))
    }

    pub fn new_from_shared(
        topology: BaseHolder<Topology, K>,
        state: BaseHolder<State, K>,
    ) -> Result<Source<MutableSerial>, SelectionError> {
        check_topology_state_sizes(&topology, &state)?;
        Ok(Source {
            topology: topology.into_arc(),
            state: state.into_arc(),
            used: Default::default(),
            _no_send: Default::default(),
            _kind: Default::default(),
        })
    }

    pub fn set_shared_topology(
        &mut self,
        topology: BaseHolder<Topology, K>,
    ) -> Result<BaseHolder<Topology, K>, SelectionError> {
        if !self.topology.interchangeable(&topology) {
            return Err(SelectionError::SetTopology);
        }
        Ok(BaseHolder::from_arc(std::mem::replace(
            &mut self.topology,
            topology.into_arc(),
        )))
    }

    pub fn set_shared_state(
        &mut self,
        state: BaseHolder<State, K>,
    ) -> Result<BaseHolder<State, K>, SelectionError> {
        if !self.state.interchangeable(&state) {
            return Err(SelectionError::SetState);
        }
        Ok(BaseHolder::from_arc(std::mem::replace(
            &mut self.state,
            state.into_arc(),
        )))
    }
}

// Specific methods of builder source
impl Source<BuilderSerial> {    
    pub fn append(&mut self, data: &(impl PosProvider + AtomsProvider)) {
        self.system.topology.get_storage_mut().add_atoms(data.iter_atoms());
        self.system.state.get_storage_mut().add_coords(data.iter_pos());
    }

    pub fn remove(&mut self, to_remove: &impl IndexProvider) -> Result<(), SelectionError> {
        // We are checking index validity inside remove methods
        self.system.topology.get_storage_mut().remove_atoms(to_remove.iter_index())?;
        self.system.state.get_storage_mut().remove_coords(to_remove.iter_index())?;
        Ok(())
    }
}

//------------------
// IO traits impls
//-----------------
impl<K: SelectionKind> TopologyProvider for Source<K> {
    fn num_atoms(&self) -> usize {
        self.system.state.num_coords()
    }
}

impl<K: SelectionKind> AtomsProvider for Source<K> {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        self.system.topology.iter_atoms()       
    }
}

impl<K: SelectionKind> StateProvider for Source<K> {
    fn get_time(&self) -> f32 {
        self.system.state.get_time()
    }

    fn num_coords(&self) -> usize {
        self.system.state.num_coords()       
    }
}

impl<K: SelectionKind> PosProvider for Source<K> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.system.state.iter_pos()
    }
}

impl<K: SelectionKind> BoxProvider for Source<K> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.system.state.get_box()
    }
}

impl<K: SelectionKind> WritableToFile for Source<K> {}
