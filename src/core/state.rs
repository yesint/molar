use crate::prelude::*;
use std::marker::PhantomData;
use sync_unsafe_cell::SyncUnsafeCell;
use triomphe::Arc;

#[doc(hidden)]
#[derive(Debug, Default,Clone)]
pub(crate) struct StateStorage {
    pub coords: Vec<Pos>,
    pub time: f32,
    pub pbox: Option<PeriodicBox>,
}
impl StateStorage {
    pub fn add_coords<'a>(&mut self, pos: impl super::PosIterator<'a>) {
        self.coords.extend(pos.cloned());
    }

    pub fn remove_coords(&mut self, removed: impl Iterator<Item = usize>) -> Result<(),BuilderError> {
        let mut ind = removed.collect::<Vec<_>>();
        if ind.len()==0 {
            return Ok(());
        }
        ind.sort_unstable();
        ind.dedup();
        if ind[0] >= self.coords.len() || ind[ind.len()-1] >= self.coords.len() {
            return Err(BuilderError::RemoveIndexes(ind[0],ind[ind.len()-1],self.coords.len()));
        }

        for i in ind.iter().rev().cloned() {
            self.coords.remove(i);
        }
        Ok(())
    }
}

/// State of molecular system including its coordinates, time stamp 
/// and [periodic box](super::PeriodicBox).
/// 
/// [State] is typically read from structure of trajectory file and is not intended
/// to be manipulated directly by the user. Insead [State] and [Topology](super::Topology)
/// are used to create atom selections, which give an access to the properties of
/// individual atoms and allow to query various properties.
#[derive(Default)]
pub struct State<K=()>{
    pub(crate) arc: Arc<SyncUnsafeCell<StateStorage>>,
    _marker: PhantomData<K>,
}

impl<K> Clone for State<K> {
    fn clone(&self) -> Self {
        Self{
            arc: Arc::new(SyncUnsafeCell::new(self.get_storage().clone())),
            _marker: Default::default(),
        }
    }
}

fn convert<K1,K2>(value: State<K1>) -> Result<State<K2>, SelectionError> {
    if Arc::count(&value.arc) == 1 {
        Ok(State::<K2> {
            arc: value.arc,
            _marker: Default::default(),
        })
    } else {
        Err(SelectionError::Release)
    }
}

impl TryFrom<State<()>> for State<MutableSerial> {
    type Error = SelectionError;
    fn try_from(value: State<()>) -> Result<Self, Self::Error> {
        convert(value)
    }
}

impl TryFrom<State<()>> for State<ImmutableParallel> {
    type Error = SelectionError;
    fn try_from(value: State<()>) -> Result<Self, Self::Error> {
        convert(value)
    }
}

impl TryFrom<State<()>> for State<MutableParallel> {
    type Error = SelectionError;
    fn try_from(value: State<()>) -> Result<Self, Self::Error> {
        convert(value)
    }
}

impl TryFrom<State<()>> for State<BuilderSerial> {
    type Error = SelectionError;
    fn try_from(value: State<()>) -> Result<Self, Self::Error> {
        convert(value)
    }
}

impl From<StateStorage> for State<()> {
    fn from(value: StateStorage) -> Self {
        State {
            arc: Arc::new(SyncUnsafeCell::new(value)),
            _marker: Default::default(),
        }    
    }
}

impl<K> State<K> {
    pub fn new_arc<K2: SelectionKind>(&self) -> State<K2> {
        State::<K2>{
            arc: Arc::clone(&self.arc),
            _marker: Default::default(),
        }
    }

    // Private convenience accessors
    #[inline(always)]
    pub(crate) fn get_storage(&self) -> &StateStorage {
        unsafe {&*self.arc.get()}
    }

    #[inline(always)]
    pub(crate) fn get_storage_mut(&self) -> &mut StateStorage {
        unsafe {&mut *self.arc.get()}
    }
    
    #[inline(always)]
    pub unsafe fn nth_pos_unchecked(&self, i: usize) -> &Pos {
        self.get_storage().coords.get_unchecked(i)
    }

    #[inline(always)]
    pub unsafe fn nth_pos_unchecked_mut(&self, i: usize) -> &mut Pos {
        self.get_storage_mut().coords.get_unchecked_mut(i)
    }

    #[inline(always)]
    pub fn nth_pos(&self, i: usize) -> Option<&Pos> {
        self.get_storage().coords.get(i)
    }

    #[inline(always)]
    pub fn nth_pos_mut(&self, i: usize) -> Option<&mut Pos> {
        self.get_storage_mut().coords.get_mut(i)
    }

    #[inline(always)]
    pub fn get_box(&self) -> Option<&PeriodicBox> {
        self.get_storage().pbox.as_ref()
    }

    #[inline(always)]
    pub fn get_box_mut(&self) -> Option<&mut PeriodicBox> {
        self.get_storage_mut().pbox.as_mut()
    }

    pub fn interchangeable(&self, other: &State<K>) -> bool {
        self.get_storage().coords.len() == other.get_storage().coords.len()
    }
}

// Impls for State itself
impl<K> StateProvider for State<K> {
    fn get_time(&self) -> f32 {
        self.get_storage().time
    }

    fn num_coords(&self) -> usize {
        self.get_storage().coords.len()
    }
}

impl<K> PosProvider for State<K> {
    fn iter_pos(&self) -> impl super::PosIterator<'_> {
        self.get_storage().coords.iter()
    }
}

impl<K> BoxProvider for State<K> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.get_storage().pbox.as_ref()
    }
}
