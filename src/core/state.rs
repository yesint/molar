use std::ops::Deref;
use sync_unsafe_cell::SyncUnsafeCell;

use crate::io::StateProvider;
use super::{providers::{BoxProvider, PosProvider}, PeriodicBox, Pos};
//use super::handle::{SharedHandle, Handle};


#[doc(hidden)]
#[derive(Debug, Default,Clone)]
pub struct StateStorage {
    pub coords: Vec<Pos>,
    pub time: f32,
    pub pbox: Option<PeriodicBox>,
}

/// State of molecular system including its coordinates, time stamp 
/// and [periodic box](super::PeriodicBox).
/// 
/// [State] is typically read from structure of trajectory file and is not intended
/// to be manipulated directly by the user. Insead [State] and [Topology](super::Topology)
/// are used to create atom selections, which give an access to the properties of
/// individual atoms and allow to query various properties.
pub struct State(SyncUnsafeCell<StateStorage>);

impl Clone for State {
    fn clone(&self) -> Self {
        Self(SyncUnsafeCell::new(self.get().clone()))
    }
}

impl From<StateStorage> for State {
    fn from(value: StateStorage) -> Self {
        Self(SyncUnsafeCell::new(value))
    }
}

impl State {
    // Private convenience accessors
    #[inline(always)]
    fn get(&self) -> &StateStorage {
        unsafe {&*self.0.get()}
    }

    #[inline(always)]
    fn get_mut(&self) -> &mut StateStorage {
        unsafe {&mut *self.0.get()}
    }

    //-----------------------------------------
    pub fn to_rc(self) -> triomphe::UniqueArc<Self> {
        triomphe::UniqueArc::new(self)
    }

    #[inline(always)]
    pub unsafe fn nth_pos_unchecked(&self, i: usize) -> &Pos {
        self.get().coords.get_unchecked(i)
    }

    #[inline(always)]
    pub unsafe fn nth_pos_unchecked_mut(&self, i: usize) -> &mut Pos {
        self.get_mut().coords.get_unchecked_mut(i)
    }

    #[inline(always)]
    pub fn nth_pos(&self, i: usize) -> Option<&Pos> {
        self.get().coords.get(i)
    }

    #[inline(always)]
    pub fn nth_pos_mut(&self, i: usize) -> Option<&mut Pos> {
        self.get_mut().coords.get_mut(i)
    }

    #[inline(always)]
    pub fn get_box(&self) -> Option<&PeriodicBox> {
        self.get().pbox.as_ref()
    }

    #[inline(always)]
    pub fn get_box_mut(&self) -> Option<&mut PeriodicBox> {
        self.get_mut().pbox.as_mut()
    }

    pub fn interchangeable(&self, other: &State) -> bool {
        self.get().coords.len() == other.get().coords.len()
    }
    
}

// Impls for smart pointers
impl<T: Deref<Target=State>> StateProvider for T {
    fn get_time(&self) -> f32 {
        self.get().time
    }

    fn num_coords(&self) -> usize {
        self.get().coords.len()
    }
}

impl<T: Deref<Target=State>> PosProvider for T {
    fn iter_pos(&self) -> impl super::PosIterator<'_> {
        self.get().coords.iter()
    }
}

impl<T: Deref<Target=State>> BoxProvider for T {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.get().pbox.as_ref()
    }
}

// Impls for State itself
impl StateProvider for State {
    fn get_time(&self) -> f32 {
        self.get().time
    }

    fn num_coords(&self) -> usize {
        self.get().coords.len()
    }
}

impl PosProvider for State {
    fn iter_pos(&self) -> impl super::PosIterator<'_> {
        self.get().coords.iter()
    }
}

impl BoxProvider for State {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.get().pbox.as_ref()
    }
}
