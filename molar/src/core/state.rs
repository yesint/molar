use sync_unsafe_cell::SyncUnsafeCell;
use crate::prelude::*;
//use super::handle::{SharedHandle, Handle};


#[doc(hidden)]
#[derive(Debug, Default,Clone)]
pub(crate) struct StateStorage {
    pub coords: Vec<Pos>,
    pub time: f32,
    pub pbox: Option<PeriodicBox>,
}

impl StateStorage {
    pub fn add_coords<'a>(&mut self, pos: impl Iterator<Item = Pos>) {
        self.coords.extend(pos);
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

        // for i in ind.iter().rev().cloned() {
        //     self.coords.remove(i);
        // }

        let mut it = ind.iter().cloned();
        let mut to_remove = it.next().unwrap_or(usize::MAX);
        let mut i = 0;
        self.coords.retain(|_| {
            let ok = i != to_remove;
            i += 1;
            if !ok {
                to_remove = it.next().unwrap_or(usize::MAX);
            }
            ok
        });

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
pub struct State(SyncUnsafeCell<StateStorage>);

impl Clone for State {
    fn clone(&self) -> Self {
        Self(SyncUnsafeCell::new(self.get_storage().clone()))
    }
}

impl From<StateStorage> for State {
    fn from(value: StateStorage) -> Self {
        State(SyncUnsafeCell::new(value))
    }
}

impl State {
    // Private convenience accessors
    #[inline(always)]
    pub(crate) fn get_storage(&self) -> &StateStorage {
        unsafe {&*self.0.get()}
    }

    #[inline(always)]
    pub(crate) fn get_storage_mut(&self) -> &mut StateStorage {
        unsafe {&mut *self.0.get()}
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

    pub fn interchangeable(&self, other: &State) -> bool {
        self.get_storage().coords.len() == other.get_storage().coords.len()
    }

    pub fn new_fake(n: usize) -> Self {
        Self(SyncUnsafeCell::new(StateStorage{
            coords: vec![Pos::origin();n],
            pbox: None,
            time: 0.0,
        }))
    }
}

//------------------------
macro_rules! impl_state_traits {
    ( $t:ty ) => {
        impl StateProvider for $t {
            fn get_time(&self) -> f32 {
                self.get_storage().time
            }
        }
        
        impl PosProvider for $t {
            fn iter_pos(&self) -> impl super::PosIterator<'_> {
                self.get_storage().coords.iter()
            }
        }

        impl RandomPosProvider for $t {
            fn num_coords(&self) -> usize {
                self.get_storage().coords.len()
            }

            unsafe fn nth_pos_unchecked(&self, i: usize) -> &Pos {
                self.get_storage().coords.get_unchecked(i)
            }
        }
        
        impl BoxProvider for $t {
            fn get_box(&self) -> Option<&PeriodicBox> {
                self.get_storage().pbox.as_ref()
            }
        }
        
        impl PosMutProvider for $t {
            fn iter_pos_mut(&self) -> impl super::PosMutIterator<'_> {
                self.get_storage_mut().coords.iter_mut()
            }
        }

        impl LenProvider for $t {
            fn len(&self) -> usize {
                self.get_storage().coords.len()
            }
        }

        impl RandomPosMut for $t {
            fn nth_pos_mut(&self, i: usize) -> Option<&mut Pos> {
                self.get_storage_mut().coords.get_mut(i)
            }

            unsafe fn nth_pos_mut_unchecked(&self, i: usize) -> &mut Pos {
                self.get_storage_mut().coords.get_unchecked_mut(i)
            }
        }

        impl MeasurePos for $t {}
        impl MeasureRandomAccess for $t {}
    }
}

// Impls for State itself
impl_state_traits!(State);
// Impls for smart pointers
impl_state_traits!(Holder<State,MutableSerial>);
impl_state_traits!(Holder<State,BuilderSerial>);
impl_state_traits!(Holder<State,MutableParallel>);
impl_state_traits!(Holder<State,ImmutableParallel>);