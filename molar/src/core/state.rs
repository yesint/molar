use crate::prelude::*;
use sync_unsafe_cell::SyncUnsafeCell;

#[doc(hidden)]
#[derive(Debug, Default, Clone)]
pub(crate) struct StateStorage {
    pub coords: Vec<Pos>,
    pub time: f32,
    pub pbox: Option<PeriodicBox>,
}

impl StateStorage {
    pub fn add_coords<'a>(&mut self, pos: impl Iterator<Item = Pos>) {
        self.coords.extend(pos);
    }

    pub fn remove_coords(
        &mut self,
        removed: impl Iterator<Item = usize>,
    ) -> Result<(), BuilderError> {
        let mut ind = removed.collect::<Vec<_>>();
        if ind.len() == 0 {
            return Ok(());
        }
        ind.sort_unstable();
        ind.dedup();
        if ind[0] >= self.coords.len() || ind[ind.len() - 1] >= self.coords.len() {
            return Err(BuilderError::RemoveIndexes(
                ind[0],
                ind[ind.len() - 1],
                self.coords.len(),
            ));
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

impl std::fmt::Debug for State {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = self.get_storage();
        f.debug_struct("State")
            .field("coords", &s.coords)
            .field("time", &s.time)
            .field("pbox", &s.pbox)
            .finish()
    }
}

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
        unsafe { &*self.0.get() }
    }

    #[inline(always)]
    pub(crate) fn get_storage_mut(&self) -> &mut StateStorage {
        unsafe { &mut *self.0.get() }
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
            // && (
            //     (self.get_storage().pbox.is_none() && other.get_storage().pbox.is_none())
            //     || 
            //     (self.get_storage().pbox.is_some() && other.get_storage().pbox.is_some())
            // )
    }

    pub fn new_fake(n: usize) -> Self {
        Self(SyncUnsafeCell::new(StateStorage {
            coords: vec![Pos::origin(); n],
            pbox: None,
            time: 0.0,
        }))
    }
}

//------------------------
macro_rules! impl_state_traits {
    ( $t:ty ) => {
        impl StateIoProvider for $t {}
        
        impl TimeProvider for $t {
            fn get_time(&self) -> f32 {
                self.get_storage().time
            }
        }

        impl TimeMutProvider for $t {
            fn set_time(&self, t: f32) {
                self.get_storage_mut().time = t;
            }
        }

        impl PosIterProvider for $t {
            fn iter_pos(&self) -> impl super::PosIterator<'_> + Clone {
                self.get_storage().coords.iter()
            }
        }

        impl LenProvider for $t {
            fn len(&self) -> usize {
                self.get_storage().coords.len()
            }
        }

        impl RandomPosProvider for $t {
            unsafe fn get_pos_unchecked(&self, i: usize) -> &Pos {
                self.get_storage().coords.get_unchecked(i)
            }
        }

        impl BoxProvider for $t {
            fn get_box(&self) -> Option<&PeriodicBox> {
                self.get_storage().pbox.as_ref()
            }
        }

        impl PosIterMutProvider for $t {
            fn iter_pos_mut(&self) -> impl super::PosMutIterator<'_> {
                self.get_storage_mut().coords.iter_mut()
            }
        }

        impl RandomPosMutProvider for $t {
            fn get_pos_mut(&self, i: usize) -> Option<&mut Pos> {
                self.get_storage_mut().coords.get_mut(i)
            }

            unsafe fn get_pos_mut_unchecked(&self, i: usize) -> &mut Pos {
                self.get_storage_mut().coords.get_unchecked_mut(i)
            }
        }

        impl MeasurePos for $t {}
        impl MeasureRandomAccess for $t {}
    };
}

// Impls for State itself
impl_state_traits!(State);
// Impls for smart pointers
impl_state_traits!(triomphe::Arc<State>);
