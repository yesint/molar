use std::sync::{RwLock, RwLockReadGuard, RwLockWriteGuard};

use crate::prelude::*;

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
pub struct State(RwLock<StateStorage>);

pub struct StateReadGuard<'a>(pub(crate) RwLockReadGuard<'a, StateStorage>);
pub struct StateWriteGuard<'a>(pub(crate) RwLockWriteGuard<'a, StateStorage>);

impl State {
    pub fn read(&self) -> StateReadGuard<'_> {
        StateReadGuard(self.0.read().unwrap())
    }

    pub fn write(&self) -> StateWriteGuard<'_> {
        StateWriteGuard(self.0.write().unwrap())
    }

    #[inline(always)]
    pub fn interchangeable(&self, other: &State) -> bool {
        self.read().0.coords.len() == other.read().0.coords.len()
        // && (
        //     (self.get_storage().pbox.is_none() && other.get_storage().pbox.is_none())
        //     ||
        //     (self.get_storage().pbox.is_some() && other.get_storage().pbox.is_some())
        // )
    }

    pub fn new_fake(n: usize) -> Self {
        Self(RwLock::new(StateStorage {
            coords: vec![Pos::origin(); n],
            pbox: None,
            time: 0.0,
        }))
    }
}

//------------------------
macro_rules! impl_state_traits {
    ( $t:ty ) => {
        //impl StateWrite for $t {}

        impl TimeProvider for $t {
            fn get_time(&self) -> f32 {
                self.0.time
            }
        }

        impl PosIterProvider for $t {
            fn iter_pos(&self) -> impl super::PosIterator<'_> + Clone {
                self.0.coords.iter()
            }
        }

        impl LenProvider for $t {
            fn len(&self) -> usize {
                self.0.coords.len()
            }
        }

        impl RandomPosProvider for $t {
            unsafe fn get_pos_unchecked(&self, i: usize) -> &Pos {
                self.0.coords.get_unchecked(i)
            }
        }

        impl BoxProvider for $t {
            fn get_box(&self) -> Option<&PeriodicBox> {
                self.0.pbox.as_ref()
            }
        }

        impl MeasurePos for $t {}
        impl MeasureRandomAccess for $t {}
    };
}

macro_rules! impl_state_mut_traits {
    ( $t:ty ) => {
        impl TimeMutProvider for $t {
            fn set_time(&mut self, t: f32) {
                self.0.time = t;
            }
        }
        
        impl PosIterMutProvider for $t {
            fn iter_pos_mut(&mut self) -> impl super::PosMutIterator<'_> {
                self.0.coords.iter_mut()
            }
        }

        impl RandomPosMutProvider for $t {
            fn get_pos_mut(&mut self, i: usize) -> Option<&mut Pos> {
                self.0.coords.get_mut(i)
            }

            unsafe fn get_pos_mut_unchecked(&mut self, i: usize) -> &mut Pos {
                self.0.coords.get_unchecked_mut(i)
            }
        }
    };
}

// Impls for State itself
impl_state_traits!(StateReadGuard<'_>);
impl_state_traits!(StateWriteGuard<'_>);
impl_state_mut_traits!(StateWriteGuard<'_>);
// Impls for smart pointers
//impl_state_traits!(triomphe::Arc<State>);
