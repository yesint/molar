use thiserror::Error;
use crate::prelude::*;

/// Error type for missing optional state data (velocities or forces)
#[derive(Error, Debug)]
pub enum StateError {
    #[error("velocities are not present in the current state")]
    NoVelocities,
    #[error("forces are not present in the current state")]
    NoForces,
}

/// State of molecular system including its coordinates, time stamp
/// and [periodic box](super::PeriodicBox).
///
/// [State] is typically read from structure of trajectory file and is not intended
/// to be manipulated directly by the user. Insead [State] and [Topology](super::Topology)
/// are used to create atom selections, which give an access to the properties of
/// individual atoms and allow to query various properties.

#[derive(Debug, Default, Clone)]
pub struct State {
    pub coords: Vec<Pos>,
    pub velocities: Vec<Vel>,
    pub forces: Vec<Force>,
    pub time: Float,
    pub pbox: Option<PeriodicBox>,
}

impl State {
    pub fn has_pos(&self) -> bool { !self.coords.is_empty() }
    pub fn has_vel(&self) -> bool { !self.velocities.is_empty() }
    pub fn has_force(&self) -> bool { !self.forces.is_empty() }
}

impl State {
    pub fn add_coords(&mut self, pos: impl Iterator<Item = Pos>) {
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

        macro_rules! retain_not_in {
            ($vec:expr) => {
                if !$vec.is_empty() {
                    let mut it = ind.iter().cloned();
                    let mut to_remove = it.next().unwrap_or(usize::MAX);
                    let mut i = 0usize;
                    $vec.retain(|_| {
                        let ok = i != to_remove;
                        i += 1;
                        if !ok { to_remove = it.next().unwrap_or(usize::MAX); }
                        ok
                    });
                }
            };
        }

        retain_not_in!(self.coords);
        retain_not_in!(self.velocities);
        retain_not_in!(self.forces);

        Ok(())
    }
}

impl State {
    pub fn interchangeable(&self, other: &State) -> bool {
        self.coords.len() == other.coords.len()
    }

    pub fn new_fake(n: usize) -> Self {
        Self {
            coords: vec![Pos::origin(); n],
            velocities: Vec::new(),
            forces: Vec::new(),
            pbox: None,
            time: 0.0,
        }
    }
}

impl SaveState for State {
    fn iter_pos_dyn<'a>(&'a self) -> Box<dyn ExactSizeIterator<Item = &'a Pos> + 'a> {
        Box::new(self.coords.iter())
    }

    fn iter_vel_dyn<'a>(&'a self) -> Box<dyn ExactSizeIterator<Item = &'a Vel> + 'a> {
        Box::new(self.velocities.iter())
    }

    fn iter_force_dyn<'a>(&'a self) -> Box<dyn ExactSizeIterator<Item = &'a Force> + 'a> {
        Box::new(self.forces.iter())
    }
}

impl TimeProvider for State {
    fn get_time(&self) -> Float {
        self.time
    }
}

impl TimeMutProvider for State {
    fn set_time(&mut self, t: Float) {
        self.time = t;
    }
}

impl LenProvider for State {
    fn len(&self) -> usize {
        let n = self.coords.len();
        if n > 0 { return n; }
        let n = self.velocities.len();
        if n > 0 { return n; }
        self.forces.len()
    }
}

/// Identity index provider for State (index i → position i)
impl IndexProvider for State {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        i
    }

    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize> {
        0..self.coords.len()
    }
}

impl PosProvider for State {
    unsafe fn coords_ptr(&self) -> *const Pos {
        self.coords.as_ptr()
    }
}

impl PosMutProvider for State {
    unsafe fn coords_ptr_mut(&mut self) -> *mut Pos {
        self.coords.as_mut_ptr()
    }
}

impl BoxProvider for State {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.pbox.as_ref()
    }
}

impl BoxMutProvider for State {
    fn get_box_mut(&mut self) -> Option<&mut PeriodicBox> {
        self.pbox.as_mut()
    }
}

impl VelProvider for State {
    unsafe fn vel_ptr(&self) -> *const Vel {
        if self.velocities.is_empty() { std::ptr::null() } else { self.velocities.as_ptr() }
    }
}

impl VelMutProvider for State {
    unsafe fn vel_ptr_mut(&mut self) -> *mut Vel {
        if self.velocities.is_empty() { std::ptr::null_mut() } else { self.velocities.as_mut_ptr() }
    }
}

impl ForceProvider for State {
    unsafe fn force_ptr(&self) -> *const Force {
        if self.forces.is_empty() { std::ptr::null() } else { self.forces.as_ptr() }
    }
}

impl ForceMutProvider for State {
    unsafe fn force_ptr_mut(&mut self) -> *mut Force {
        if self.forces.is_empty() { std::ptr::null_mut() } else { self.forces.as_mut_ptr() }
    }
}

// MeasurePos and MeasureRandomAccess are provided by blanket impls in traits.rs
