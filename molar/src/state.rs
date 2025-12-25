use crate::prelude::*;

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
    pub time: f32,
    pub pbox: Option<PeriodicBox>,
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

impl State {
    pub fn interchangeable(&self, other: &State) -> bool {
        self.coords.len() == other.coords.len()
        // && (
        //     (self.get_storage().pbox.is_none() && other.get_storage().pbox.is_none())
        //     ||
        //     (self.get_storage().pbox.is_some() && other.get_storage().pbox.is_some())
        // )
    }

    pub fn new_fake(n: usize) -> Self {
        Self {
            coords: vec![Pos::origin(); n],
            pbox: None,
            time: 0.0,
        }
    }
}

impl SaveState for State {}

impl TimeProvider for State {
    fn get_time(&self) -> f32 {
        self.time
    }
}

impl TimeMutProvider for State {
    fn set_time(&mut self, t: f32) {
        self.time = t;
    }
}

impl PosIterProvider for State {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.coords.iter()
    }
}

impl LenProvider for State {
    fn len(&self) -> usize {
        self.coords.len()
    }
}

impl RandomPosProvider for State {
    unsafe fn get_pos_unchecked(&self, i: usize) -> &Pos {
        self.coords.get_unchecked(i)
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

impl PosIterMutProvider for State {
    fn iter_pos_mut(&mut self) -> impl PosMutIterator<'_> {
        self.coords.iter_mut()
    }
}

impl RandomPosMutProvider for State {
    fn get_pos_mut(&mut self, i: usize) -> Option<&mut Pos> {
        self.coords.get_mut(i)
    }

    unsafe fn get_pos_mut_unchecked(&mut self, i: usize) -> &mut Pos {
        self.coords.get_unchecked_mut(i)
    }
}

impl MeasurePos for State {}
impl MeasureRandomAccess for State {}
