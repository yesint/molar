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
    }

    pub fn new_fake(n: usize) -> Self {
        Self {
            coords: vec![Pos::origin(); n],
            pbox: None,
            time: 0.0,
        }
    }
}

impl SaveState for State {
    fn iter_pos_dyn<'a>(&'a self) -> Box<dyn ExactSizeIterator<Item = &'a Pos> + 'a> {
        Box::new(self.coords.iter())
    }
}

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

impl LenProvider for State {
    fn len(&self) -> usize {
        self.coords.len()
    }
}

// Identity indexing for State - index i maps to storage position i
impl IndexProvider for State {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        i
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> {
        0..self.coords.len()
    }
}

impl IndexParProvider for State {
    fn par_iter_index(&self) -> impl IndexedParallelIterator<Item = usize> {
        use rayon::iter::IntoParallelIterator;
        (0..self.coords.len()).into_par_iter()
    }
}

impl PosProvider for State {
    unsafe fn pos_unchecked(&self, i: usize) -> &Pos {
        self.coords.get_unchecked(i)
    }
}

impl PosMutProvider for State {
    unsafe fn pos_mut_unchecked(&mut self, i: usize) -> &mut Pos {
        self.coords.get_unchecked_mut(i)
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

impl PosIterProvider for State {
    fn iter_pos(&self) -> impl Iterator<Item = &Pos> {
        self.coords.iter()
    }
}

impl MeasurePos for State {}
impl MeasureRandomAccess for State {}
