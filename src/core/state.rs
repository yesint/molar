use super::{IndexIterator, PeriodicBox, Pos, PosIterator};
use nalgebra::Point3;

#[derive(Debug, Default)]
pub struct State {
    pub coords: Vec<Point3<f32>>,
    pub time: f32,
    pub box_: Option<PeriodicBox>,
}

impl State {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn subset<T>(&self, index: T) -> SubState<T>
    where
        T: IndexIterator,
    {
        SubState::new(self, index)
    }
}

//======================================================

pub trait IdPosIterator<'a>: ExactSizeIterator<Item = (usize, &'a Pos)> {}
impl<'a, T> IdPosIterator<'a> for T where T: ExactSizeIterator<Item = (usize, &'a Pos)> {}

//======================================================

/// Substate encapsulating state and an iterator over selected indexes
/// It doesn't store indexes internally thus not indexible
pub struct SubState<'a, T>
where
    T: IndexIterator,
{
    state: &'a State,
    index_iter: T,
    pub time: &'a f32,
    pub box_: &'a Option<PeriodicBox>,
}

impl<'a, T> SubState<'a, T>
where
    T: IndexIterator,
{
    fn new(state: &'a State, index_iter: T) -> Self {
        Self {
            state,
            index_iter,
            box_: &state.box_,
            time: &state.time,
        }
    }

    // Iterator over selected positions
    fn iter_pos(&self) -> impl PosIterator<'a> {
        self.index_iter.clone().map(|i| &self.state.coords[i])
    }

    // Iterator over selected positions with their indexes
    fn iter_index_pos(&self) -> impl IdPosIterator<'a> {
        self.index_iter.clone().map(|i| (i,&self.state.coords[i]))
    }
}

/// Indexible substate encapsulating state and an iterator over selected indexes
/// It stores indexes internally in a vector taking extra memory
/// in exchange of possibility of indexing.
pub struct SubStateIndexible<'a> {
    state: &'a State,
    index: Vec<usize>,
    pub time: &'a f32,
    pub box_: &'a Option<PeriodicBox>,
}

impl<'a> SubStateIndexible<'a> {
    fn new(state: &'a State, index_iter: impl IndexIterator) -> Self {
        Self {
            state,
            index: index_iter.collect(),
            box_: &state.box_,
            time: &state.time,
        }
    }

    // Iterator over selected positions
    fn iter_pos(&self) -> impl PosIterator<'a> + '_ {
        self.index.iter().cloned().map(|i| &self.state.coords[i])
    }

    // Iterator over selected positions with their indexes
    fn iter_index_pos(&self) -> impl IdPosIterator<'a> + '_ {
        self.index.iter().cloned().map(|i| (i,&self.state.coords[i]))
    }
}

impl<'a> std::ops::Index<usize> for SubStateIndexible<'a> {
    type Output = Pos;
    fn index(&self, i: usize) -> &Self::Output {
        &self.state.coords[self.index[i]]
    }
}

