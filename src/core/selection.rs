use super::{IndexIterator, State, Structure};

#[derive(Debug)]
struct Selection<'a,I: IndexIterator> {
    index_iter: I,
    structure: Option<&'a Structure>,
    state: Option<&'a State>,
}
