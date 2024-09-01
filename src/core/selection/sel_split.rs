use crate::prelude::*;

//-------------------------------------------------------
// Splitting iterator
//-------------------------------------------------------

struct SplitData<RT, F> {
    func: F,
    counter: usize,
    id: RT,
}

/// Iterator over contiguous pieces of selection returned by [Sel::split_contig] 
/// and `various Sel::split_contig_*` convenience methods.
/// 
/// This iterator keeps the parent selection alive and yelds selections of the same kind 
/// as sub-selections of the parent [Sel].
pub struct SelectionSplitIterator<'a, RT, F, S> {
    sel: &'a Sel<S>,
    data: SplitData<RT, F>,
}

impl<S> SelectionSplitIterator<'_, (), (), S> {
    pub fn new<RT, F>(sel: &Sel<S>, func: F) -> SelectionSplitIterator<'_, RT, F, S>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(Particle) -> RT,
    {
        SelectionSplitIterator {
            sel,
            data: SplitData {
                func,
                counter: 0,
                id: RT::default(),
            },
        }
    }
}

impl<RT, F, S> Iterator for SelectionSplitIterator<'_, RT, F, S>
where
    RT: Default + std::cmp::PartialEq,
    F: Fn(Particle) -> RT,
    S: SelectionKind,
{
    type Item = Sel<S>;

    fn next(&mut self) -> Option<Self::Item> {
        next_split(&mut self.data, self.sel)
    }
}

//--------------------------------------------
/// Iterator over contiguous pieces of selection returned by [Sel::into_split_contig] 
/// and `various Sel::into_split_contig_*` convenience methods.
/// 
/// This iterator consumes the parent selection and yelds selections of the same kind 
/// as the parent [Sel].
pub struct IntoSelectionSplitIterator<RT, F, S> {
    sel: Sel<S>,
    data: SplitData<RT, F>,
}

impl<S> IntoSelectionSplitIterator<(), (), S> {
    pub fn new<RT, F>(sel: Sel<S>, func: F) -> IntoSelectionSplitIterator<RT, F, S>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(usize, &Atom, &Pos) -> RT,
    {
        IntoSelectionSplitIterator {
            sel,
            data: SplitData {
                func,
                counter: 0,
                id: RT::default(),
            },
        }
    }
}

impl<RT, F, S> Iterator for IntoSelectionSplitIterator<RT, F, S>
where
    RT: Default + std::cmp::PartialEq,
    F: Fn(Particle) -> RT,
    S: SelectionKind,
{
    // Consuming splitter always return the same selection kind as parent
    type Item = Sel<S>;

    fn next(&mut self) -> Option<Self::Item> {
        next_split(&mut self.data, &self.sel)
    }
}

// Actual function that does the splitting
fn next_split<RT, F, S, SR>(data: &mut SplitData<RT, F>, sel: &Sel<S>) -> Option<Sel<SR>>
where
    RT: Default + std::cmp::PartialEq,
    F: Fn(Particle) -> RT,
    S: SelectionKind,
    SR: SelectionKind,
{
    let mut index = Vec::<usize>::new();
    while data.counter < sel.len() {
        let p = unsafe { sel.nth_particle_unchecked(data.counter) };
        let i = p.id;
        let id = (data.func)(p);

        if id == data.id {
            // Current selection continues. Add current index
            index.push(i);
        } else if index.is_empty() {
            // The very first id is not default, this is Ok, add index
            // and update self.id
            data.id = id;
            index.push(i);
        } else {
            // The end of current selection
            data.id = id; // Update self.id for the next selection
            return unsafe { Some(sel.subsel_from_vec_unchecked(index).unwrap()) };
        }
        // Next element
        data.counter += 1;
    }

    // Return any remaining index as last selection
    if !index.is_empty() {
        return unsafe { Some(sel.subsel_from_vec_unchecked(index).unwrap()) };
    }

    // If we are here stop iterating
    None
}