use sorted_vec::SortedSet;

use crate::prelude::*;

//-------------------------------------------------------
// Splitting iterators
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
pub struct SelectionSplitIterator<'a, RT, F, K: SelectionKind> {
    sel: &'a Sel<K>,
    data: SplitData<RT, F>,
}

impl<K: SelectionKind> SelectionSplitIterator<'_, (), (), K> {
    pub fn new<RT, F>(sel: &Sel<K>, func: F) -> SelectionSplitIterator<'_, RT, F, K>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(Particle) -> Option<RT>,
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
    F: Fn(Particle) -> Option<RT>,
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
pub struct IntoFragmentsIterator<RT, F, K: SelectionKind> {
    sel: Sel<K>,
    data: SplitData<RT, F>,
}

impl<RT,F,K: SelectionKind> Drop for IntoFragmentsIterator<RT, F, K> {
    fn drop(&mut self) {
        // If stored selection is MutableParallel it will clear used indexes
        // when dropped. This will invalidate used indexes because they are already
        // used by the fragments created by iterator.
        // To avoid this we clear index so that nothing is cleared upon dropping selection.
        unsafe {
            self.sel.clear_index_before_drop();
        }
        // self.sel is now Ok to drop
    }
}

impl<K: SelectionKind> IntoFragmentsIterator<(), (), K> {
    pub fn new<RT, F>(sel: Sel<K>, func: F) -> IntoFragmentsIterator<RT, F, K>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(Particle) -> Option<RT>,
    {
        IntoFragmentsIterator {
            sel,
            data: SplitData {
                func,
                counter: 0,
                id: RT::default(),
            },
        }
    }
}

impl<RT, F, S> Iterator for IntoFragmentsIterator<RT, F, S>
where
    RT: Default + std::cmp::PartialEq,
    F: Fn(Particle) -> Option<RT>,
    S: SelectionKind,
{
    // Consuming splitter always return the same selection kind as parent
    type Item = Sel<S>;

    fn next(&mut self) -> Option<Self::Item> {
        next_split(&mut self.data, &self.sel)
    }
}

// Actual function that does the splitting
fn next_split<RT, F, K>(data: &mut SplitData<RT, F>, sel: &Sel<K>) -> Option<Sel<K>>
where
    RT: Default + std::cmp::PartialEq,
    F: Fn(Particle) -> Option<RT>,
    K: SelectionKind,
{
    let mut index = Vec::<usize>::new();
    while data.counter < sel.len() {
        let p = unsafe { sel.nth_particle_unchecked(data.counter) };
        let i = p.id;
        let id = (data.func)(p);

        if let Some(id) = id {
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
                return unsafe { Some(sel.subsel_from_sorted_vec_unchecked(index).unwrap()) };
            }
        }
        // Next element
        data.counter += 1;
    }

    // Return any remaining index as last selection
    if !index.is_empty() {
        return unsafe { Some(sel.subsel_from_sorted_vec_unchecked(index).unwrap()) };
    }

    // If we are here stop iterating
    None
}