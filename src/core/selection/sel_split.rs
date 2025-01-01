use std::marker::PhantomData;

use rayon::iter::IntoParallelRefMutIterator;
use crate::prelude::*;

//-------------------------------------------------------
// Splitting iterators
//-------------------------------------------------------

struct SplitData<RT, F> {
    func: F,
    counter: usize,
    val: RT,
}

/// Iterator over contiguous pieces of selection returned by [Sel::split_contig]
/// and `various Sel::split_contig_*` convenience methods.
///
/// This iterator keeps the parent selection alive and yelds selections of the same kind
/// as sub-selections of the parent [Sel].
pub struct FragmentsIterator<'a, K, I, RT> 
where
    K: SelectionKind,
{
    sel: &'a Sel<K>,
    iter: I,
    cur: (usize,RT),
}

impl<'a,K: SelectionKind> FragmentsIterator<'a, K, (), ()> {
    pub fn new<RT, F>(sel: &'a Sel<K>, func: F) -> Result<FragmentsIterator<'a, K, impl Iterator<Item = (usize, RT)> + 'a, RT>, SelectionError>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(Particle) -> Option<RT> + 'a,
    {
        // Iterator oved valid particles for which split_n returns Some
        let mut iter = sel.iter_particle().filter_map(move |p| {
            let i = p.id;
            func(p).map(|val| (i, val))
        });

        let cur = iter.next().ok_or_else(|| SelectionError::EmptySplit)?;

        Ok(FragmentsIterator {
            sel,
            iter,
            cur,
        })
    }
}

impl<K, I, RT> Iterator for FragmentsIterator<'_, K, I, RT>
where
    K: UserCreatableKind,
    I: Iterator<Item = (usize, RT)>,
    RT: PartialEq,
{
    type Item = Sel<K>;

    fn next(&mut self) -> Option<Self::Item> {
        // Iterate until we get something different or the end (indicated by None)
        let e = self.iter.find_map(|el| if el.1 != self.cur.1 { Some(el) } else { None });
        let bi = self.cur.0;
        // Create selection (b:e) or (b:end)
        if let Some(e) = e {
            let ei = e.0;
            self.cur = e;
            Some(unsafe{
                self.sel.subsel_from_sorted_vec_unchecked((bi..ei).collect()).unwrap()
            })
        } else if bi < self.sel.last_index() {
            Some(unsafe{
                self.cur.0 = self.sel.last_index(); // Indicate end
                self.sel.subsel_from_sorted_vec_unchecked((bi..self.sel.last_index()).collect()).unwrap()
            })
        } else {
            None
        }
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
                val: RT::default(),
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
        let val = (data.func)(p);

        if let Some(val) = val {
            if val == data.val {
                // Current selection continues. Add current index
                index.push(i);
            } else if index.is_empty() {
                // The very first id is not default, this is Ok, add index
                // and update self.id
                data.val = val;
                index.push(i);
            } else {
                // The end of current selection
                data.val = val; // Update self.id for the next selection
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


//===============================
// Container for parallel split
//===============================

pub struct ParallelSplit {
    pub(super) parts: Vec<Sel<MutableParallel>>,
    pub(super) _marker: PhantomData<*const ()>,
}

impl ParallelSplit {
    pub fn par_iter(&mut self) -> rayon::slice::IterMut<'_, Sel<MutableParallel>> {
        self.parts.par_iter_mut()
    }
}
