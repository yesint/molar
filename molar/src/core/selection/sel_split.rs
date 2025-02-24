use std::marker::PhantomData;

use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator};
use crate::prelude::*;

//-------------------------------------------------------
// Splitting iterators
//-------------------------------------------------------

pub(super) struct SplitData<RT, F> {
    pub(super) func: F,
    pub(super) cur: usize,
    pub(super) val: RT,
}

/// Iterator over contiguous pieces of selection.
/// This iterator keeps the parent selection alive.
pub struct FragmentsIterator<'a, RT, F, K> {
    sel: &'a Sel<K>,
    data: SplitData<RT, F>,
}

impl<K: UserCreatableKind> FragmentsIterator<'_, (), (), K> {
    pub fn new<RT, F,>(sel: &Sel<K>, func: F) -> FragmentsIterator<'_, RT, F, K>
    where
        RT: Default + std::cmp::PartialEq,
        F: Fn(Particle) -> Option<RT>,
    {
        FragmentsIterator {
            sel,
            data: SplitData {
                func,
                cur: 0,
                val: RT::default(),
            },
        }
    }
}

impl<RT, F, S> Iterator for FragmentsIterator<'_, RT, F, S>
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

//------------------------------------------------------------
/// Iterator over contiguous pieces of selection.
/// This iterator consumes the parent selection.
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
                cur: 0,
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
pub(super) fn next_split<RT, F, K, KO>(data: &mut SplitData<RT, F>, sel: &Sel<K>) -> Option<Sel<KO>>
where
    RT: Default + std::cmp::PartialEq,
    F: Fn(Particle) -> Option<RT>,
    K: SelectionKind,
    KO: SelectionKind,
{
    let mut index = Vec::<usize>::new();
    while data.cur < sel.len() {
        let p = unsafe { sel.nth_particle_unchecked(data.cur) };
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
        data.cur += 1;
    }

    // Return any remaining index as last selection
    if !index.is_empty() {
        return unsafe { Some(sel.subsel_from_sorted_vec_unchecked(index).unwrap()) };
    }

    // If we are here stop iterating
    None
}

/// Opaque container for parallel selections produced 
/// by `par_collect` and `par_collect_contig` methods.
pub struct ParallelSplit {
    pub(super) parts: Vec<Sel<MutableParallel>>,
    pub(super) _marker: PhantomData<*const ()>,
}

impl ParallelSplit {
    /// Returns parallel iterator over stored parallel selections.
    pub fn par_iter(&mut self) -> rayon::slice::Iter<'_, Sel<MutableParallel>> {
        self.parts.par_iter()
    }

    /// Returns parallel mutable iterator over stored parallel selections.
    pub fn par_iter_mut(&mut self) -> rayon::slice::IterMut<'_, Sel<MutableParallel>> {
        self.parts.par_iter_mut()
    }
}

/// Iterator over molecules
pub struct MoleculesIterator<'a,K> {
    pub(super) sel: &'a Sel<K>,
    pub(super) cur: usize,
}

impl<'a,S: SelectionKind> Iterator for MoleculesIterator<'a,S> {
    type Item = Sel<S>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.cur < self.sel.num_molecules() {
            let [i,j] = unsafe{self.sel.nth_molecule_unchecked(self.cur).to_owned()};
            let index = (i..=j).collect();
            self.cur += 1;
            unsafe{ Some(self.sel.subsel_from_sorted_vec_unchecked(index).unwrap()) }
        } else {
            None
        }
    }
}