use std::ops::Range;
use sorted_vec::SortedSet;
use crate::prelude::*;

pub(crate) fn check_topology_state_sizes(topology: &Topology, state: &State) -> Result<(),super::TopologyStateSizes> {
    let n1 = topology.num_atoms();
    let n2 = state.num_coords();
    if n1 != n2 { Err(super::TopologyStateSizes(n1,n2))? }
    Ok(())
}

pub(super) fn index_from_all(n: usize) -> SortedSet<usize> {
    unsafe { SortedSet::from_sorted((0..n).collect()) }
}

pub(super) fn index_from_expr(expr: &mut SelectionExpr, topology: &Topology, state: &State) -> Result<SortedSet<usize>, SelectionError> {
    let index = expr.apply_whole(&topology, &state)?;
    if index.len() > 0 {
        Ok(index)
    } else {
        Err(SelectionError::FromExpr {
            expr_str: expr.get_str().into(), 
            source: SelectionIndexError::IndexEmpty, 
        })
    }
}

pub(super) fn index_from_expr_sub(expr: &mut SelectionExpr, topology: &Topology, state: &State, subset: &SortedSet<usize>) -> Result<SortedSet<usize>, SelectionError> {
    let index = expr.apply_subset(&topology, &state, subset.iter().cloned())?;
    if index.len() > 0 {
        Ok(index)
    } else {
        Err(SelectionError::FromExpr {
            expr_str: expr.get_str().into(), 
            source: SelectionIndexError::IndexEmpty, 
        })
    }
}

pub(super) fn index_from_str(selstr: &str, topology: &Topology, state: &State) -> Result<SortedSet<usize>, SelectionError> {
    let index = SelectionExpr::try_from(selstr)?.apply_whole(&topology, &state)?;
    if index.len() > 0 {
        Ok(index)
    } else {
        Err(SelectionError::FromExpr {
            expr_str: selstr.into(), 
            source: SelectionIndexError::IndexEmpty, 
        })
    }
}

pub(super) fn index_from_range(range: Range<usize>, n: usize) -> Result<SortedSet<usize>, SelectionError> {
    if range.start > n || range.end > n {
        Err(SelectionError::FromRange {
            first: range.start,
            last: range.end, 
            source: SelectionIndexError::IndexOutOfBounds(range.start, range.end, n), 
        })
    } else if range.len() > 0 {
        unsafe { Ok(SortedSet::from_sorted(range.clone().collect())) }
    } else {
        Err(SelectionError::FromRange {
            first: range.start,
            last: range.end, 
            source: SelectionIndexError::IndexEmpty, 
        })
    }
}

pub(super) fn index_from_vec(vec: Vec<usize>, n: usize) -> Result<SortedSet<usize>, SelectionError> {
    let ind = SortedSet::from_unsorted(vec);
    if ind.is_empty() {        
        Err(SelectionError::FromVec {
            first: ind[0],
            last: ind[ind.len()-1],
            size: ind.len(), 
            source: SelectionIndexError::IndexEmpty, 
        })
    } else if ind[0] > n || ind[ind.len()-1] > n {
        Err(SelectionError::FromVec {
            first: ind[0],
            last: ind[ind.len()-1],
            size: ind.len(), 
            source: SelectionIndexError::IndexOutOfBounds(ind[0], ind[ind.len()-1], n), 
        })
    } else {
        Ok(ind)
    }
}

pub(super) fn index_from_iter(it: impl Iterator<Item = usize>, n: usize) -> Result<SortedSet<usize>, SelectionError> {
    index_from_vec(it.collect(), n)   
}

// Macro for implementing traits
macro_rules! impl_read_only_source_traits {
    ( $t:ty ) => {
        impl TopologyProvider for $t {
            fn num_atoms(&self) -> usize {
                self.state.num_coords()
            }
        }
        
        impl AtomsProvider for $t {
            fn iter_atoms(&self) -> impl AtomIterator<'_> {
                self.topology.iter_atoms()
            }
        }
        
        impl StateProvider for $t {
            fn get_time(&self) -> f32 {
                self.state.get_time()
            }
        
            fn num_coords(&self) -> usize {
                self.state.num_coords()
            }
        }
        
        impl PosProvider for $t {
            fn iter_pos(&self) -> impl PosIterator<'_> {
                self.state.iter_pos()
            }
        }
        
        impl BoxProvider for $t {
            fn get_box(&self) -> Option<&PeriodicBox> {
                self.state.get_box()
            }
        }
        
        impl WritableToFile for $t {}

        impl ParticleProvider for $t {
            fn iter_particle(&self) -> impl ExactSizeIterator<Item = Particle<'_>> {
                self.iter_index().map(|i| Particle {
                    id: i,
                    atom: unsafe{self.topology.nth_atom_unchecked(i)},
                    pos: unsafe{self.state.nth_pos_unchecked(i)},
                })
            }
        }

        impl LenProvider for $t {
            fn len(&self) -> usize {
                self.state.num_coords()
            }
        }
        
        impl RandomPos for $t {
            unsafe fn nth_pos_unchecked(&self, i: usize) -> &Pos {
                self.state.nth_pos_unchecked(i)
            }
        }
        
    };
}

macro_rules! impl_read_write_source_traits {
    ( $t:ty ) => {
        impl_read_only_source_traits!($t);

        impl PosMutProvider for $t {
            fn iter_pos_mut(&self) -> impl PosMutIterator<'_> {
                self.state.iter_pos_mut()
            }
        }

        impl RandomPosMut for $t {
            fn nth_pos_mut(&self, i: usize) -> Option<&mut Pos> {
                self.state.nth_pos_mut(i)
            }

            unsafe fn nth_pos_mut_unchecked(&self, i: usize) -> &mut Pos {
                self.state.nth_pos_mut_unchecked(i)
            }
        }

        impl AtomsMutProvider for $t {
            fn iter_atoms_mut(&self) -> impl AtomMutIterator<'_> {
                self.topology.iter_atoms_mut()
            }
        }
        
        impl ParticleMutProvider for $t {
            fn iter_particle_mut(&self) -> impl ExactSizeIterator<Item = ParticleMut<'_>> {
                self.iter_index().map(|i| ParticleMut {
                    id: i,
                    atom: unsafe{self.topology.nth_atom_unchecked_mut(i)},
                    pos: unsafe{self.state.nth_pos_unchecked_mut(i)},
                })
            }
        }
    }
}

pub(crate) use impl_read_only_source_traits;
pub(crate) use impl_read_write_source_traits;

//------------------------------
// Operations on sorted vectors
//------------------------------

pub fn union_sorted<T: Ord + Clone + Copy>(lhs: &SortedSet<T>, rhs: &SortedSet<T>) -> SortedSet<T> {
    let mut l = 0;
    let mut r = 0;

    let mut ret = vec![];

    while l<lhs.len() && r<rhs.len() {
        if lhs[l] < rhs[r] {
            ret.push(lhs[l]);
            l += 1;
        } else if lhs[l] > rhs[r]{
            ret.push(rhs[r]);
            r += 1;
        } else {
            ret.push(lhs[l]);
            l += 1;
            r += 1;
        }
    }

    if l<lhs.len() {
        ret.extend(lhs[l..].into_iter().cloned());
    } else if r<rhs.len() {
        ret.extend(rhs[r..].into_iter().cloned());
    }

    unsafe {SortedSet::from_sorted(ret)}
}

pub fn intersection_sorted<T: Ord + Clone + Copy>(lhs: &SortedSet<T>, rhs: &SortedSet<T>) -> SortedSet<T> {
    let mut l = 0;
    let mut r = 0;

    let mut ret = vec![];

    while l<lhs.len() && r<rhs.len() {
        if lhs[l] < rhs[r] {
            while l<lhs.len() && lhs[l] < rhs[r] {l+=1}
        } else if lhs[l] > rhs[r]{
            while r<lhs.len() && lhs[l] > rhs[r] {r+=1}
        } else {
            ret.push(lhs[l]);
            l += 1;
            r += 1;
        }
    }

    unsafe {SortedSet::from_sorted(ret)}
}

pub fn difference_sorted<T: Ord + Clone + Copy>(lhs: &SortedSet<T>, rhs: &SortedSet<T>) -> SortedSet<T> {
    let mut l = 0;
    let mut r = 0;

    let mut ret = vec![];

    while l<lhs.len() && r<rhs.len() {
        if lhs[l] < rhs[r] {
            while l<lhs.len() && lhs[l] < rhs[r] {
                ret.push(lhs[l]);
                l+=1;
            }
        } else if lhs[l] > rhs[r]{
            while r<lhs.len() && lhs[l] > rhs[r] {r+=1}
        } else {
            l += 1;
            r += 1;
        }
    }

    if l<lhs.len() {
        ret.extend(lhs[l..].into_iter().cloned());
    }

    unsafe {SortedSet::from_sorted(ret)}
}


#[cfg(test)]
mod tests {
    use super::*;
    use sorted_vec::SortedSet;

    #[test]
    fn test_union_sorted() {
        let set1: SortedSet<usize> = SortedSet::from_unsorted(vec![1, 3, 5, 7, 10, 12]);
        let set2: SortedSet<usize> = SortedSet::from_unsorted(vec![2, 3, 6, 8]);
        let result = union_sorted(&set1, &set2);
        let expected: SortedSet<usize> = SortedSet::from_unsorted(vec![1, 2, 3, 5, 6, 7, 8, 10, 12]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_union_sorted_empty() {
        let set1: SortedSet<usize> = SortedSet::from_unsorted(vec![]);
        let set2: SortedSet<usize> = SortedSet::from_unsorted(vec![]);
        let result = union_sorted(&set1, &set2);
        let expected: SortedSet<usize> = SortedSet::from_unsorted(vec![]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_union_sorted_disjoint() {
        let set1: SortedSet<usize> = SortedSet::from_unsorted(vec![2, 3]);
        let set2: SortedSet<usize> = SortedSet::from_unsorted(vec![4, 5, 6]);
        let result = union_sorted(&set1, &set2);
        let expected: SortedSet<usize> = SortedSet::from_unsorted(vec![2, 3, 4, 5, 6]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_union_sorted_subset() {
        let set1: SortedSet<usize> = SortedSet::from_unsorted(vec![1, 2, 3]);
        let set2: SortedSet<usize> = SortedSet::from_unsorted(vec![2, 3]);
        let result = union_sorted(&set1, &set2);
        let expected: SortedSet<usize> = SortedSet::from_unsorted(vec![1, 2, 3]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_intersection_sorted() {
        let set1: SortedSet<usize> = SortedSet::from_unsorted(vec![1, 3, 5, 7, 10, 12]);
        let set2: SortedSet<usize> = SortedSet::from_unsorted(vec![2, 3, 6, 7, 10]);
        let result = intersection_sorted(&set1, &set2);
        let expected: SortedSet<usize> = SortedSet::from_unsorted(vec![3, 7, 10]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_intersection_sorted_empty() {
        let set1: SortedSet<usize> = SortedSet::from_unsorted(vec![1, 2, 3]);
        let set2: SortedSet<usize> = SortedSet::from_unsorted(vec![4, 5, 6]);
        let result = intersection_sorted(&set1, &set2);
        let expected: SortedSet<usize> = SortedSet::from_unsorted(vec![]);
        println!("{:?}",result);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_intersection_sorted_subset() {
        let set1: SortedSet<usize> = SortedSet::from_unsorted(vec![1, 2, 3, 4, 5]);
        let set2: SortedSet<usize> = SortedSet::from_unsorted(vec![2, 3]);
        let result = intersection_sorted(&set1, &set2);
        let expected: SortedSet<usize> = SortedSet::from_unsorted(vec![2, 3]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_intersection_sorted_identical() {
        let set1: SortedSet<usize> = SortedSet::from_unsorted(vec![1, 2, 3]);
        let set2: SortedSet<usize> = SortedSet::from_unsorted(vec![1, 2, 3]);
        let result = intersection_sorted(&set1, &set2);
        let expected: SortedSet<usize> = SortedSet::from_unsorted(vec![1, 2, 3]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_difference_sorted() {
        let set1: SortedSet<usize> = SortedSet::from_unsorted(vec![3, 5, 7, 10, 12, 13]);
        let set2: SortedSet<usize> = SortedSet::from_unsorted(vec![2, 3, 6, 7, 10]);
        let result = difference_sorted(&set1, &set2);
        let expected: SortedSet<usize> = SortedSet::from_unsorted(vec![5, 12, 13]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_difference_sorted_empty() {
        let set1: SortedSet<usize> = SortedSet::from_unsorted(vec![1, 2, 3]);
        let set2: SortedSet<usize> = SortedSet::from_unsorted(vec![]);
        let result = difference_sorted(&set1, &set2);
        let expected: SortedSet<usize> = SortedSet::from_unsorted(vec![1, 2, 3]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_difference_sorted_subset() {
        let set1: SortedSet<usize> = SortedSet::from_unsorted(vec![1, 2, 3, 4, 5]);
        let set2: SortedSet<usize> = SortedSet::from_unsorted(vec![2, 3]);
        let result = difference_sorted(&set1, &set2);
        let expected: SortedSet<usize> = SortedSet::from_unsorted(vec![1, 4, 5]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_difference_sorted_identical() {
        let set1: SortedSet<usize> = SortedSet::from_unsorted(vec![1, 2, 3]);
        let set2: SortedSet<usize> = SortedSet::from_unsorted(vec![1, 2, 3, 5]);
        let result = difference_sorted(&set1, &set2);
        let expected: SortedSet<usize> = SortedSet::from_unsorted(vec![]);
        assert_eq!(result, expected);
    }
}

