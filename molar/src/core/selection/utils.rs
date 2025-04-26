use sorted_vec::SortedSet;
use crate::prelude::*;

/// Verifies that topology and state have matching number of atoms
pub(crate) fn check_topology_state_sizes(topology: &Topology, state: &State) -> Result<(),super::TopologyStateSizes> {
    let n1 = topology.num_atoms();
    let n2 = state.num_pos();
    if n1 != n2 { Err(super::TopologyStateSizes(n1,n2))? }
    Ok(())
}

//------------------------------
// Operations on sorted vectors
//------------------------------

/// Computes the union of two sorted sets
/// 
/// Returns a new set containing all elements that are in either set
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

/// Computes the intersection of two sorted sets
/// 
/// Returns a new set containing elements that appear in both sets
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

/// Computes the difference between two sorted sets
/// 
/// Returns a new set containing elements from the first set that do not appear in the second set
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

// Convert local selection indexes to global
pub(super) fn local_to_global(
    local: impl Iterator<Item = usize>,
    subset: &[usize],
) -> Result<SortedSet<usize>, SelectionError> {
    Ok(local
        .map(|i| {
            subset
                .get(i)
                .cloned()
                .ok_or_else(|| SelectionError::LocalToGlobal(i, subset.len()-1))
        })
        .collect::<Result<Vec<_>, _>>()?
        .into())
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

