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

pub(super) fn index_from_expr(expr: &SelectionExpr, topology: &Topology, state: &State) -> Result<SortedSet<usize>, SelectionError> {
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

pub(super) fn index_from_range(range: &Range<usize>, n: usize) -> Result<SortedSet<usize>, SelectionError> {
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

pub(super) fn index_from_vec(vec: &Vec<usize>, n: usize) -> Result<SortedSet<usize>, SelectionError> {
    let ind = SortedSet::from_unsorted(vec.clone());
    if ind.is_empty() {        
        Err(SelectionError::FromVec {
            first: vec[0],
            last: vec[vec.len()-1],
            size: vec.len(), 
            source: SelectionIndexError::IndexEmpty, 
        })
    } else if ind[0] > n || ind[ind.len()-1] > n {
        Err(SelectionError::FromVec {
            first: vec[0],
            last: vec[vec.len()-1],
            size: vec.len(), 
            source: SelectionIndexError::IndexOutOfBounds(vec[0], vec[vec.len()-1], n), 
        })
    } else {
        Ok(ind)
    }
}

pub(super) fn index_from_iter(it: impl Iterator<Item = usize>, n: usize) -> Result<SortedSet<usize>, SelectionError> {
    index_from_vec(&it.collect(), n)   
}
