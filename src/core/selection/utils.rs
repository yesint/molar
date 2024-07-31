use std::ops::Range;
use anyhow::{Result,bail};
use sorted_vec::SortedSet;
use crate::prelude::*;

pub(super) fn check_sizes(topology: &Topology, state: &State) -> Result<()> {
    let n1 = topology.num_atoms();
    let n2 = state.num_coords();
    match n1 == n2 {
        true => Ok(()),
        false => bail!("Structure and State are incompatible (sizes {n1} and {n2})"),
    }
}

pub(super) fn index_from_all(n: usize) -> SortedSet<usize> {
    unsafe { SortedSet::from_sorted((0..n).collect()) }
}

pub(super) fn index_from_expr(expr: &SelectionExpr, topology: &Topology, state: &State) -> Result<SortedSet<usize>> {
    let index = expr.apply_whole(&topology, &state)?;
    if index.len() > 0 {
        Ok(index)
    } else {
        bail!(format!("Selection constructed from expr '{}' is empty",expr.get_str()))
    }
}

pub(super) fn index_from_str(selstr: &str, topology: &Topology, state: &State) -> Result<SortedSet<usize>> {
    let index = SelectionExpr::try_from(selstr)?.apply_whole(&topology, &state)?;
    if index.len() > 0 {
        Ok(index)
    } else {
        bail!(format!("Selection constructed from string '{selstr}' is empty"))
    }
}

pub(super) fn index_from_range(range: &Range<usize>, n: usize) -> Result<SortedSet<usize>> {
    if range.start > n || range.end > n {
        bail!(
            "Range {}:{} is invalid, 0:{} is allowed for constructing selection",
            range.start,
            range.end,
            n
        );
    }
    if range.len() > 0 {
        unsafe { Ok(SortedSet::from_sorted(range.clone().collect())) }
    } else {
        bail!(format!("Selection constructed from range {}:{} is empty",range.start,range.end))
    }
}

pub(super) fn index_from_vec(vec: &Vec<usize>, n: usize) -> Result<SortedSet<usize>> {
    let ind = SortedSet::from_unsorted(vec.clone());
    if ind.is_empty() {        
        bail!(format!("Selection constructed from vector is empty"))
    }
    if ind[0] > n || ind[ind.len()-1] > n {
        bail!(
            "Vector with indexes {}:{} is invalid, 0:{} is allowed for constructing selection",
            ind[0],ind[ind.len()-1],n
        );
    }
    Ok(ind)
}


pub(super) fn index_from_iter(it: impl Iterator<Item = usize>, n: usize) -> Result<SortedSet<usize>> {
    let index = SortedSet::from_unsorted(it.collect());
    if index.is_empty() {
        bail!("Iterator is empty, which results in empty selection")
    }
    if index[0] > n - 1 || index[index.len() - 1] > n - 1 {
        bail!(
            "Selection iterator range {}:{} is invalid, 0:{} is allowed",
            index[0],
            index[index.len() - 1],
            n
        );
    }
    Ok(index)
}
