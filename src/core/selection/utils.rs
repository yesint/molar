use std::ops::Range;
use anyhow::{Result,bail};
use crate::prelude::*;
use itertools::Itertools;

pub(super) fn check_sizes(topology: &Topology, state: &State) -> Result<()> {
    let n1 = topology.num_atoms();
    let n2 = state.num_coords();
    match n1 == n2 {
        true => Ok(()),
        false => bail!("Structure and State are incompatible (sizes {n1} and {n2})"),
    }
}

pub(super) fn index_from_all(n: usize) -> Vec<usize> {
    (0..n).collect()
}

pub(super) fn index_from_expr(expr: &SelectionExpr, topology: &Topology, state: &State) -> Result<Vec<usize>> {
    let index = expr.apply_whole(&topology, &state)?;
    if index.len() > 0 {
        Ok(index)
    } else {
        bail!(format!("Selection constructed from expr '{}' is empty",expr.get_str()))
    }
}

pub(super) fn index_from_str(selstr: &str, topology: &Topology, state: &State) -> Result<Vec<usize>> {
    let index = SelectionExpr::try_from(selstr)?.apply_whole(&topology, &state)?;
    if index.len() > 0 {
        Ok(index)
    } else {
        bail!(format!("Selection constructed from string '{selstr}' is empty"))
    }
}

pub(super) fn index_from_range(range: &Range<usize>, n: usize) -> Result<Vec<usize>> {
    if range.start > n - 1 || range.end > n - 1 {
        bail!(
            "Range {}:{} is invalid, 0:{} is allowed for constructing selection",
            range.start,
            range.end,
            n
        );
    }
    if range.len() > 0 {
        Ok(range.clone().collect())
    } else {
        bail!(format!("Selection constructed from range {}:{} is empty",range.start,range.end))
    }
}

pub(super) fn index_from_iter(it: impl Iterator<Item = usize>, n: usize) -> Result<Vec<usize>> {
    let index: Vec<usize> = it.sorted().dedup().collect();
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
