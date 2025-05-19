use std::{collections::HashMap, path::Path};

use crate::prelude::*;
use itertools::Itertools;
use sorted_vec::SortedSet;
use super::utils::local_to_global;

/// Trait for selection definitions. 
/// Implementors could be used as argumnets for selection creation methods.
pub trait SelectionDef {
    /// All errors, including bounds checks, have to be captured, 
    /// caller assumes that the index is correct.
    fn into_sel_index(
        self,
        top: &Topology,
        st: &State,
        subset: Option<&[usize]>,
    ) -> Result<SortedSet<usize>, SelectionError>;
}

impl SelectionDef for &SelectionExpr {
    fn into_sel_index(
        self,
        top: &Topology,
        st: &State,
        subset: Option<&[usize]>,
    ) -> Result<SortedSet<usize>, SelectionError> {        
        let ind = match subset {
            None => self.apply_whole(top, st)?,
            Some(sub) => self.apply_subset(top, st, sub)?,
        };
        if ind.is_empty() {
            Err(SelectionError::EmptyExpr(self.get_str().to_string()))
        } else {
            Ok(ind)
        }
    }
}

impl SelectionDef for &str {
    fn into_sel_index(
        self,
        top: &Topology,
        st: &State,
        subset: Option<&[usize]>,
    ) -> Result<SortedSet<usize>, SelectionError> {
        SelectionExpr::new(self)?.into_sel_index(top, st, subset)
    }
}

impl SelectionDef for String {
    fn into_sel_index(
        self,
        top: &Topology,
        st: &State,
        subset: Option<&[usize]>,
    ) -> Result<SortedSet<usize>, SelectionError> {
        self.as_str().into_sel_index(top, st, subset)
    }
}

impl SelectionDef for &String {
    fn into_sel_index(
        self,
        top: &Topology,
        st: &State,
        subset: Option<&[usize]>,
    ) -> Result<SortedSet<usize>, SelectionError> {
        self.as_str().into_sel_index(top, st, subset)
    }
}

impl SelectionDef for std::ops::Range<usize> {
    fn into_sel_index(
        self,
        top: &Topology,
        _st: &State,
        subset: Option<&[usize]>,
    ) -> Result<SortedSet<usize>, SelectionError> {
        let n = match subset {
            None => top.num_atoms(),
            Some(sub) => sub.len(),
        };

        if self.start >= n || self.end > n {
            return Err(SelectionError::IndexValidation(self.start, self.end, n-1));
        }

        match subset {
            None => Ok(unsafe { SortedSet::from_sorted(self.collect::<Vec<_>>()) }),
            Some(sub) => {                
                Ok(self.map(|i| sub[i]).collect::<Vec<_>>().into())
            }
        }
    }
}

impl SelectionDef for std::ops::RangeInclusive<usize> {
    fn into_sel_index(
        self,
        top: &Topology,
        st: &State,
        subset: Option<&[usize]>,
    ) -> Result<SortedSet<usize>, SelectionError> {
        if self.is_empty() {
            Err(SelectionError::EmptyRange(*self.start(),*self.end()))
        } else {
            let (b, e) = self.into_inner();
            (b..e + 1).into_sel_index(top, st, subset)
        }
    }
}

impl SelectionDef for &[usize] {
    fn into_sel_index(
        self,
        top: &Topology,
        _st: &State,
        subset: Option<&[usize]>,
    ) -> Result<SortedSet<usize>, SelectionError> {
        if self.is_empty() {
            Err(SelectionError::EmptySlice)
        } else {            
            match subset {
                None => {
                    let v: SortedSet<usize> = self.to_vec().into();
                    let n = top.num_atoms();
                    if v[0] >= n || v[v.len()-1] >= n {
                        Err(SelectionError::IndexValidation(v[0], v[v.len()-1], n-1))
                    } else {
                        Ok(v)
                    }
                },
                Some(sub) => local_to_global(self.into_iter().cloned(), sub),
            }            
        }
    }
}

impl SelectionDef for &Vec<usize> {
    fn into_sel_index(
        self,
        top: &Topology,
        st: &State,
        subset: Option<&[usize]>,
    ) -> Result<SortedSet<usize>, SelectionError> {
        self.as_slice().into_sel_index(top, st, subset)
    }
}

impl SelectionDef for Vec<usize> {
    fn into_sel_index(
        self,
        top: &Topology,
        st: &State,
        subset: Option<&[usize]>,
    ) -> Result<SortedSet<usize>, SelectionError> {
        self.as_slice().into_sel_index(top, st, subset)
    }
}

impl SelectionDef for SVec {
    fn into_sel_index(
        self,
        top: &Topology,
        _st: &State,
        subset: Option<&[usize]>,
    ) -> Result<SortedSet<usize>, SelectionError> {
        if self.is_empty() {
            Err(SelectionError::EmptySlice)
        } else {            
            match subset {
                None => {                    
                    let n = top.num_atoms();
                    if self[0] >= n || self[self.len()-1] >= n {
                        Err(SelectionError::IndexValidation(self[0], self[self.len()-1], n-1))
                    } else {
                        Ok(self)
                    }
                },
                Some(sub) => local_to_global(self.iter().cloned(), sub),
            }            
        }
    }
}

impl SelectionDef for &SVec {
    fn into_sel_index(
        self,
        top: &Topology,
        st: &State,
        subset: Option<&[usize]>,
    ) -> Result<SortedSet<usize>, SelectionError> {
        self.clone().into_sel_index(top, st, subset)
    }
}

impl<K: UserCreatableKind> SelectionDef for &Sel<K> {
    fn into_sel_index(
        self,
        top: &Topology,
        _st: &State,
        subset: Option<&[usize]>,
    ) -> Result<SortedSet<usize>, SelectionError> {
        if let Some(_) = subset {
            return Err(SelectionError::SelDefInSubsel)
        }

        let n = top.num_atoms();
        if self.first_index() >= n || self.last_index() >= n {
            return Err(SelectionError::IndexValidation(self.first_index(), self.last_index(), n-1));
        }

        Ok(self.get_index_vec().clone())
    }
}

/// Representation of Gromacs index files
pub struct NdxFile {
    groups: HashMap<String, SortedSet<usize>>,
}

impl NdxFile {
    /// Creates a new NdxFile by parsing a Gromacs index file
    pub fn new(path: impl AsRef<Path>) -> Result<Self, SelectionError> {
        let path = path.as_ref();
        let ndx_str = std::fs::read_to_string(path)
            .map_err(|e| NdxError::NdxIo(path.to_owned(), e))?;
        
        let mut groups = HashMap::new();
        let mut current_group = None;
        let mut current_numbers = Vec::new();

        for line in ndx_str.lines() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }

            if line.starts_with('[') && line.ends_with(']') {
                // Store previous group if exists
                if let Some(group_name) = current_group.take() {
                    if !current_numbers.is_empty() {
                        groups.insert(group_name, current_numbers.into());
                    }
                }
                // Start new group
                let group_name = line[1..line.len()-1].trim().to_string();
                current_group = Some(group_name);
                current_numbers = Vec::new();
            } else if let Some(group_name) = &current_group {
                // Parse numbers for current group
                current_numbers.extend(
                    line.split_whitespace()
                        .map(|s| s.parse::<usize>())
                        .map_ok(|i| i - 1) // Convert to zero-based
                        .collect::<Result<Vec<_>, _>>()
                        .map_err(|e| NdxError::Parse(group_name.clone(), e))?,
                );
            } else {
                return Err(NdxError::MalformedNdxFile(path.to_owned()))?;
            }
        }

        // Store last group if exists
        if let Some(group_name) = current_group {
            if !current_numbers.is_empty() {
                groups.insert(group_name, current_numbers.into());
            }
        }

        Ok(Self { groups })
    }

    /// Get an index group by name
    pub fn get_group(&self, name: impl AsRef<str>) -> Result<&SortedSet<usize>, SelectionError> {
        let gr = name.as_ref();
        Ok(self.groups.get(gr).ok_or_else(|| NdxError::NoGroup(gr.to_owned()))?)
    }
}
