use crate::prelude::*;
use sorted_vec::SortedSet;
use super::utils::local_to_global;

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
        match subset {
            None => Ok(self.apply_whole(top, st)?),
            Some(sub) => {
                let ind = self.apply_subset(top, st, sub)?;
                if ind.is_empty() {
                    Err(SelectionError::Empty)
                } else {
                    Ok(ind)
                }
            }
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

        if self.start >= n || self.end >= n {
            return Err(SelectionError::IndexCheck(self.start, self.end, n));
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
            Err(SelectionError::Empty)
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
            Err(SelectionError::Empty)
        } else {            
            match subset {
                None => {
                    let v: SortedSet<usize> = self.to_vec().into();
                    let n = top.num_atoms();
                    if v[0] >= n || v[v.len()-1] >= n {
                        Err(SelectionError::IndexCheck(v[0], v[v.len()-1], n))
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
            Err(SelectionError::Empty)
        } else {            
            match subset {
                None => {                    
                    let n = top.num_atoms();
                    if self[0] >= n || self[self.len()-1] >= n {
                        Err(SelectionError::IndexCheck(self[0], self[self.len()-1], n))
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
