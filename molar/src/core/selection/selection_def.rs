use crate::prelude::*;
use sorted_vec::SortedSet;
use super::utils::local_to_global;

/// Trait for selection definitions. 
/// Implementors could be used as argumnets for selection creation methods.
pub trait SelectionDef {
    /// All errors, including bounds checks, have to be captured, 
    /// caller assumes that returned index is correct and does no additional checks..
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
            None => top.len(),
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
                    let n = top.len();
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
                    let n = top.len();
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

impl SelectionDef for &SubSystem<'_> {
    fn into_sel_index(
        self,
        top: &Topology,
        _st: &State,
        subset: Option<&[usize]>,
    ) -> Result<SortedSet<usize>, SelectionError> {
        if let Some(_) = subset {
            return Err(SelectionError::SelDefInSubsel)
        }

        let n = top.len();
        let first_index = unsafe{self.get_index_unchecked(0)};
        let last_index = unsafe{self.get_index_unchecked(self.len()-1)};
        if first_index >= n || last_index >= n {
            return Err(SelectionError::IndexValidation(first_index, last_index, n-1));
        }

        Ok( unsafe {SVec::from_sorted( self.iter_index().collect() )} )
    }
}

// impl SelectionDef for &SelSerial {
//     fn into_sel_index(
//         self,
//         top: &Topology,
//         _st: &State,
//         subset: Option<&[usize]>,
//     ) -> Result<SortedSet<usize>, SelectionError> {
//         if let Some(_) = subset {
//             return Err(SelectionError::SelDefInSubsel)
//         }

//         Ok(self.)
//     }
// }
