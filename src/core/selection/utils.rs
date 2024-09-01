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

pub(super) fn index_from_expr_sub(expr: &SelectionExpr, topology: &Topology, state: &State, subset: &SortedSet<usize>) -> Result<SortedSet<usize>, SelectionError> {
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

        impl RandomPosMutProvider for $t {
            unsafe fn nth_pos_unchecked_mut(&self, i: usize) -> &mut Pos {
                self.state.nth_pos_unchecked_mut(i)
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