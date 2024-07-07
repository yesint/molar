use std::marker::PhantomData;
use crate::prelude::*;
use anyhow::{bail, Result};
use super::utils::*;

/// Source of serial selections.
/// 
/// [Source] produces mutable selections, which may overlap and are can only be
/// used from the same thread where they are created.
/// 
/// [Source] takes ownership of [Topology] and [State] so that after creating a [Source] they are no longer accessible outside it.
/// This guarantees correct access without data integrity issues.

#[derive(Default)]
pub struct Source {
    topology: TopologyArc,
    state: StateArc,
    _marker: PhantomData<*const ()>,
}

impl Source {
    /// Create the [Source] producing mutable selections that may overlap accessible from a single thread.
    pub fn new(
        topology: TopologyUArc,
        state: StateUArc,
    ) -> Result<Self> {
        check_sizes(&topology, &state)?;
        Ok(Source {
            topology: topology.shareable(),
            state: state.shareable(),
            _marker: Default::default(),
        })
    }

    pub(super) unsafe fn new_from_arc(
        topology: TopologyArc,
        state: StateArc,
    ) -> Self {        
        Source {
            topology,
            state,
            _marker: Default::default(),
        }
    }

    /// Release and return [Topology] and [State]. Fails if any selections created from this [Source] are still alive.
    pub fn release(self) -> anyhow::Result<(TopologyUArc, StateUArc)> {
        Ok((
            triomphe::Arc::try_unique(self.topology)
                .or_else(|_| bail!("Can't release topology: multiple references are active!"))?,
            triomphe::Arc::try_unique(self.state)
                .or_else(|_| bail!("Can't release state: multiple references are active!"))?,
        ))
    }

    //---------------------------------
    // Creating selections
    //---------------------------------

    
    /// Get a pointer to contained [Topology]
    pub(super) fn get_topology(&self) -> &TopologyArc {
        &self.topology
    }
    
    /// Get a pointer to contained [State]
    pub(super) fn get_state(&self) -> &StateArc {
        &self.state
    }
    
    /// Creates new selection from an iterator of indexes. Indexes are bound checked, sorted and duplicates are removed.
    /// If any index is out of bounds the error is returned.
    pub fn select_from_iter(
        &self,
        iter: impl Iterator<Item = usize>,
    ) -> anyhow::Result<Sel<MutableSerial>> {
        let vec = index_from_iter(iter, self.topology.num_atoms())?;
        Ok(Sel::new(
            triomphe::Arc::clone(&self.topology),
            triomphe::Arc::clone(&self.state),
            vec,            
        ))
    }

    /// Selects all
    pub fn select_all(&self) -> anyhow::Result<Sel<MutableSerial>> {
        let vec = index_from_all(self.topology.num_atoms());
        Ok(Sel::new(
            triomphe::Arc::clone(&self.topology),
            triomphe::Arc::clone(&self.state),
            vec,
        ))
    }

    /// Creates new selection from a selection expression string. Selection expression is constructed internally but
    /// can't be reused. Consider using [select_expr](Self::select_expr) if you already have selection expression.
    pub fn select_str(&self, selstr: &str) -> anyhow::Result<Sel<MutableSerial>> {
        let vec = index_from_str(selstr, &self.topology, &self.state)?;
        Ok(Sel::new(
            triomphe::Arc::clone(&self.topology),
            triomphe::Arc::clone(&self.state),
            vec,            
        ))
    }

    /// Creates new selection from an existing selection expression.
    pub fn select_expr(&self, expr: &SelectionExpr) -> anyhow::Result<Sel<MutableSerial>> {
        let vec = index_from_expr(expr, &self.topology, &self.state)?;
        Ok(Sel::new(
            triomphe::Arc::clone(&self.topology),
            triomphe::Arc::clone(&self.state),
            vec,
        ))
    }

    /// Creates new selection from a range of indexes.
    /// If rangeis out of bounds the error is returned.
    pub fn select_range(&self, range: &std::ops::Range<usize>) -> anyhow::Result<Sel<MutableSerial>> {
        let vec = index_from_range(range, self.topology.num_atoms())?;
        Ok(Sel::new(
            triomphe::Arc::clone(&self.topology),
            triomphe::Arc::clone(&self.state),
            vec,
        ))
    }

    /// Sets new [State] in this [Source]. All selections created from this [Source] will automatically view
    /// the new state. 
    /// 
    /// New state should be compatible with the old one (have the same number of atoms). If not, the error is returned.
    /// 
    /// Returns unique pointer to the old state, so it could be reused if needed.
    pub fn set_state(&mut self, state: StateUArc) -> Result<StateUArc> {
        // Check if the states are compatible
        if !self.state.interchangeable(&state) {
            bail!("Can't set state: states are incompatible!")
        }

        let ret = state.shareable();
        unsafe {
            std::ptr::swap(
                self.state.as_ptr() as *mut State,
                ret.as_ptr() as *mut State,
            )
        }
        triomphe::Arc::try_unique(ret).or_else(|_| bail!("Can't set state: multiple references are active!"))
    }

    /// Sets new topology in this [Source]. All selections created from this [Source] will automatically view
    /// the new [Topology]. 
    /// 
    /// New [Topology] should be compatible with the old one (have the same number of atoms, bonds, molecules, etc.). If not, the error is returned.
    /// 
    /// Returns unique pointer to the old topology, so it could be reused if needed.
    pub fn set_topology(&mut self, topology: TopologyUArc) -> Result<TopologyUArc> {
        // Check if the states are compatible
        if !self.topology.interchangeable(&topology) {
            bail!("Can't set topology: topologies are incompatible!")
        }

        let ret = topology.shareable();
        unsafe {
            std::ptr::swap(
                self.topology.as_ptr() as *mut State,
                ret.as_ptr() as *mut State,
            )
        }
        triomphe::Arc::try_unique(ret).or_else(|_| bail!("Can't set topology: multiple references are active!"))
    }
}