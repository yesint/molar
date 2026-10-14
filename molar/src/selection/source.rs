#![doc = include_str!("../../../docs/selection-sources.md")]

use crate::prelude::*;
use std::{
    cell::RefCell,
    collections::HashMap,
    sync::atomic::{AtomicU64, Ordering},
};
use thiserror::Error;

/// Identity of a prepared-reference table. Keep it when creating a binding adapter.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct SelectionResolverId(u64);
impl SelectionResolverId {
    pub fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(1);
        Self(NEXT.fetch_add(1, Ordering::Relaxed))
    }
    /// Make a key for one entry in this resolver's table.
    pub fn reference(self, slot: usize) -> SelectionReferenceId {
        SelectionReferenceId {
            resolver: self,
            slot,
        }
    }
}
impl Default for SelectionResolverId {
    fn default() -> Self {
        Self::new()
    }
}

/// Opaque prepared-reference key. Only the caller assigns meaning to `slot`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct SelectionReferenceId {
    resolver: SelectionResolverId,
    slot: usize,
}
impl SelectionReferenceId {
    pub fn slot(self) -> usize {
        self.slot
    }
    pub fn resolver(self) -> SelectionResolverId {
        self.resolver
    }
}

/// Caller-defined preparation or binding failure.
#[derive(Debug, Error)]
#[error("{message}")]
pub struct SelectionResolveError {
    pub message: String,
}
impl SelectionResolveError {
    pub fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
        }
    }
}

/// Resolve names in two stages. Implementors own the prepared-reference table.
///
/// `prepare` is called once for each distinct variable string during construction.
/// `bind` is called at most once per reference during an evaluation. Each subsequent
/// evaluation binds again. Adapters over different scenes or frames must retain the
/// identity and table of the resolver used to construct the expression.
pub trait SelectionResolver {
    fn identity(&self) -> SelectionResolverId;
    fn prepare(&mut self, variable: &str) -> Result<SelectionReferenceId, SelectionResolveError>;
    fn bind(
        &self,
        reference: SelectionReferenceId,
    ) -> Result<SelectionSourceView<'_>, SelectionResolveError>;
}

/// Adapter from a selection-like object to a concrete borrowed source view.
pub trait SelectionSource {
    fn selection_view(&self) -> SelectionSourceView<'_>;
}

/// Checked borrowed source data. Indices address the full atom/coordinate arrays.
/// Empty domains are allowed for reference operands. Coordinates use nanometers.
/// Source identity is the shared atom storage plus coordinate buffer; different
/// selections over that same data have the same identity. All compared sources
/// must use the same coordinate frame.
#[derive(Clone, Copy, Debug)]
pub struct SelectionSourceView<'a> {
    atoms: &'a AtomStorage,
    coords: &'a [Pos],
    indices: Option<&'a [usize]>,
    bonds: Option<&'a BondStorage>,
    pbox: Option<&'a PeriodicBox>,
    velocities: &'a [Vel],
    forces: &'a [Force],
}
impl<'a> SelectionSourceView<'a> {
    pub fn from_buffers(
        atoms: &'a AtomStorage,
        coords: &'a [Pos],
        indices: Option<&'a [usize]>,
    ) -> Result<Self, SelectionError> {
        if atoms.len() != coords.len() {
            return Err(TopologyStateSizesError(atoms.len(), coords.len()).into());
        }
        if let Some(ind) = indices {
            if ind.windows(2).any(|w| w[0] >= w[1]) || ind.last().is_some_and(|&i| i >= atoms.len())
            {
                return Err(SelectionError::InvalidSourceIndices);
            }
        }
        Ok(Self {
            atoms,
            coords,
            indices,
            bonds: None,
            pbox: None,
            velocities: &[],
            forces: &[],
        })
    }
    pub fn from_parts(
        topology: &'a Topology,
        state: &'a State,
        indices: Option<&'a [usize]>,
    ) -> Result<Self, SelectionError> {
        Self::from_buffers(&topology.atoms, &state.coords, indices)?
            .with_bonds(&topology.bonds)?
            .with_box(state.get_box())
            .with_velocities(&state.velocities)?
            .with_forces(&state.forces)
    }
    pub fn with_box(mut self, pbox: Option<&'a PeriodicBox>) -> Self {
        self.pbox = pbox;
        self
    }
    pub fn with_bonds(mut self, bonds: &'a BondStorage) -> Result<Self, SelectionError> {
        if bonds
            .iter()
            .any(|b| b.pair().iter().any(|&i| i >= self.atoms.len()))
        {
            return Err(SelectionError::InvalidSourceIndices);
        }
        self.bonds = Some(bonds);
        Ok(self)
    }
    pub fn with_velocities(mut self, velocities: &'a [Vel]) -> Result<Self, SelectionError> {
        if !velocities.is_empty() && velocities.len() != self.coords.len() {
            return Err(SelectionError::InvalidSourceIndices);
        }
        self.velocities = velocities;
        Ok(self)
    }
    pub fn with_forces(mut self, forces: &'a [Force]) -> Result<Self, SelectionError> {
        if !forces.is_empty() && forces.len() != self.coords.len() {
            return Err(SelectionError::InvalidSourceIndices);
        }
        self.forces = forces;
        Ok(self)
    }
    pub fn indices(&self) -> Option<&'a [usize]> {
        self.indices
    }
    pub(super) fn same_source<S: AtomProvider + PosProvider>(&self, other: &S) -> bool {
        std::ptr::eq(self.atoms, other.atom_storage())
            && self.coords.as_ptr() == unsafe { other.coords_ptr() }
    }
}
impl SelectionSource for SelectionSourceView<'_> {
    fn selection_view(&self) -> SelectionSourceView<'_> {
        *self
    }
}
impl LenProvider for SelectionSourceView<'_> {
    fn len(&self) -> usize {
        self.indices.map_or(self.coords.len(), |i| i.len())
    }
}
impl IndexProvider for SelectionSourceView<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        self.indices
            .map_or(i, |ind| unsafe { *ind.get_unchecked(i) })
    }
}
impl AtomProvider for SelectionSourceView<'_> {
    fn atom_storage(&self) -> &AtomStorage {
        self.atoms
    }
}
impl PosProvider for SelectionSourceView<'_> {
    unsafe fn coords_ptr(&self) -> *const Pos {
        self.coords.as_ptr()
    }
}
impl BoxProvider for SelectionSourceView<'_> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.pbox
    }
}
impl VelProvider for SelectionSourceView<'_> {
    unsafe fn vel_ptr(&self) -> *const Vel {
        if self.velocities.is_empty() {
            std::ptr::null()
        } else {
            self.velocities.as_ptr()
        }
    }
}
impl ForceProvider for SelectionSourceView<'_> {
    unsafe fn force_ptr(&self) -> *const Force {
        if self.forces.is_empty() {
            std::ptr::null()
        } else {
            self.forces.as_ptr()
        }
    }
}
impl BondProvider for SelectionSourceView<'_> {
    fn num_bonds(&self) -> usize {
        self.bonds.map_or(0, |b| b.len())
    }
    unsafe fn get_bond_unchecked(&self, i: usize) -> BondRef<'_> {
        unsafe { self.bonds.unwrap().get_unchecked(i) }
    }
    fn iter_bonds(&self) -> impl Iterator<Item = BondRef<'_>> {
        self.bonds.into_iter().flat_map(|b| b.iter())
    }
}

impl SelectionSource for System {
    fn selection_view(&self) -> SelectionSourceView<'_> {
        // System and bound selections already validate these arrays and indices.
        SelectionSourceView {
            atoms: &self.top.atoms,
            coords: &self.st.coords,
            indices: None,
            bonds: Some(&self.top.bonds),
            pbox: self.st.get_box(),
            velocities: &self.st.velocities,
            forces: &self.st.forces,
        }
    }
}
impl<T: SystemProvider + IndexSliceProvider> SelectionSource for T {
    fn selection_view(&self) -> SelectionSourceView<'_> {
        let mut view = unsafe { &*self.get_system_ptr() }.selection_view();
        view.indices = Some(self.get_index_slice());
        view
    }
}

/// A cache lives only as long as one evaluation and borrows the binding resolver.
pub(super) struct EvaluationSources<'a> {
    resolver: &'a dyn SelectionResolver,
    cache: RefCell<HashMap<SelectionReferenceId, SelectionSourceView<'a>>>,
}
impl<'a> EvaluationSources<'a> {
    pub fn new(resolver: &'a dyn SelectionResolver) -> Self {
        Self {
            resolver,
            cache: RefCell::new(HashMap::new()),
        }
    }
    pub fn bind(
        &self,
        id: SelectionReferenceId,
        name: &str,
        span: &std::ops::Range<usize>,
    ) -> Result<SelectionSourceView<'a>, SelectionParserError> {
        if let Some(view) = self.cache.borrow().get(&id) {
            return Ok(*view);
        }
        let view = self
            .resolver
            .bind(id)
            .map_err(|error| SelectionParserError::Variable {
                name: name.into(),
                span: span.clone(),
                error,
            })?;
        self.cache.borrow_mut().insert(id, view);
        Ok(view)
    }
}

// Erase the cache's lifetime so evaluation contexts can borrow shorter subsets.
pub(super) trait SourceBindings {
    fn bind<'a>(
        &'a self,
        id: SelectionReferenceId,
        name: &str,
        span: &std::ops::Range<usize>,
    ) -> Result<SelectionSourceView<'a>, SelectionParserError>;
}
impl SourceBindings for EvaluationSources<'_> {
    fn bind<'a>(
        &'a self,
        id: SelectionReferenceId,
        name: &str,
        span: &std::ops::Range<usize>,
    ) -> Result<SelectionSourceView<'a>, SelectionParserError> {
        EvaluationSources::bind(self, id, name, span)
    }
}
