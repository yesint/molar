//! The owned bond row.
//!
//! [`Bond`] is the detached construction/interchange type for a chemical bond, mirroring the
//! role [`Atom`](crate::Atom) plays for [`AtomStorage`](crate::AtomStorage): IO readers and
//! builders assemble `Bond`s and [`BondStorage::push`](crate::BondStorage::push) scatters each
//! one into the columns. Bonds already in storage are read back through the
//! [`BondRef`](crate::BondRef) proxy — there is no `&Bond` to borrow.

use serde::{Deserialize, Serialize};

/// Chemical bond order. File formats that don't record it (PDB, GRO, XYZ, …) yield
/// [`BondOrder::Unspecified`]; formats that do (e.g. SDF) set the concrete order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum BondOrder {
    /// The source didn't record an order (the common case for guessed/PDB bonds).
    #[default]
    Unspecified,
    Single,
    Double,
    Triple,
    Aromatic,
}

/// A bond between two atoms (by global index) with an optional chemical order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Bond {
    pub i1: usize,
    pub i2: usize,
    pub order: BondOrder,
}

impl Bond {
    /// A bond between `i1` and `i2` with an unspecified order.
    pub fn new(i1: usize, i2: usize) -> Self {
        Self { i1, i2, order: BondOrder::Unspecified }
    }

    /// A bond with an explicit order.
    pub fn with_order(i1: usize, i2: usize, order: BondOrder) -> Self {
        Self { i1, i2, order }
    }

    /// The two atom indices as a pair (for code that just needs the endpoints).
    pub fn pair(&self) -> [usize; 2] {
        [self.i1, self.i2]
    }

    /// Whether `idx` is one of this bond's endpoints.
    pub fn contains(&self, idx: usize) -> bool {
        self.i1 == idx || self.i2 == idx
    }
}
