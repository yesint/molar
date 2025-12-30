use thiserror::Error;

use crate::{io::FileIoError, prelude::{BuilderError, NdxError, PeriodicBoxError}, selection::ast::SelectionParserError};
//############################################################
//#  Error enums
//############################################################

/// Error for different sizes of topology and state
#[derive(Error, Debug)]
#[error("topology and state have different sizes ({0},{1})")]
pub struct TopologyStateSizesError(pub(crate) usize, pub(crate) usize);

/// Error related to creation of selections
#[derive(Error, Debug)]
pub enum SelectionError {
    #[error("selection parser failed")]
    Parser(#[from] SelectionParserError),

    #[error(transparent)]
    DifferentSizes(#[from] TopologyStateSizesError),

    #[error("creating selection from expression {expr_str}")]
    FromExpr {
        expr_str: String,
        source: SelectionIndexError,
    },

    #[error("creating selection from range {first}:{last}")]
    FromRange {
        first: usize,
        last: usize,
        source: SelectionIndexError,
    },

    #[error("invalid local sub-range: {0}:{1}, valid range: 0:{2}")]
    LocalRange(usize, usize, usize),

    #[error("creating selection from vec {first}..{last} of size {size}")]
    FromVec {
        first: usize,
        last: usize,
        size: usize,
        source: SelectionIndexError,
    },

    #[error("index {0} is beyond the allowed range 0:{1}")]
    OutOfBounds(usize, usize),

    #[error("local index {0} is beyond the allowed range 0:{1}")]
    LocalToGlobal(usize, usize),

    #[error("selection index {0}:{1} is outside the source range: 0:{2}")]
    IndexValidation(usize, usize, usize),

    #[error(transparent)]
    FileIo(#[from] FileIoError),

    #[error(transparent)]
    Builder(#[from] BuilderError),

    #[error("can't set incompatible state")]
    IncompatibleState,

    #[error("can't set incompatible topology")]
    IncompatibleTopology,

    #[error("can't release source: multiple references are active")]
    Release,

    #[error("selection from vector slice is empty")]
    EmptySlice,

    #[error("selection range {0}:{1} is empty")]
    EmptyRange(usize, usize),

    #[error("selection '{0}' is empty")]
    EmptyExpr(String),

    #[error("splitting produced no selections")]
    EmptySplit,

    #[error("selection intersection is empy")]
    EmptyIntersection,

    #[error("selection difference is empy")]
    EmptyDifference,

    #[error("selection complement is empy")]
    EmptyComplement,

    #[error(transparent)]
    PeriodicBox(#[from] PeriodicBoxError),

    #[error("no molecules in topology")]
    NoMolecules,

    #[error("gromacs ndx error")]
    Ndx(#[from] NdxError),

    #[error("selection expr is not allowed as a definition for subselecting")]
    SelDefInSubsel,

    #[error("can't make par_split from overlapping selections")]
    ParSplitOverlap,

    #[error("can't make par_split from different systems")]
    ParSplitDifferentSystems,
}

/// Errors related to accessing selection indexes
#[derive(Error, Debug)]
pub enum SelectionIndexError {
    #[error("selection index is empty")]
    IndexEmpty,
    #[error("selection indeces {0}:{1} are out of allowed range 0:{2}")]
    IndexOutOfBounds(usize, usize, usize),
}