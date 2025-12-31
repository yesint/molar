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

    #[error("local index {0} is beyond the allowed range 0:{1}")]
    LocalToGlobal(usize, usize),

    #[error("selection {0}:{1} is out of the source bounds 0:{2}")]
    IndexValidation(usize, usize, usize),

    #[error(transparent)]
    FileIo(#[from] FileIoError),

    #[error(transparent)]
    Builder(#[from] BuilderError),

    #[error("can't set incompatible state")]
    IncompatibleState,

    #[error("can't set incompatible topology")]
    IncompatibleTopology,

    #[error("selection slice is empty")]
    EmptySlice,

    #[error("selection range is empty")]
    EmptyRange,

    #[error("selection '{0}' is empty")]
    EmptyExpr(String),

    #[error("splitting produced no selections")]
    EmptySplit,

    #[error("selection intersection is empy")]
    EmptyIntersection,

    #[error("selection difference is empy")]
    EmptyDifference,

    #[error("selection inversion is empy")]
    EmptyInvert,

    #[error(transparent)]
    PeriodicBox(#[from] PeriodicBoxError),

    #[error("gromacs ndx error")]
    Ndx(#[from] NdxError),

    #[error("selection expr is not allowed as a definition for subselecting")]
    SelDefInSubsel,
}
