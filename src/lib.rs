#![cfg_attr(not(doctest), doc = include_str!("../README.md"))]

pub mod core;
pub mod io;
pub mod distance_search;
pub mod voronoi_cell;

pub mod prelude {
    pub use crate::core::*;
    pub use crate::io::*;
    pub use crate::distance_search::*;
    pub use rayon::iter::ParallelIterator;
}

pub const VERSION: &str = env!("CARGO_PKG_VERSION");