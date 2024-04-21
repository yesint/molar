#![doc = include_str!("../README.md")]

mod core;
mod io;
mod distance_search;

pub mod voro2d {
    pub mod voronoi_cell;
}

pub mod prelude {
    pub use crate::core::*;
    pub use crate::io::*;
    pub use crate::distance_search::*;
}
