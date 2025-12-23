//! MolAR is a library for molecular modeling and analysis written in Rust with an emphasis on memory safety and performance.
//! Molar is designed to simplify the analysis of molecular dynamics trajectories and to implement new analysis algorithms. Molar is intended to provide facilities, which are routinely used in all molecular analysis programs, namely input/output of popular file formats, powerful and flexible atom selections, geometry transformations, RMSD fitting and alignment, etc.
//! MolAR is a logical successor of [Pteros](https://github.com/yesint/pteros) molecular modeling library, which is written in C++ and become hard to develop and maintain due to all C++ idiosyncrasies.

pub mod aliases;
pub mod analysis_task;
pub mod io;
pub mod voronoi_cell;

mod atom;
mod connectivity;
mod distance_search;
mod measure;
mod modify;
mod ndx_file;
mod particle;
mod periodic_box;
mod periodic_table;
mod providers;
mod selection;
mod state;
mod topology;

/// Most useful public imports exposed to the users
pub mod prelude {
    pub use rayon::iter::IndexedParallelIterator;
    pub use rayon::iter::ParallelIterator;
    pub use crate::{
        aliases::*, analysis_task::*, atom::*, connectivity::*,
        distance_search::*, io::*, measure::*, modify::*,
        ndx_file::*, particle::*, periodic_box::*, providers::*,
        selection::*, state::*,
        topology::*,
    };
}

/// Prints a welcome message for MolAR with package information and the specified tool name
/// # Example
/// ```
/// use molar::greeting;
/// greeting("analysis");
/// ```
pub fn greeting(tool: impl AsRef<str>) {
    use comfy_table::modifiers::UTF8_ROUND_CORNERS;
    use comfy_table::presets::UTF8_FULL;
    use comfy_table::{Attribute, Cell, Table};

    let mut table = Table::new();
    table
        .load_preset(UTF8_FULL)
        .apply_modifier(UTF8_ROUND_CORNERS)
        .add_row(vec![
            Cell::new("MolAR - Molecular Analysis for Rust").add_attributes(vec![Attribute::Bold])
        ])
        .add_row(vec![format!(
            "{}\n{}",
            env!("CARGO_PKG_HOMEPAGE"),
            env!("CARGO_PKG_AUTHORS")
        )])
        .add_row(vec![format!(
            "MolAR version: {}",
            env!("CARGO_PKG_VERSION")
        )])
        .add_row(vec![format!("Tool: {}", tool.as_ref())]);
    println!("{table}");
}

// Test code in README
#[cfg(doctest)]
#[doc = include_str!("../../README.md")]
struct _ReadMe;
