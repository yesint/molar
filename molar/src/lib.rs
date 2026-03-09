//! MolAR is a library for molecular modeling and analysis written in Rust with an emphasis on memory safety and performance.
//! Molar is designed to simplify the analysis of molecular dynamics trajectories and to implement new analysis algorithms. Molar is intended to provide facilities, which are routinely used in all molecular analysis programs, namely input/output of popular file formats, powerful and flexible atom selections, geometry transformations, RMSD fitting and alignment, etc.
//! MolAR is a logical successor of [Pteros](https://github.com/yesint/pteros) molecular modeling library, which is written in C++ and become hard to develop and maintain due to all C++ idiosyncrasies.

mod aliases;
mod analysis_task;
mod atom;
mod connectivity;
mod distance_search;
mod dssp;
mod measure;
mod modify;
mod ndx_file;
mod particle;
mod periodic_box;
mod periodic_table;
mod providers;
mod sasa;
mod selection;
mod state;
mod topology;
mod seq_align;

pub mod io;
pub mod voronoi_cell;

/// Most useful public imports exposed to the users
pub mod prelude {
    pub use crate::{
        aliases::*, analysis_task::*, atom::*, connectivity::*, distance_search::*, dssp::*,
        io::*, measure::*, modify::*, ndx_file::*, particle::*, periodic_box::*, providers::*,
        sasa::*, selection::*, state::*, topology::*,
    };
    pub use rayon::iter::{IndexedParallelIterator, ParallelIterator};
}

pub use crate::{
    aliases::*, analysis_task::*, atom::*, connectivity::*, distance_search::*, dssp::*, io::*,
    measure::*, modify::*, ndx_file::*, particle::*, periodic_box::*, providers::*, sasa::*,
    selection::*, state::*, topology::*,
};

const BOLD: &str = "\x1b[1m";
const RESET: &str = "\x1b[0m";

/// Prints a welcome message for MolAR with package information and the specified tool name
/// # Example
/// ```
/// use molar::greeting;
/// greeting("analysis");
/// ```
pub fn greeting(tool: impl AsRef<str>) {
    const TITLE: &str = "MolAR - Molecular Analysis for Rust";

    let version = format!("MolAR version: {}", env!("CARGO_PKG_VERSION"));
    let tool    = format!("Tool: {}", tool.as_ref());

    // Content width = widest row; bar = content + one space of padding on each side
    let w = [env!("CARGO_PKG_HOMEPAGE"), env!("CARGO_PKG_AUTHORS"), &version, &tool]
        .iter()
        .map(|s| s.len())
        .max()
        .unwrap_or(0)
        .max(TITLE.len());
    let bar = "─".repeat(w + 2);

    println!("╭{bar}╮");
    println!("│ {BOLD}{TITLE:<w$}{RESET} │");
    println!("├{bar}┤");
    println!("│ {:<w$} │", env!("CARGO_PKG_HOMEPAGE"));
    println!("│ {:<w$} │", env!("CARGO_PKG_AUTHORS"));
    println!("├{bar}┤");
    println!("│ {version:<w$} │");
    println!("├{bar}┤");
    println!("│ {tool:<w$} │");
    println!("╰{bar}╯");
}

// Test code in README
#[cfg(doctest)]
#[doc = include_str!("../../README.md")]
struct _ReadMe;
