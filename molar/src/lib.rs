//#![doc = include_str!("../../README.md")]

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


pub fn greeting(tool: &str) {
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
        .add_row(vec![format!("MolAR version: {}", env!("CARGO_PKG_VERSION"))])
        .add_row(vec![format!(
            "Tool: {tool}"
        )]);
    println!("{table}");
}
