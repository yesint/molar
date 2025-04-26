pub mod core;
pub mod io;
pub mod voronoi_cell;
pub mod analysis_task;

pub mod prelude {
    pub use crate::core::*;
    pub use crate::io::*;
    pub use crate::analysis_task::*;
    pub use rayon::iter::ParallelIterator;
}

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
        .add_row(vec![format!("MolAR version: {}", env!("CARGO_PKG_VERSION"))])
        .add_row(vec![format!(
            "Tool: {}",tool.as_ref()
        )]);
    println!("{table}");
}

// Test code in README
#[cfg(doctest)]
#[doc = include_str!("../../README.md")]
struct _ReadMe;