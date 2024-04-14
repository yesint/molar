//! # MolAR
//! **MolAR** is a **Mol**ecular **A**nalysis and modeling library for **R**ust.
//! 
//! ## What is Molar?
//! Molar is a rust library for molecular analysis and modeling. 
//! It is started as a successor of [Pteros](https://github.com/yesint/pteros) 
//! molecular modeling library, which is written in C++ and become hard to 
//! develop and maintain due to all C++ idiosyncrasies. 
//! Eventually molar may become a "Pteros 3.0".
//! 
//! Molar is designed to simplify the analysis of molecular dynamics trajectories 
//! and to implement new analysis algorithms. 
//! Molar is intended to provide facilities, which are routinely used in all
//!  molecular analysis programs, namely input/output of popular file formats, 
//! powerful and flexible atom selections, 
//! geometry transformations, RMSD fitting and alignment, etc.
//! 
//! ## Features
//! * Reading and writing PDB, GRO, XTC files.
//! * Reading Gromacs TPR files.
//! * Selections using the syntaxis similar to VMD and Pteros.
//! * Subselections and splitting selections.
//! * SASA calculations with the fastest PowerSasa method.
//! * RMSD fitting and alignment
//! * Basic algorithm (center of mass, gyration radius, etc.)
//! * Seamless PBC treatment
//! * Automatic seamless compiling and linking with VMD molfile,
//! xdrfile libraries and the git version of Gromacs.
//! 
//! # Design decisions
//! Molecular analysis typically involves multiple views of the arrays of atoms and
//! coordinates (aka Selections). Selections could overlap arbitraryly - for example
//! one selection may represent the whole protein and the other only the resiudes of its
//! active site. Selections are mutable in the sense that their atoms could be 
//! manipulated - translated, rotated, aligned, renamed, etc. Such cnahges are expected 
//! to be immediately picked up by all other selections involving affected atoms. 
//! 
//! This concept doesn't play well with Rust ownership rules, where either a single
//! exclusive reference (`mut&`) _or_ multiple immutable references (`&`) could exists
//! at the same time. If one selection holds a `mut&` of the underlying array of atoms 
//! then no other selection is allowed to access it neither for reading nor for writing.
//! 
//! If we want to obey the single ownership rule we need to create and drop selections 
//! for every single operation to make sure that the atoms array is never aliased mutably
//! by several selections. This appears to be extremely inconveniet in practice and
//! contrasts with _all_ existing molecular analysis software.
//! 
//! Thus Molar uses unsafe rust and interior mutability pattern 
//! internally to sidestep aliasing restriction while still having safe user 
//! accessible API.
//! In Molar one can have multiple immutable references `&Sel` to selections
//! while still being able to invoke selection methods that mutate underlying arrays
//! of atoms and coordinates, such as `sel.translate()` or `sel.rotate()`. Changes
//! made by one selection are immediately visible to all other selections that
//! point to the same atoms.
//! 
//! ### Safety guarantees
//! For each distinc pair of [Topology](crate::core::Topology) and [State](crate::core::State) _exactly one_
//! access mode for selection may exists at any point in program runtime:
//! 1. Overlapping mutable sequential
//!     * Overlapping selections may only mutate atoms and coordinates sequentially
//!     from the same thread where they are created.
//! 1. Non-overlapping mutable parallel
//!     * Non-overlapping selections can mutate atoms and coordinates in parallel 
//!     from multiple threads. 
//!     * All parallel processing threads are guaranteed to be joined 
//!     and synchronized before the next processing block starts.
//!     * Any sub-selection made from such selection can only leave in the same
//!     thread as a parent selection.
//!     * Data races are not possible.
//! 1. Overlapping immutable parallel
//!     * Read-only access to atoms and coordinates is possible in parallel from
//!     multiple threads.

pub mod core;
pub mod io;
pub mod distance_search;

pub mod voro2d {
    pub mod voronoi_cell;
}
