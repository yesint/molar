//! Core data structures and traits: [Atom], [State], [Topology], [selections](Sel), distance search, etc.

mod atom;
mod topology;
mod state;
mod periodic_box;
mod selection;
mod providers;
mod measure;
mod modify;
mod periodic_table;
mod particle;
#[allow(dead_code)]
mod selection_parser;
mod distance_search;
mod connectivity;
mod ndx_file;
//mod system;

pub use {
    atom::*, 
    topology::*,
    state::*,
    periodic_box::*,
    selection::*,
    providers::*,
    measure::*,
    modify::*,
    particle::*,
    selection_parser::SelectionExpr,
    distance_search::*,
    connectivity::*,
    ndx_file::*,
    
};

// Aliases for vector and points
/// Convenience alias for 3D vector
pub type Vector3f = nalgebra::Vector3<f32>;
/// Convenience alias for 3x3 matrix
pub type Matrix3f = nalgebra::Matrix3<f32>;
/// Atom position
pub type Pos = nalgebra::Point3<f32>;

// Define alias traits for iterators to make it less verbose
/// Convenience alias for iterator over indices
pub trait IndexIterator: Iterator<Item = usize> {}
impl<T> IndexIterator for T where T: Iterator<Item = usize> {}

/// Convenience alias for iterator over atoms
pub trait AtomIterator<'a>: Iterator<Item = &'a Atom> {}
impl<'a, T> AtomIterator<'a> for T where T: Iterator<Item = &'a Atom> {}

/// Convenience alias for mutable iterator over atoms
pub trait AtomMutIterator<'a>: Iterator<Item = &'a mut Atom> {}
impl<'a, T> AtomMutIterator<'a> for T where T: Iterator<Item = &'a mut Atom> {}

/// Convenience alias for iterator over positions
pub trait PosIterator<'a>: Iterator<Item = &'a Pos> + Clone {}
impl<'a, T> PosIterator<'a> for T where T: Iterator<Item = &'a Pos> + Clone {}

/// Convenience alias for mutable iterator over positions  
pub trait PosMutIterator<'a>: Iterator<Item = &'a mut Pos> {}
impl<'a, T> PosMutIterator<'a> for T where T: Iterator<Item = &'a mut Pos> {}

/// Convenience alias for iterator over index-position pairs
pub trait IdPosIterator<'a>: Iterator<Item = (usize, &'a Pos)> {}
impl<'a, T> IdPosIterator<'a> for T where T: Iterator<Item = (usize, &'a Pos)> {}

/// Convenience alias for mutable iterator over index-position pairs
pub trait IdPosMutIterator<'a>: Iterator<Item = (usize, &'a mut Pos)> {}
impl<'a, T> IdPosMutIterator<'a> for T where T: Iterator<Item = (usize, &'a mut Pos)> {}