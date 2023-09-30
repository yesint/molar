mod atom;
mod structure;
mod state;
mod periodic_box;
#[allow(dead_code)]
mod selection_parser;
mod selection;

pub use {
    atom::Atom, 
    structure::*, 
    state::*,
    periodic_box::*,
}; 

// Aliases for vector and points
pub type Vector3f = nalgebra::Vector3<f32>;
pub type Matrix3f = nalgebra::Matrix3<f32>;
pub type Pos = nalgebra::Point3<f32>; // Atom position

// Define alias traits for iterators to make it less verbose
pub trait IndexIterator: ExactSizeIterator<Item = usize> + Clone {}
impl<T> IndexIterator for T where T: ExactSizeIterator<Item = usize> + Clone {}

pub trait PosIterator<'a>: ExactSizeIterator<Item = &'a Pos> + Clone {}
impl<'a,T> PosIterator<'a> for T where T: ExactSizeIterator<Item = &'a Pos> + Clone {}