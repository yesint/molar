mod atom;
mod structure;
mod state;
mod periodic_box;
#[allow(dead_code)]
mod selection;

pub use {
    atom::Atom, 
    structure::Structure, 
    state::State,
    periodic_box::*,
}; 

// Aliases for vector and points
pub type Vector3f = nalgebra::Vector3<f32>;
pub type Matrix3f = nalgebra::Matrix3<f32>;
pub type Pos = nalgebra::Point3<f32>; // Atom position

// Define alias traits for iterators to make it less verbose
pub trait IndexIterator: ExactSizeIterator<Item = usize> {}
impl<T> IndexIterator for T where T: ExactSizeIterator<Item = usize> {}

pub trait PosIterator<'a>: ExactSizeIterator<Item = &'a Pos> {}
impl<'a,T> PosIterator<'a> for T where T: ExactSizeIterator<Item = &'a Pos> {}