mod atom;
mod topology;
mod state;
mod periodic_box;
#[allow(dead_code)]
mod selection_parser;
mod selection;
mod algorithms;
mod particle;

pub use {
    atom::Atom, 
    topology::*, 
    state::*,
    periodic_box::*,
    algorithms::*,
    selection::*,
    particle::*,
}; 

// Aliases for vector and points
pub type Vector3f = nalgebra::Vector3<f32>;
pub type Matrix3f = nalgebra::Matrix3<f32>;
pub type Pos = nalgebra::Point3<f32>; // Atom position

// Define alias traits for iterators to make it less verbose
pub trait IndexIterator: ExactSizeIterator<Item = usize> {}
impl<T> IndexIterator for T where T: ExactSizeIterator<Item = usize> {}

pub trait PosIterator<'a>: ExactSizeIterator<Item = &'a Pos> {}
impl<'a, T> PosIterator<'a> for T where T: ExactSizeIterator<Item = &'a Pos> {}

pub trait IdPosIterator<'a>: ExactSizeIterator<Item = (usize, &'a Pos)> {}
impl<'a, T> IdPosIterator<'a> for T where T: ExactSizeIterator<Item = (usize, &'a Pos)> {}

pub trait IdPosIteratorMut<'a>: ExactSizeIterator<Item = (usize, &'a mut Pos)> {}
impl<'a, T> IdPosIteratorMut<'a> for T where T: ExactSizeIterator<Item = (usize, &'a mut Pos)> {}
