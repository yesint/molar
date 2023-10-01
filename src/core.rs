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

pub trait PosIteratorMut<'a>: ExactSizeIterator<Item = &'a mut Pos> {}
impl<'a,T> PosIteratorMut<'a> for T where T: ExactSizeIterator<Item = &'a mut Pos> {}

pub trait IdPosIterator<'a>: ExactSizeIterator<Item = (usize, &'a Pos)> {}
impl<'a, T> IdPosIterator<'a> for T where T: ExactSizeIterator<Item = (usize, &'a Pos)> {}

pub trait IdPosIteratorMut<'a>: ExactSizeIterator<Item = (usize, &'a mut Pos)> {}
impl<'a, T> IdPosIteratorMut<'a> for T where T: ExactSizeIterator<Item = (usize, &'a mut Pos)> {}

pub trait AtomIterator<'a>: ExactSizeIterator<Item = &'a Atom> {}
impl<'a, T> AtomIterator<'a> for T where T: ExactSizeIterator<Item = &'a Atom> {}

pub trait AtomIteratorMut<'a>: ExactSizeIterator<Item = &'a mut Atom> {}
impl<'a, T> AtomIteratorMut<'a> for T where T: ExactSizeIterator<Item = &'a mut Atom> {}

pub trait IdAtomIterator<'a>: ExactSizeIterator<Item = (usize, &'a Atom)> {}
impl<'a, T> IdAtomIterator<'a> for T where T: ExactSizeIterator<Item = (usize, &'a Atom)> {}

pub trait IdAtomIteratorMut<'a>: ExactSizeIterator<Item = (usize, &'a mut Atom)> {}
impl<'a, T> IdAtomIteratorMut<'a> for T where T: ExactSizeIterator<Item = (usize, &'a mut Atom)> {}

pub trait IdAtomPosIterator<'a>: ExactSizeIterator<Item = (usize, &'a Atom, &'a Pos)> {}
impl<'a, T> IdAtomPosIterator<'a> for T where T: ExactSizeIterator<Item = (usize, &'a Atom, &'a Pos)> {}

pub trait IdAtomPosIteratorMut<'a>: ExactSizeIterator<Item = (usize, &'a mut Atom, &'a mut Pos)> {}
impl<'a, T> IdAtomPosIteratorMut<'a> for T where T: ExactSizeIterator<Item = (usize, &'a mut Atom, &'a mut Pos)> {}
