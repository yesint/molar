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

pub use {
    atom::*, 
    topology::*,
    state::*,
    periodic_box::*,
    selection::*,
    providers::*,
    measure::*,
    modify::*,
    periodic_table::*,
    particle::*,
    selection_parser::SelectionExpr,
};

// Aliases for vector and points
pub type Vector3f = nalgebra::Vector3<f32>;
pub type Matrix3f = nalgebra::Matrix3<f32>;
pub type Pos = nalgebra::Point3<f32>; // Atom position

// Define alias traits for iterators to make it less verbose
pub trait IndexIterator: ExactSizeIterator<Item = usize> {}
impl<T> IndexIterator for T where T: ExactSizeIterator<Item = usize> {}

pub trait AtomIterator<'a>: Iterator<Item = &'a Atom> {}
impl<'a, T> AtomIterator<'a> for T where T: Iterator<Item = &'a Atom> {}

pub trait AtomMutIterator<'a>: Iterator<Item = &'a mut Atom> {}
impl<'a, T> AtomMutIterator<'a> for T where T: Iterator<Item = &'a mut Atom> {}

pub trait PosIterator<'a>: Iterator<Item = &'a Pos> {}
impl<'a, T> PosIterator<'a> for T where T: Iterator<Item = &'a Pos> {}

pub trait PosMutIterator<'a>: Iterator<Item = &'a mut Pos> {}
impl<'a, T> PosMutIterator<'a> for T where T: Iterator<Item = &'a mut Pos> {}

pub trait IdPosIterator<'a>: ExactSizeIterator<Item = (usize, &'a Pos)> {}
impl<'a, T> IdPosIterator<'a> for T where T: ExactSizeIterator<Item = (usize, &'a Pos)> {}

pub trait IdPosMutIterator<'a>: ExactSizeIterator<Item = (usize, &'a mut Pos)> {}
impl<'a, T> IdPosMutIterator<'a> for T where T: ExactSizeIterator<Item = (usize, &'a mut Pos)> {}