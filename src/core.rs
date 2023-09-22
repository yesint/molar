use nalgebra::{Vector3,Point3};

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

pub type Vector3f = Vector3<f32>;
pub type Point3f = Point3<f32>;