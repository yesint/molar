use crate::prelude::*;

/// Alias to sorted vector
pub type SVec = sorted_vec::SortedSet<usize>;

/// Scalar floating-point type used throughout MolAR.
///
/// Defaults to [`f32`]. When the `f64` feature is enabled, the entire
/// workspace switches to [`f64`] precision.
#[cfg(not(feature = "f64"))]
pub type Float = f32;
#[cfg(feature = "f64")]
pub type Float = f64;

// Aliases for vector and points
/// Convenience alias for 2D vector
pub type Vector2f = nalgebra::Vector2<Float>;
/// Convenience alias for 3D vector
pub type Vector3f = nalgebra::Vector3<Float>;
/// Convenience alias for 3x3 matrix
pub type Matrix3f = nalgebra::Matrix3<Float>;
/// Atom position
pub type Pos = nalgebra::Point3<Float>;
/// Atom velocity (nm/ps)
pub type Vel = nalgebra::Vector3<Float>;
/// Atom force (kJ/mol/nm)
pub type Force = nalgebra::Vector3<Float>;

// Define alias traits for iterators to make it less verbose
/// Convenience alias for iterator over indices
pub trait IndexIterator: ExactSizeIterator<Item = usize> {}
impl<T> IndexIterator for T where T: ExactSizeIterator<Item = usize> {}

/// Convenience alias for iterator over atoms
pub trait AtomIterator<'a>: ExactSizeIterator<Item = &'a Atom> {}
impl<'a, T> AtomIterator<'a> for T where T: ExactSizeIterator<Item = &'a Atom> {}

/// Convenience alias for mutable iterator over atoms
pub trait AtomMutIterator<'a>: ExactSizeIterator<Item = &'a mut Atom> {}
impl<'a, T> AtomMutIterator<'a> for T where T: ExactSizeIterator<Item = &'a mut Atom> {}

/// Convenience alias for iterator over positions
pub trait PosIterator<'a>: ExactSizeIterator<Item = &'a Pos> {}
impl<'a, T> PosIterator<'a> for T where T: ExactSizeIterator<Item = &'a Pos> {}

/// Convenience alias for mutable iterator over positions
pub trait PosMutIterator<'a>: ExactSizeIterator<Item = &'a mut Pos> {}
impl<'a, T> PosMutIterator<'a> for T where T: ExactSizeIterator<Item = &'a mut Pos> {}

/// Convenience alias for iterator over index-position pairs
pub trait IdPosIterator<'a>: ExactSizeIterator<Item = (usize, &'a Pos)> {}
impl<'a, T> IdPosIterator<'a> for T where T: ExactSizeIterator<Item = (usize, &'a Pos)> {}

/// Convenience alias for mutable iterator over index-position pairs
pub trait IdPosMutIterator<'a>: ExactSizeIterator<Item = (usize, &'a mut Pos)> {}
impl<'a, T> IdPosMutIterator<'a> for T where T: ExactSizeIterator<Item = (usize, &'a mut Pos)> {}

/// Convenience alias for iterator over particles
pub trait ParticleIterator<'a>: ExactSizeIterator<Item = Particle<'a>> {}
impl<'a, T> ParticleIterator<'a> for T where T: ExactSizeIterator<Item = Particle<'a>> {}

/// Convenience alias for mutable iterator over particles
pub trait ParticleMutIterator<'a>: ExactSizeIterator<Item = ParticleMut<'a>> {}
impl<'a, T> ParticleMutIterator<'a> for T where T: ExactSizeIterator<Item = ParticleMut<'a>> {}

/// Convenience alias for iterator over velocities
pub trait VelIterator<'a>: ExactSizeIterator<Item = &'a Vel> {}
impl<'a, T> VelIterator<'a> for T where T: ExactSizeIterator<Item = &'a Vel> {}

/// Convenience alias for mutable iterator over velocities
pub trait VelMutIterator<'a>: ExactSizeIterator<Item = &'a mut Vel> {}
impl<'a, T> VelMutIterator<'a> for T where T: ExactSizeIterator<Item = &'a mut Vel> {}

/// Convenience alias for iterator over forces
pub trait ForceIterator<'a>: ExactSizeIterator<Item = &'a Force> {}
impl<'a, T> ForceIterator<'a> for T where T: ExactSizeIterator<Item = &'a Force> {}

/// Convenience alias for mutable iterator over forces
pub trait ForceMutIterator<'a>: ExactSizeIterator<Item = &'a mut Force> {}
impl<'a, T> ForceMutIterator<'a> for T where T: ExactSizeIterator<Item = &'a mut Force> {}