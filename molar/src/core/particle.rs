use super::{Atom, Pos};

/// Holds immutable reference to [Atom] and [Pos] and particle id.
/// Usually created indirectly by types implementing [ParticleIterProvider](crate::core::ParticleIterProvider).
#[derive(Debug)]
pub struct Particle<'a> {
    pub id: usize,
    pub atom: &'a Atom,
    pub pos: &'a Pos,
}

/// Holds mutable reference to [Atom] and [Pos] and particle id.
/// Usually created indirectly by types implementing [ParticleIterMutProvider](crate::core::ParticleIterMutProvider).
#[derive(Debug)]
pub struct ParticleMut<'a> {
    pub id: usize,
    pub atom: &'a mut Atom,
    pub pos: &'a mut Pos,
}

impl<'a> From<ParticleMut<'a>> for Particle<'a> {
    fn from(p: ParticleMut<'a>) -> Self {
        Particle{
            id: p.id,
            atom: p.atom,
            pos: p.pos,
        }
    }
}