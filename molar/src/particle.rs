use crate::prelude::*;

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
        Particle {
            id: p.id,
            atom: p.atom,
            pos: p.pos,
        }
    }
}

//------------------------------------------------------

impl LenProvider for Particle<'_> {
    fn len(&self) -> usize {
        1
    }
}

impl PosIterProvider for Particle<'_> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        std::iter::once(self.pos)
    }
}

impl AtomIterProvider for Particle<'_> {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        std::iter::once(self.atom)
    }
}

impl IndexProvider for Particle<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        if i > 0 {
            panic!("single particle can only be accessed with id=0, not {i}")
        } else {
            self.id
        }
    }

    fn iter_index(&self) -> impl Iterator<Item = usize> + Clone {
        std::iter::once(self.id)
    }
}

impl ParticleIterProvider for Particle<'_> {
    fn iter_particle(&self) -> impl Iterator<Item = Particle<'_>> {
        std::iter::once(Particle {
            id: self.id,
            pos: self.pos,
            atom: self.atom,
        })
    }
}

//-----------------------------------------------------------------
