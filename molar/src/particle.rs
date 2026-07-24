use crate::prelude::*;

/// Holds a read-only column proxy for an atom together with its position and particle id.
#[derive(Debug)]
pub struct Particle<'a> {
    pub id: usize,
    pub atom: AtomRef<'a>,
    pub pos: &'a Pos,
}

/// Holds a mutable column proxy for an atom together with its position and particle id.
#[derive(Debug)]
pub struct ParticleMut<'a> {
    pub id: usize,
    pub atom: AtomRefMut<'a>,
    pub pos: &'a mut Pos,
}

impl<'a> From<ParticleMut<'a>> for Particle<'a> {
    fn from(p: ParticleMut<'a>) -> Self {
        Particle {
            id: p.id,
            atom: p.atom.into(),
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

impl IndexProvider for Particle<'_> {
    unsafe fn get_index_unchecked(&self, i: usize) -> usize {
        if i > 0 {
            panic!("single particle can only be accessed with id=0, not {i}")
        } else {
            self.id
        }
    }

    fn iter_index(&self) -> impl ExactSizeIterator<Item = usize> {
        std::iter::once(self.id)
    }
}

impl PosProvider for Particle<'_> {
    unsafe fn coords_ptr(&self) -> *const Pos {
        self.pos as *const Pos
    }

    // Override iter_pos to avoid index indirection (particle has a direct pos ptr)
    fn iter_pos(&self) -> impl PosIterator<'_> {
        std::iter::once(self.pos)
    }

    unsafe fn get_pos_unchecked(&self, _i: usize) -> &Pos {
        self.pos
    }
}

impl AtomProvider for Particle<'_> {
    fn atom_storage(&self) -> &AtomStorage {
        self.atom.storage()
    }

    // Override iter_atoms to avoid index indirection
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        std::iter::once(self.atom)
    }

    unsafe fn get_atom_unchecked(&self, _i: usize) -> AtomRef<'_> {
        self.atom
    }
}

//-----------------------------------------------------------------
