use super::{Atom, Pos};

#[derive(Debug, Clone)]
pub struct Particle<'a> {
    pub id: usize,
    pub atom: &'a Atom,
    pub pos: &'a Pos,
}

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

/// Iterator over particles
pub trait ParticleIterator<'a>: ExactSizeIterator<Item = Particle<'a>> {}
impl<'a, T> ParticleIterator<'a> for T where T: ExactSizeIterator<Item = Particle<'a>> {}

/// Mutable iterator over particles
pub trait ParticleMutIterator<'a>: ExactSizeIterator<Item = ParticleMut<'a>> {}
impl<'a, T> ParticleMutIterator<'a> for T where T: ExactSizeIterator<Item = ParticleMut<'a>> {}

//--------------------------------------------------------
/// Helper struct for creating subscripted  iterator
/// from atoms and positions

#[derive(Clone)]
pub struct ParticleIteratorAdaptor<'a> {
    atom_ref: &'a Vec<Atom>,
    pos_ref: &'a Vec<Pos>,
    index_ref: &'a Vec<usize>,
    cur: usize,
}

impl<'a> ParticleIteratorAdaptor<'a> {
    pub fn new(
        atom_ref: &'a Vec<Atom>,
        pos_ref: &'a Vec<Pos>,
        index_ref: &'a Vec<usize>,
    ) -> Self {
        Self {
            atom_ref,
            pos_ref,
            index_ref,
            cur: 0,
        }
    }
}

impl<'a> Iterator for ParticleIteratorAdaptor<'a> {
    type Item = Particle<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.cur < self.index_ref.len() {
            let ind = self.index_ref[self.cur];
            self.cur += 1;
            Some(Particle {
                atom: &self.atom_ref[ind],
                pos: &self.pos_ref[ind],
                id: self.cur,
            })
        } else {
            None
        }
    }
}

impl<'a> ExactSizeIterator for ParticleIteratorAdaptor<'a> {
    fn len(&self) -> usize {
        self.index_ref.len()
    }
}

//--------------------------------------------------------
/// Helper struct for creating subscripted mutable iterator
/// from iterators over atoms and positions
/// IMPORTANT! Only works for **sorted** indexes!
//--------------------------------------------------------
pub struct ParticleMutIteratorAdaptor<'a> {
    atom_ref: &'a mut Vec<Atom>,
    pos_ref: &'a mut Vec<Pos>,
    index_ref: &'a Vec<usize>,
    cur: usize,
}

impl<'a> ParticleMutIteratorAdaptor<'a> {
    pub fn new(atom_ref: &'a mut Vec<Atom>, pos_ref: &'a mut Vec<Pos>, index_ref: &'a Vec<usize>) -> Self {
        Self {
            atom_ref,
            pos_ref,
            index_ref,
            cur: 0,
        }
    }
}

impl<'a> Iterator for ParticleMutIteratorAdaptor<'a> {
    type Item = ParticleMut<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.cur < self.index_ref.len() {
            let ind = self.index_ref[self.cur];
            self.cur += 1;
            // We need to use unsafe pointers here because the borrow checker
            // won't guarantee that an element is not borrowed as mut
            // multiple time if there are duplicated indexes
            unsafe {
                let p_atom = self.atom_ref.as_mut_ptr().add(ind);
                let p_pos = self.pos_ref.as_mut_ptr().add(ind);
                Some(ParticleMut {
                    atom: &mut *p_atom,
                    pos: &mut *p_pos,
                    id: self.cur-1,
                })
            }
        } else {
            None
        }
    }
}

impl<'a> ExactSizeIterator for ParticleMutIteratorAdaptor<'a> {
    fn len(&self) -> usize {
        self.index_ref.len()
    }
}
