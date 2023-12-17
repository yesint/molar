use super::{Atom, IndexIterator, Pos};

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
        if self.cur >= self.index_ref.len() {
            None
        } else {
            let ind = self.index_ref[self.cur];
            self.cur += 1;
            Some(Particle {
                atom: &self.atom_ref[ind],
                pos: &self.pos_ref[ind],
                id: ind,
            })
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
#[derive(Clone)]
pub struct ParticleMutIteratorAdaptor<'a, AtomI, PosI, IndexI>
where
    AtomI: Iterator<Item = &'a mut Atom>, // iterator over atoms
    PosI: Iterator<Item = &'a mut Pos>,   // iterator over positions
    IndexI: IndexIterator,                // Index iterator
{
    atom_iter: AtomI,
    pos_iter: PosI,
    index_iter: IndexI,
    cur: usize,
}

impl<'a, AtomI, PosI, IndexI> ParticleMutIteratorAdaptor<'a, AtomI, PosI, IndexI>
where
    AtomI: Iterator<Item = &'a mut Atom>, // iterator over atoms
    PosI: Iterator<Item = &'a mut Pos>,   // iterator over positions
    IndexI: IndexIterator,                // Index iterator
{
    pub fn new(atom_iter: AtomI, pos_iter: PosI, index_iter: IndexI) -> Self {
        Self {
            atom_iter,
            pos_iter,
            index_iter,
            cur: 0,
        }
    }
}

impl<'a, AtomI, PosI, IndexI> Iterator for ParticleMutIteratorAdaptor<'a, AtomI, PosI, IndexI>
where
    AtomI: Iterator<Item = &'a mut Atom>, // iterator over atoms
    PosI: Iterator<Item = &'a mut Pos>,   // iterator over positions
    IndexI: IndexIterator,                // Index iterator
{
    type Item = ParticleMut<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        self.index_iter.next().map(|ind| {
            // Advance iterators by offset and yield
            let atom = self.atom_iter.nth(ind - self.cur)?;
            let pos = self.pos_iter.nth(ind - self.cur)?;
            // Advance current position
            self.cur = ind + 1;
            Some(ParticleMut { atom, pos, id: ind })
        })?
    }
}

impl<'a, AtomI, PosI, IndexI> ExactSizeIterator
    for ParticleMutIteratorAdaptor<'a, AtomI, PosI, IndexI>
where
    AtomI: Iterator<Item = &'a mut Atom>, // iterator over atoms
    PosI: Iterator<Item = &'a mut Pos>,   // iterator over positions
    IndexI: IndexIterator,                // Index iterator
{
    fn len(&self) -> usize {
        self.index_iter.len()
    }
}
