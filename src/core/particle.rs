use super::{Atom,Pos, IndexIterator, PosIterator};

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

pub trait ParticleIterator<'a>: ExactSizeIterator<Item = Particle<'a>> {}
impl<'a, T> ParticleIterator<'a> for T where T: ExactSizeIterator<Item = Particle<'a>> {}

pub trait ParticleMutIterator<'a>: ExactSizeIterator<Item = ParticleMut<'a>> {}
impl<'a, T> ParticleMutIterator<'a> for T where T: ExactSizeIterator<Item = ParticleMut<'a>> {}

//--------------------------------------------------------
// Helper struct for creating subscripted  iterator
// from iterators over atoms and positions
// IMPORTANT! Only works for **sorted** indexes!
//--------------------------------------------------------
#[derive(Clone)]
pub struct ParticleIteratorAdaptor<'a, AtomI, PosI, IndexI>
where
    AtomI: Iterator<Item = &'a Atom>, // iterator over atoms
    PosI: Iterator<Item = &'a Pos>,   // iterator over positions
    IndexI: IndexIterator,                // Index iterator
{
    atom_iter: AtomI,
    pos_iter: PosI,
    index_iter: IndexI,
    cur: usize,
}

impl<'a, AtomI, PosI, IndexI> ParticleIteratorAdaptor<'a, AtomI, PosI, IndexI>
where
    AtomI: Iterator<Item = &'a Atom>, // iterator over atoms
    PosI: Iterator<Item = &'a Pos>,   // iterator over positions
    IndexI: IndexIterator,                // Index iterator
{
    pub fn new(atom_iter: AtomI, pos_iter: PosI, index_iter: IndexI) -> Self {
        Self{atom_iter,pos_iter,index_iter,cur: 0}
    }
}

impl<'a, AtomI, PosI, IndexI> Iterator for ParticleIteratorAdaptor<'a, AtomI, PosI, IndexI>
where
    AtomI: Iterator<Item = &'a Atom>, // iterator over atoms
    PosI: Iterator<Item = &'a Pos>,   // iterator over positions
    IndexI: IndexIterator,                // Index iterator
{
    type Item = Particle<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.index_iter.next() {
            Some(id) => {
                // Advance iterators by offset and yield
                let atom = self.atom_iter.nth(id - self.cur)?;
                let pos = self.pos_iter.nth(id - self.cur)?;
                // Advance current position
                self.cur = id + 1;
                Some(Particle { atom, pos, id })
            }
            None => None,
        }
    }
}

impl<'a, AtomI, PosI, IndexI> ExactSizeIterator
    for ParticleIteratorAdaptor<'a, AtomI, PosI, IndexI>
where
    AtomI: Iterator<Item = &'a Atom>, // iterator over atoms
    PosI: Iterator<Item = &'a Pos>,   // iterator over positions
    IndexI: IndexIterator,                // Index iterator
{
    fn len(&self) -> usize {
        self.index_iter.len()
    }
}

//--------------------------------------------------------
// Helper struct for creating subscripted mutable iterator
// from iterators over atoms and positions
// IMPORTANT! Only works for **sorted** indexes!
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
        Self{atom_iter,pos_iter,index_iter,cur: 0}
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
        match self.index_iter.next() {
            Some(id) => {
                // Advance iterators by offset and yield
                let atom = self.atom_iter.nth(id - self.cur)?;
                let pos = self.pos_iter.nth(id - self.cur)?;
                // Advance current position
                self.cur = id + 1;
                Some(ParticleMut { atom, pos, id })
            }
            None => None,
        }
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