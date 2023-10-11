use std::marker::PhantomData;

use super::{
    selection_parser::{apply_ast_whole, generate_ast, SelectionAst},
    Atom, AtomIterator, AtomIteratorMut, IdAtomIterator, IdAtomIteratorMut, IdAtomPosIterator,
    IdAtomPosIteratorMut, IdPosIterator, IdPosIteratorMut, IndexIterator, Pos, PosIterator,
    PosIteratorMut, State, Structure,
};
use anyhow::{bail, Result};
use nalgebra::Point3;

// There are following distinct states of the Selection:
//  Selection::from_iter(0..10).iter(structure,state)
//  Selection::from_expr("name_CA").iter(structure,state)

struct Particle<'a> {
    id: usize,
    atom: &'a Atom,
    pos: &'a Pos,
}

struct ParticleMut<'a> {
    id: usize,
    atom: &'a mut Atom,
    pos: &'a mut Pos,
}


pub trait ParticleIterator<'a>: ExactSizeIterator<Item = Particle<'a>> {}
impl<'a, T> ParticleIterator<'a> for T where T: ExactSizeIterator<Item = Particle<'a>> {}

pub trait ParticleMutIterator<'a>: ExactSizeIterator<Item = ParticleMut<'a>> {}
impl<'a, T> ParticleMutIterator<'a> for T where T: ExactSizeIterator<Item = ParticleMut<'a>> {}


struct Selection {}
struct SelectionWithAst(SelectionAst);
struct SelectionWithIter<I: IndexIterator>(I);

impl Selection {
    fn from_expr(sel_str: &str) -> Result<SelectionWithAst> {
        let ast = generate_ast(sel_str)?;
        Ok(SelectionWithAst(ast))
    }

    fn from_iter<I: IndexIterator>(iter: I) -> SelectionWithIter<I> {
        SelectionWithIter(iter)
    }
}

impl SelectionWithAst {
    fn apply<'a>(&self,structure: &'a Structure, state: &'a State) -> Result<impl ParticleIterator<'a>> {
        let vec = apply_ast_whole(&self.0, structure, state)?;
        Ok(vec.into_iter().map(|i| Particle{
            id: i,
            atom: &structure.atoms[i],
            pos: &state.coords[i]
        }))
    }

    fn apply_mut<'a>(&self,structure: &'a mut Structure, state: &'a mut State) -> Result<impl ParticleMutIterator<'a>> {
        let vec = apply_ast_whole(&self.0, structure, state)?;
        Ok(ParticleMutIteratorAdaptor {
            atom_iter: structure.atoms.iter_mut(),
            pos_iter: state.coords.iter_mut(),
            index_iter: vec.into_iter(),
            cur: 0,
        })
    }
}

impl<I: IndexIterator> SelectionWithIter<I> {
    fn apply<'a>(&self,structure: &'a Structure, state: &'a State) -> Result<impl ParticleIterator<'a>> {
        Ok(self.0.clone().map(|i| Particle{
            id: i,
            atom: &structure.atoms[i],
            pos: &state.coords[i]
        }))
    }

    fn apply_mut<'a>(&self,structure: &'a mut Structure, state: &'a mut State) -> Result<impl ParticleMutIterator<'a>> {
        Ok(ParticleMutIteratorAdaptor {
            atom_iter: structure.atoms.iter_mut(),
            pos_iter: state.coords.iter_mut(),
            index_iter: self.0.clone(),
            cur: 0,
        })
    }
}

//---------------------------------------


// Helper struct for creating subscripted mutable iterator
// from container iterator and index iterator
// IMPORTANT! Only works for **sorted** indexes!
struct ParticleMutIteratorAdaptor<'a, AtomI, PosI, IndexI>
where
    AtomI: Iterator<Item = &'a mut Atom>, // iterator over atoms
    PosI: Iterator<Item = &'a mut Pos>, // iterator over positions
    IndexI: IndexIterator,              // Index iterator
{
    atom_iter: AtomI,
    pos_iter: PosI,
    index_iter: IndexI,
    cur: usize,
}

impl<'a, AtomI, PosI, IndexI> Iterator for ParticleMutIteratorAdaptor<'a, AtomI, PosI, IndexI>
where
    AtomI: Iterator<Item = &'a mut Atom>, // iterator over atoms
    PosI: Iterator<Item = &'a mut Pos>, // iterator over positions
    IndexI: IndexIterator,              // Index iterator
{
    type Item = ParticleMut<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.index_iter.next() {
            Some(i) => {
                // Advance pos_iter by offset and yield
                let at = self.atom_iter.nth(i - self.cur)?;
                let pos = self.pos_iter.nth(i - self.cur)?;
                // Advance current position
                self.cur = i + 1;
                Some(ParticleMut{
                    atom: at,
                    pos: pos,
                    id: i,
                })
            }
            None => None,
        }
    }
}

impl<'a, AtomI, PosI, IndexI> ExactSizeIterator for ParticleMutIteratorAdaptor<'a, AtomI, PosI, IndexI>
where
    AtomI: Iterator<Item = &'a mut Atom>, // iterator over atoms
    PosI: Iterator<Item = &'a mut Pos>, // iterator over positions
    IndexI: IndexIterator,              // Index iterator
{
    fn len(&self) -> usize {
        self.index_iter.len()
    }
}
