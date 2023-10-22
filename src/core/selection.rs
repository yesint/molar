use std::borrow::Cow;

use crate::distance_search::search::{self, SearchConnectivity, SearcherSingleGrid};

use super::{
    selection_parser::{apply_ast_whole, generate_ast, SelectionAst},
    Atom, IndexIterator, PbcDims, PeriodicBox, Pos, State, Structure,
};
use anyhow::{anyhow, Result};

#[derive(Debug, Clone)]
pub struct Particle<'a> {
    pub id: usize,
    pub atom: &'a Atom,
    pub pos: &'a Pos,
}

pub struct ParticleMut<'a> {
    pub id: usize,
    pub atom: &'a mut Atom,
    pub pos: &'a mut Pos,
}

pub trait ParticleIterator<'a>: ExactSizeIterator<Item = Particle<'a>> + Clone {}
impl<'a, T> ParticleIterator<'a> for T where T: ExactSizeIterator<Item = Particle<'a>> + Clone {}

pub trait ParticleMutIterator<'a>: ExactSizeIterator<Item = ParticleMut<'a>> {}
impl<'a, T> ParticleMutIterator<'a> for T where T: ExactSizeIterator<Item = ParticleMut<'a>> {}

//-----------------------------------------------------------------

struct SelectionWithAst(SelectionAst);
struct SelectionWithIter<I: IndexIterator>(I);
struct SelectionAll {}

impl SelectionWithAst {
    fn apply<'a>(
        &self,
        structure: &'a Structure,
        state: &'a State,
    ) -> Result<impl ParticleIterator<'a>> {
        let vec = apply_ast_whole(&self.0, structure, state)?;
        Ok(vec.into_iter().map(|i| Particle {
            id: i,
            atom: &structure.atoms[i],
            pos: &state.coords[i],
        }))
    }

    fn apply_mut<'a>(
        &self,
        structure: &'a mut Structure,
        state: &'a mut State,
    ) -> Result<impl ParticleMutIterator<'a>> {
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
    fn apply<'a>(
        &self,
        structure: &'a Structure,
        state: &'a State,
    ) -> Result<impl ParticleIterator<'a>> {
        Ok(self.0.clone().map(|i| Particle {
            id: i,
            atom: &structure.atoms[i],
            pos: &state.coords[i],
        }))
    }

    fn apply_mut<'a>(
        &self,
        structure: &'a mut Structure,
        state: &'a mut State,
    ) -> Result<impl ParticleMutIterator<'a>> {
        Ok(ParticleMutIteratorAdaptor {
            atom_iter: structure.atoms.iter_mut(),
            pos_iter: state.coords.iter_mut(),
            index_iter: self.0.clone(),
            cur: 0,
        })
    }
}

impl SelectionAll {
    fn apply<'a>(
        &self,
        structure: &'a Structure,
        state: &'a State,
    ) -> Result<impl ParticleIterator<'a>> {
        Ok((0..state.coords.len()).map(|i| Particle {
            id: i,
            atom: &structure.atoms[i],
            pos: &state.coords[i],
        }))
    }

    fn apply_mut<'a>(
        &self,
        structure: &'a mut Structure,
        state: &'a mut State,
    ) -> Result<impl ParticleMutIterator<'a>> {
        let n = state.coords.len();
        Ok(ParticleMutIteratorAdaptor {
            atom_iter: structure.atoms.iter_mut(),
            pos_iter: state.coords.iter_mut(),
            index_iter: 0..n,
            cur: 0,
        })
    }
}

//---------------------------------------

struct Selection<'a> {
    structure: &'a Structure,
    state: &'a State,
    //index: It,
}

/*
impl<'a> Iterator for Selection<'a> {
    type Item = Particle<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        self.index.next().map(|id| Particle {id, atom: &self.structure.atoms[id], pos: &self.state.coords[id]})
    }
}

impl<'a> ExactSizeIterator for Selection<'a> {
    fn len(&self) -> usize {
        self.index.len()
    }
}

*/
impl<'a> Selection<'a> {
    fn from_expr(sel_str: &str) -> Result<SelectionWithAst> {
        let ast = generate_ast(sel_str)?;
        Ok(SelectionWithAst(ast))
    }

    fn from_iter<I: IndexIterator>(iter: I) -> SelectionWithIter<I> {
        SelectionWithIter(iter)
    }

    fn all() -> SelectionAll {
        SelectionAll {}
    }
}

// Helper struct for creating subscripted mutable iterator
// from iterators over atoms and positions
// IMPORTANT! Only works for **sorted** indexes!
#[derive(Clone)]
struct ParticleMutIteratorAdaptor<'a, AtomI, PosI, IndexI>
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

//##############################
//#  Tests
//##############################

#[cfg(test)]
mod tests {
    use super::Selection;
    use crate::{core::State, core::Structure, io::*};
    use lazy_static::lazy_static;

    fn read_test_pdb() -> (Structure, State) {
        let mut h = FileHandler::new_reader("tests/triclinic.pdb").unwrap();
        let structure = h.read_structure().unwrap();
        let state = h.read_next_state().unwrap().unwrap();
        (structure, state)
    }

    // Read the test PDB file once and provide the content for tests
    lazy_static! {
        static ref SS: (Structure, State) = read_test_pdb();
    }

    #[test]
    fn test_sel1() {
        let sel = Selection::from_expr("name CA").unwrap();
        let particles = sel.apply(&SS.0, &SS.1).unwrap();
        println!("sz: {}", particles.len());
        for p in particles {
            println!("{:?}", p);
        }
    }
}
