use std::str::FromStr;

use dyn_clone::clone_box;


use crate::distance_search::search::{self, SearchConnectivity, SearcherSingleGrid};

use super::{
    selection_parser::SelectionExpr,
    Atom, IndexIterator, PbcDims, PeriodicBox, Pos, State, Structure, PosIterator,
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

pub trait ParticleIterator<'a>: ExactSizeIterator<Item = Particle<'a>> {}
impl<'a, T> ParticleIterator<'a> for T where T: ExactSizeIterator<Item = Particle<'a>> {}

pub trait ParticleMutIterator<'a>: ExactSizeIterator<Item = ParticleMut<'a>> {}
impl<'a, T> ParticleMutIterator<'a> for T where T: ExactSizeIterator<Item = ParticleMut<'a>> {}

//-----------------------------------------------------------------

trait Select {
    fn select<'a>(&self,structure: &'a Structure, state: &'a State) -> Selection<'a>;
    fn select_mut<'a>(&self,structure: &'a mut Structure, state: &'a mut State) -> SelectionMut<'a>;
}

//-----------------------------------------------------------------

struct SelectionAll {}

impl Select for SelectionAll {    
    fn select<'a>(&self,structure: &'a Structure, state: &'a State) -> Selection<'a> {
        Selection { structure, state, index: (0..state.coords.len()).collect()}
    }

    fn select_mut<'a>(&self,structure: &'a mut Structure, state: &'a mut State) -> SelectionMut<'a> {
        let index = (0..state.coords.len()).collect();
        SelectionMut { structure, state, index}
    }
}

impl Select for SelectionExpr {
    fn select<'a>(&self,structure: &'a Structure, state: &'a State) -> Selection<'a> {
        let index = self.apply_whole(structure, state).unwrap(); 
        Selection {structure,state,index}
    }

    fn select_mut<'a>(&self,structure: &'a mut Structure, state: &'a mut State) -> SelectionMut<'a> {
        let index = self.apply_whole(structure, state).unwrap();
        SelectionMut {structure,state,index}
    }
}

impl Select for str {
    fn select<'a>(&self,structure: &'a Structure, state: &'a State) -> Selection<'a> {
        let index = SelectionExpr::try_from(self).unwrap().apply_whole(structure, state).unwrap();
        Selection {structure,state,index}
    }

    fn select_mut<'a>(&self,structure: &'a mut Structure, state: &'a mut State) -> SelectionMut<'a> {
        let index = SelectionExpr::try_from(self).unwrap().apply_whole(structure, state).unwrap();
        SelectionMut {structure,state,index}
    }
}

impl Select for std::ops::Range<usize> {
    fn select<'a>(&self,structure: &'a Structure, state: &'a State) -> Selection<'a> {
        let index = self.clone().collect();
        Selection {structure,state,index}
    }

    fn select_mut<'a>(&self,structure: &'a mut Structure, state: &'a mut State) -> SelectionMut<'a> {
        let index = self.clone().collect();
        SelectionMut {structure,state,index}
    }
}

/*
impl<T> Select for T 
where
    T: Iterator<Item=usize>+Clone
{
    fn select<'a>(&self,structure: &'a Structure, state: &'a State) -> Selection<'a> {
        let index = self.clone().collect();
        Selection {structure,state,index}
    }

    fn select_mut<'a>(&self,structure: &'a mut Structure, state: &'a mut State) -> SelectionMut<'a> {
        let index = self.clone().collect();
        SelectionMut {structure,state,index}
    }
}
*/



//---------------------------------------
struct Selection<'a> {
    structure: &'a Structure,
    state: &'a State,
    index: Vec<usize>,
}

struct SelectionMut<'a> {
    structure: &'a mut Structure,
    state: &'a mut State,
    index: Vec<usize>,
}

 

impl<'a> SelectionMut<'a> {
    fn iter_particle(&'a mut self) -> Result<impl ParticleMutIterator<'a>+'_> {        
        Ok(ParticleMutIteratorAdaptor {
            atom_iter: self.structure.atoms.iter_mut(),
            pos_iter: self.state.coords.iter_mut(),
            index_iter: self.index.iter().cloned(),
            cur: 0,
        })
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
    use std::{borrow::Borrow, fmt::Debug};

    use super::Selection;
    use crate::{
        core::State,
        core::{Pos, PosIterator, Structure, selection::Select, StructureHandle, StateHandle},
        io::*,
    };
    use lazy_static::lazy_static;

    fn read_test_pdb() -> (StructureHandle, StateHandle) {
        let mut h = FileHandler::new_reader("tests/triclinic.pdb").unwrap();
        let structure = h.read_structure().unwrap();
        let state = h.read_next_state().unwrap().unwrap();
        (structure, state)
    }

    // Read the test PDB file once and provide the content for tests
    lazy_static! {
        static ref SS: (StructureHandle, StateHandle) = read_test_pdb();
    }

    #[test]
    fn test_sel1() {
        let sel = "name CA".select(&SS.0.read(), &SS.1.read());
        /*
        let particles = sel.apply(&SS.0, &SS.1).unwrap();
        println!("sz: {}", particles.len());
        for p in particles {
            println!("{:?}", p);
        }
        */
    }

}
