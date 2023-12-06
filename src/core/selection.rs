use std::{str::FromStr, borrow::{Borrow, BorrowMut}, cell::{RefCell, Ref}, ops::Deref, rc::Rc, sync::{RwLock, RwLockReadGuard}};

use crate::distance_search::search::{self, SearchConnectivity, SearcherSingleGrid};

use super::{
    selection_parser::SelectionExpr,
    Atom, IndexIterator, PbcDims, PeriodicBox, Pos, State, Structure, PosIterator, structure,
};
use anyhow::{anyhow, Result, bail};
use itertools::Itertools;
use uni_rc_lock::*;

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
    fn select<R: UniRcLock<Structure>,S: UniRcLock<State>>(&self,structure: R, state: S) -> Result<Selection<R,S>>;
}

//-----------------------------------------------------------------
fn check_sizes(structure: &Structure, state: &State) -> Result<()> {
    let n1 = structure.atoms.len();
    let n2 = state.coords.len();
    match n1==n2 {
        true => Ok(()),
        false => bail!("Structure and State are incompatible (sizes {} and {})",n1,n2),
    }
}


struct SelectionAll {}

impl Select for SelectionAll {    
    fn select<R: UniRcLock<Structure>,S: UniRcLock<State>>(&self,structure: R, state: S) -> Result<Selection<R,S>> {
        check_sizes(&structure.read(), &state.read())?;
        let index = (0..state.read().coords.len()).collect();
        Ok(Selection { structure, state, index})
    }
}

impl Select for SelectionExpr {
    fn select<R: UniRcLock<Structure>,S: UniRcLock<State>>(&self,structure: R, state: S) -> Result<Selection<R,S>> {
        check_sizes(&structure.read(), &state.read())?;
        let index = self.apply_whole(&structure.read(), &state.read()).unwrap(); 
        Ok(Selection {structure,state,index})
    }
}

impl Select for str {
    fn select<R: UniRcLock<Structure>,S: UniRcLock<State>>(&self,structure: R, state: S) -> Result<Selection<R,S>> {
        check_sizes(&structure.read(), &state.read())?;
        let index = SelectionExpr::try_from(self).unwrap().apply_whole(&structure.read(), &state.read()).unwrap();
        Ok(Selection {structure,state,index})
    }
}

impl Select for std::ops::Range<usize> {
    fn select<R: UniRcLock<Structure>,S: UniRcLock<State>>(&self,structure: R, state: S) -> Result<Selection<R,S>> {
        check_sizes(&structure.read(), &state.read())?;
        let n = structure.read().atoms.len();
        if self.start>n-1 || self.end>n-1 {
            bail!("Index range {}:{} is invalid, 0:{} is allowed.",self.start,self.end,structure.read().atoms.len());
        }
        let index = self.clone().collect();
        Ok(Selection {structure,state,index})
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
struct Selection<R,S>
where
    R: UniRcLock<Structure>,
    S: UniRcLock<State>,
{
    structure: R,
    state: S,
    index: Vec<usize>,
}

impl<R,S> Selection<R,S> 
where
    R: UniRcLock<Structure>,
    S: UniRcLock<State>,
{
    fn subsel_from_expr(&self,expr: &SelectionExpr) -> Result<Selection<R,S>> {
        let index = expr.apply_subset(&self.structure.read(), &self.state.read(), &self.index)?;
        Ok(Selection {structure: self.structure.clone(),state: self.state.clone(),index})
    }

    fn subsel_from_str(&self,sel_str: &str) -> Result<Selection<R,S>> {
        let expr = SelectionExpr::try_from(sel_str)?;
        self.subsel_from_expr(&expr)
    }

    fn subsel_from_local_range(&self,range: std::ops::Range<usize>) -> Result<Selection<R,S>> {
        // Translate range of local indexes to global indexes
        let index: Vec<usize> = self.index.iter().cloned().skip(range.start).take(range.len()).collect();
        if index.len()==0 {
            bail!("Empty local sub-range: {}:{}, valid range: 0:{}",range.start,range.end,self.index.len()-1);
        }
        Ok(Selection {structure: self.structure.clone(),state: self.state.clone(),index})
    }

    fn subsel_from_iter(&self,iter: impl ExactSizeIterator<Item = usize>) -> Result<Selection<R,S>> {
        let index = iter.sorted().dedup().map(|i| self.index[i]).collect();
        Ok(Selection {structure: self.structure.clone(),state: self.state.clone(),index})
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
    use std::{borrow::Borrow, fmt::Debug, cell::RefCell, rc::Rc};

    use super::{Selection};
    use crate::{
        core::State,
        core::{Pos, PosIterator, Structure, selection::Select, structure},
        io::*,
    };
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
        let r = SS.0.clone().to_rc();
        let s = SS.1.clone().to_rc();
        let sel = "name CA".select(r, s);
        
        let particles = sel;
        /*
        println!("sz: {}", particles.len());
        for p in particles {
            println!("{:?}", p);
        }
        */
    }

}
