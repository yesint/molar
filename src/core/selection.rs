use std::{str::FromStr, borrow::{Borrow, BorrowMut}, cell::{RefCell, Ref}, ops::Deref, rc::Rc, sync::{RwLock, RwLockReadGuard}, marker::PhantomData};

use crate::distance_search::search::{self, SearchConnectivity, SearcherSingleGrid};

use super::{
    selection_parser::SelectionExpr,
    Atom, IndexIterator, PbcDims, PeriodicBox, Pos, State, Structure, PosIterator, structure, atom,
};
use anyhow::{anyhow, Result, bail};
use itertools::Itertools;
use ndarray::IndexLonger;
use uni_rc_lock::{UniversalRcLock, UniRcLock};

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

//-----------------------------------------------------------------

trait Select<R,S>
where
    R: UniRcLock<Structure>,
    S: UniRcLock<State>,
{
    fn select(&self,structure: R, state: S) -> Result<Selection<R,S>>;
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

impl<R,S> Select<R,S> for SelectionAll 
where
    R: for<'a> UniversalRcLock<'a,Structure>,
    S: for<'a> UniversalRcLock<'a,State>,
{    
    fn select(&self,structure: R, state: S) -> Result<Selection<R,S>> {
        let s = structure.clone();
        check_sizes(&s.read(), &state.read())?;
        let index = (0..state.read().coords.len()).collect();
        Ok(Selection { structure, state, index})
    }
}

impl<R,S> Select<R,S> for SelectionExpr 
where
    R: UniRcLock<Structure>,
    S: UniRcLock<State>,
{
    fn select(&self,structure: R, state: S) -> Result<Selection<R,S>> {
        check_sizes(&structure.read(), &state.read())?;
        let index = self.apply_whole(&structure.read(), &state.read()).unwrap(); 
        Ok(Selection {structure,state,index})
    }
}

impl<R,S> Select<R,S> for str 
where
    R: UniRcLock<Structure>,
    S: UniRcLock<State>,
{
    fn select(&self,structure: R, state: S) -> Result<Selection<R,S>> {
        check_sizes(&structure.read(), &state.read())?;
        let index = SelectionExpr::try_from(self).unwrap().apply_whole(&structure.read(), &state.read()).unwrap();
        Ok(Selection {structure,state,index})
    }
}

impl<R,S> Select<R,S> for std::ops::Range<usize>
where
    R: UniRcLock<Structure>,
    S: UniRcLock<State>,
{
    fn select(&self,structure: R, state: S) -> Result<Selection<R,S>> {
        check_sizes(&structure.read(), &state.read())?;
        let n = structure.read().atoms.len();
        if self.start>n-1 || self.end>n-1 {
            bail!("Index range {}:{} is invalid, 0:{} is allowed.",self.start,self.end,structure.read().atoms.len());
        }
        let index = self.clone().collect();
        Ok(Selection {structure,state,index})
    }
}


impl<R,S> Select<R,S> for Vec<usize>
where
    R: UniRcLock<Structure>,
    S: UniRcLock<State>,
{
    fn select(&self,structure: R, state: S) -> Result<Selection<R,S>> {
        let index = self.iter().cloned().sorted().dedup().collect();
        Ok(Selection {structure,state,index})
    }
 }


//---------------------------------------
pub struct Selection<R,S>
where
    R: UniRcLock<Structure>,
    S: UniRcLock<State>,
{
    structure: R,
    state: S,
    index: Vec<usize>,
    //ph: PhantomData<&'a R>,
}

pub struct SelectionReadGuard<'a,R,S>
where
    R: UniversalRcLock<'a,Structure>,
    S: UniversalRcLock<'a,State>,
{
    structure_ptr: <R as UniversalRcLock<'a,Structure>>::OutRead,
    state_ptr: <S as UniversalRcLock<'a,State>>::OutRead,
    index: &'a Vec<usize>,
}

pub struct SelectionWriteGuard<'a,R,S>
where
    R: UniversalRcLock<'a,Structure>,
    S: UniversalRcLock<'a,State>,
{
    structure_ptr: <R as UniversalRcLock<'a,Structure>>::OutWrite,
    state_ptr: <S as UniversalRcLock<'a,State>>::OutWrite,
    index: &'a Vec<usize>,
}


impl<R,S> Selection<R,S> 
where
    R: UniRcLock<Structure>,
    S: UniRcLock<State>,
{
    // Subselections

    pub fn subsel_from_expr(&self,expr: &SelectionExpr) -> Result<Selection<R,S>> {
        let index = expr.apply_subset(&self.structure.read(), &self.state.read(), &self.index)?;
        Ok(Selection {structure: self.structure.clone(),state: self.state.clone(),index})
    }

    pub fn subsel_from_str(&self,sel_str: &str) -> Result<Selection<R,S>> {
        let expr = SelectionExpr::try_from(sel_str)?;
        self.subsel_from_expr(&expr)
    }

    pub fn subsel_from_local_range(&self,range: std::ops::Range<usize>) -> Result<Selection<R,S>> {
        // Translate range of local indexes to global indexes
        let index: Vec<usize> = self.index.iter().cloned().skip(range.start).take(range.len()).collect();
        if index.len()==0 {
            bail!("Empty local sub-range: {}:{}, valid range: 0:{}",range.start,range.end,self.index.len()-1);
        }
        Ok(Selection {structure: self.structure.clone(),state: self.state.clone(),index})
    }

    pub fn subsel_from_iter(&self,iter: impl ExactSizeIterator<Item = usize>) -> Result<Selection<R,S>> {
        let index = iter.sorted().dedup().map(|i| self.index[i]).collect();
        Ok(Selection {structure: self.structure.clone(),state: self.state.clone(),index})
    }

    pub fn read<'a>(&'a self) -> SelectionReadGuard<'a,R,S> {
        SelectionReadGuard{
            structure_ptr: self.structure.read(),
            state_ptr: self.state.read(),
            index: &self.index,
        }
    }

    pub fn write<'a>(&'a self) -> SelectionWriteGuard<'a,R,S> {
        SelectionWriteGuard{
            structure_ptr: self.structure.write(),
            state_ptr: self.state.write(),
            index: &self.index,
        }
    }

}

impl<'a,R,S> SelectionReadGuard<'a,R,S> 
where
    R: UniversalRcLock<'a,Structure>,
    S: UniversalRcLock<'a,State>,
{
    fn iter(&self) -> impl ParticleIterator<'_> {
        ParticleIteratorAdaptor{
            atom_iter: self.structure_ptr.atoms.iter(),
            pos_iter: self.state_ptr.coords.iter(),
            index_iter: self.index.iter().cloned(),
            cur: 0,
        }
    }
}

impl<'a,R,S> SelectionWriteGuard<'a,R,S> 
where
    R: UniversalRcLock<'a,Structure>,
    S: UniversalRcLock<'a,State>,
{
    fn iter(&mut self) -> impl ParticleMutIterator<'_> {
        ParticleMutIteratorAdaptor{
            atom_iter: self.structure_ptr.atoms.iter_mut(),
            pos_iter: self.state_ptr.coords.iter_mut(),
            index_iter: self.index.iter().cloned(),
            cur: 0,
        }
    }
}

//--------------------------------------------------------
// Helper struct for creating subscripted  iterator
// from iterators over atoms and positions
// IMPORTANT! Only works for **sorted** indexes!
//--------------------------------------------------------
#[derive(Clone)]
struct ParticleIteratorAdaptor<'a, AtomI, PosI, IndexI>
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
        let sel = "name CA".select(r, s).unwrap();
        
        for p in sel.read().iter() {
            println!("{:?}",p)
        }
        
        for p in sel.write().iter() {
            println!("{:?}",p)
        }

        /*
        println!("sz: {}", particles.len());
        for p in particles {
            println!("{:?}", p);
        }
        */
    }

}
