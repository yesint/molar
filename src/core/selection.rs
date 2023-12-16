use super::{particle::*, selection_parser::SelectionExpr, Atom, PosIterator, State, Structure, Pos, Vector3f};
use anyhow::{bail, Result};
use itertools::Itertools;
use num_traits::Bounded;
use uni_rc_lock::{UniRcLock, UniversalRcLock};

//-----------------------------------------------------------------

/// Trait whic provides select method operating with generic smart pointers
trait Select<R, S>
where
    R: UniRcLock<Structure>,
    S: UniRcLock<State>,
{
    fn select(&self, structure: R, state: S) -> Result<Selection<R, S>>;
}

//-----------------------------------------------------------------
fn check_sizes(structure: &Structure, state: &State) -> Result<()> {
    let n1 = structure.atoms.len();
    let n2 = state.coords.len();
    match n1 == n2 {
        true => Ok(()),
        false => bail!(
            "Structure and State are incompatible (sizes {} and {})",
            n1,
            n2
        ),
    }
}

struct SelectionAll {}

impl<R, S> Select<R, S> for SelectionAll
where
    R: for<'a> UniversalRcLock<'a, Structure>,
    S: for<'a> UniversalRcLock<'a, State>,
{
    fn select(&self, structure: R, state: S) -> Result<Selection<R, S>> {
        let s = structure.clone();
        check_sizes(&s.read(), &state.read())?;
        let index: Vec<usize> = (0..state.read().coords.len()).collect();
        if index.len() >0 {
            Ok(Selection {
                structure,
                state,
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }
}

impl<R, S> Select<R, S> for SelectionExpr
where
    R: UniRcLock<Structure>,
    S: UniRcLock<State>,
{
    fn select(&self, structure: R, state: S) -> Result<Selection<R, S>> {
        check_sizes(&structure.read(), &state.read())?;
        let index = self.apply_whole(&structure.read(), &state.read()).unwrap();
        if index.len() >0 {
            Ok(Selection {
                structure,
                state,
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }
}

impl<R, S> Select<R, S> for str
where
    R: UniRcLock<Structure>,
    S: UniRcLock<State>,
{
    fn select(&self, structure: R, state: S) -> Result<Selection<R, S>> {
        check_sizes(&structure.read(), &state.read())?;
        let index = SelectionExpr::try_from(self)
            .unwrap()
            .apply_whole(&structure.read(), &state.read())
            .unwrap();
        if index.len() >0 {
            Ok(Selection {
                structure,
                state,
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }
}

impl<R, S> Select<R, S> for std::ops::Range<usize>
where
    R: UniRcLock<Structure>,
    S: UniRcLock<State>,
{
    fn select(&self, structure: R, state: S) -> Result<Selection<R, S>> {
        check_sizes(&structure.read(), &state.read())?;
        let n = structure.read().atoms.len();
        if self.start > n - 1 || self.end > n - 1 {
            bail!(
                "Index range {}:{} is invalid, 0:{} is allowed.",
                self.start,
                self.end,
                structure.read().atoms.len()
            );
        }
        let index: Vec<usize> = self.clone().collect();
        if index.len() >0 {
            Ok(Selection {
                structure,
                state,
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }
}

impl<R, S> Select<R, S> for Vec<usize>
where
    R: UniRcLock<Structure>,
    S: UniRcLock<State>,
{
    fn select(&self, structure: R, state: S) -> Result<Selection<R, S>> {
        let index: Vec<usize> = self.iter().cloned().sorted().dedup().collect();
        if index.len() >0 {
            Ok(Selection {
                structure,
                state,
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }
}

//---------------------------------------
pub struct Selection<R, S>
where
    R: UniRcLock<Structure>,
    S: UniRcLock<State>,
{
    structure: R,
    state: S,
    index: Vec<usize>,
}

impl<R, S> Selection<R, S>
where
    R: UniRcLock<Structure>,
    S: UniRcLock<State>,
{
    /// Subselection from expression
    pub fn subsel_from_expr(&self, expr: &SelectionExpr) -> Result<Selection<R, S>> {
        let index = expr.apply_subset(&self.structure.read(), &self.state.read(), &self.index)?;
        if index.len() >0 {
            Ok(Selection {
                structure: self.structure.clone(),
                state: self.state.clone(),
                index,
            })    
        } else {
            bail!("Selection is empty")
        }
    }

    /// Subselection from string
    pub fn subsel_from_str(&self, sel_str: &str) -> Result<Selection<R, S>> {
        let expr = SelectionExpr::try_from(sel_str)?;
        self.subsel_from_expr(&expr)
    }

    /// Subselection from the range of local selection indexes
    pub fn subsel_from_local_range(
        &self,
        range: std::ops::Range<usize>,
    ) -> Result<Selection<R, S>> {
        // Translate range of local indexes to global indexes
        let index: Vec<usize> = self
            .index
            .iter()
            .cloned()
            .skip(range.start)
            .take(range.len())
            .collect();
        if index.len() == 0 {
            bail!(
                "Empty local sub-range: {}:{}, valid range: 0:{}",
                range.start,
                range.end,
                self.index.len() - 1
            );
        }

        if index.len() >0 {
            Ok(Selection {
                structure: self.structure.clone(),
                state: self.state.clone(),
                index,
            })    
        } else {
            bail!("Selection is empty")
        }
    }

    /// Subselection from iterator that provides local selection indexes
    pub fn subsel_from_iter(
        &self,
        iter: impl ExactSizeIterator<Item = usize>,
    ) -> Result<Selection<R, S>> {
        let index: Vec<usize> = iter.sorted().dedup().map(|i| self.index[i]).collect();
        if index.len() >0 {
            Ok(Selection {
                structure: self.structure.clone(),
                state: self.state.clone(),
                index,
            })    
        } else {
            bail!("Selection is empty")
        }
    }

    pub fn read<'a>(&'a self) -> SelectionReadGuard<'a, R, S> {
        SelectionReadGuard {
            structure_ref: self.structure.read(),
            state_ref: self.state.read(),
            index: &self.index,
        }
    }

    pub fn write<'a>(&'a self) -> SelectionWriteGuard<'a, R, S> {
        SelectionWriteGuard {
            structure_ref: self.structure.write(),
            state_ref: self.state.write(),
            index: &self.index,
        }
    }

}
//----------------------------------------------------

/// Scoped guard giving read-only access to selection
pub struct SelectionReadGuard<'a, R, S>
where
    R: UniversalRcLock<'a, Structure>,
    S: UniversalRcLock<'a, State>,
{
    structure_ref: <R as UniversalRcLock<'a, Structure>>::OutRead,
    state_ref: <S as UniversalRcLock<'a, State>>::OutRead,
    index: &'a Vec<usize>,
}

impl<'a, R, S> SelectionReadGuard<'a, R, S>
where
    R: UniversalRcLock<'a, Structure>,
    S: UniversalRcLock<'a, State>,
{
    fn iter(&self) -> impl ParticleIterator<'_> {
        ParticleIteratorAdaptor::new(
            &self.structure_ref.atoms,
            &self.state_ref.coords,
            &self.index,
        )
    }

    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.iter().map(|p| p.pos)
    }

    fn iter_atoms(&self) -> impl ExactSizeIterator<Item = &'_ Atom> {
        self.iter().map(|p| p.atom)
    }

    pub fn min_max(&self) -> (Pos,Pos) {
        let mut lower = Pos::max_value();
        let mut upper = Pos::min_value();
        for p in self.iter_pos() {
            for d in 0..3 {
                if p[d] < lower[d] { lower[d] = p[d] }
                if p[d] > upper[d] { upper[d] = p[d] }
            }
        }
        (lower,upper)
    }

    
}

/// Scoped guard giving read-write access to selection
pub struct SelectionWriteGuard<'a, R, S>
where
    R: UniversalRcLock<'a, Structure>,
    S: UniversalRcLock<'a, State>,
{
    structure_ref: <R as UniversalRcLock<'a, Structure>>::OutWrite,
    state_ref: <S as UniversalRcLock<'a, State>>::OutWrite,
    index: &'a Vec<usize>,
}

impl<'a, R, S> SelectionWriteGuard<'a, R, S>
where
    R: UniversalRcLock<'a, Structure>,
    S: UniversalRcLock<'a, State>,
{
    fn iter(&mut self) -> impl ParticleMutIterator<'_> {
        ParticleMutIteratorAdaptor::new(
            self.structure_ref.atoms.iter_mut(),
            self.state_ref.coords.iter_mut(),
            self.index.iter().cloned(),
        )
    }
}

//==================================================================

//##############################
//#  Tests
//##############################

#[cfg(test)]
mod tests {

    use crate::{
        core::State,
        core::{selection::Select, Structure, Vector3f},
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

        //for p in sel.read().iter().unwrap() {
        //    println!("{:?}", p)
        //}

        //for p in sel.write().iter().unwrap() {
        //    println!("{:?}", p)
        //}

        println!("before {}",sel.read().iter_pos().next().unwrap());

        let (minv,maxv) = sel.read().min_max();
        println!("{minv}:{maxv}");

        //sel.translate(&Vector3f::new(10.0,10.0,10.0));
        println!("after {}",sel.read().iter_pos().next().unwrap());
        println!("{:?}",sel.read().min_max());

        /*
        println!("sz: {}", particles.len());
        for p in particles {
            println!("{:?}", p);
        }
        */
    }
}
