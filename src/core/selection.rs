use super::{particle::*, selection_parser::SelectionExpr, Atom, PosIterator, State, Topology, Pos, BoxProvider, PeriodicBox};
use anyhow::{bail, Result, anyhow};
use itertools::Itertools;
use num_traits::Bounded;
use uni_rc_lock::UniRcLock;

//-----------------------------------------------------------------

/// Trait whic provides select method operating with generic smart pointers
trait Select<T, S>
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{
    fn select(&self, topology: T, state: S) -> Result<Selection<T, S>>;
}

//-----------------------------------------------------------------
fn check_sizes(topology: &Topology, state: &State) -> Result<()> {
    let n1 = topology.atoms.len();
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

impl<T, S> Select<T, S> for SelectionAll
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{
    fn select(&self, topology: T, state: S) -> Result<Selection<T, S>> {
        let s = topology.clone();
        check_sizes(&s.read(), &state.read())?;
        let index: Vec<usize> = (0..state.read().coords.len()).collect();
        if index.len() >0 {
            Ok(Selection {
                topology,
                state,
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }
}

impl<T, S> Select<T, S> for SelectionExpr
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{
    fn select(&self, topology: T, state: S) -> Result<Selection<T, S>> {
        check_sizes(&topology.read(), &state.read())?;
        let index = self.apply_whole(&topology.read(), &state.read()).unwrap();
        if index.len() >0 {
            Ok(Selection {
                topology,
                state,
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }
}

impl<T, S> Select<T, S> for str
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{
    fn select(&self, topology: T, state: S) -> Result<Selection<T, S>> {
        check_sizes(&topology.read(), &state.read())?;
        let index = SelectionExpr::try_from(self)
            .unwrap()
            .apply_whole(&topology.read(), &state.read())
            .unwrap();
        if index.len() >0 {
            Ok(Selection {
                topology: topology,
                state,
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }
}

impl<T, S> Select<T, S> for std::ops::Range<usize>
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{
    fn select(&self, topology: T, state: S) -> Result<Selection<T, S>> {
        check_sizes(&topology.read(), &state.read())?;
        let n = topology.read().atoms.len();
        if self.start > n - 1 || self.end > n - 1 {
            bail!(
                "Index range {}:{} is invalid, 0:{} is allowed.",
                self.start,
                self.end,
                topology.read().atoms.len()
            );
        }
        let index: Vec<usize> = self.clone().collect();
        if index.len() >0 {
            Ok(Selection {
                topology,
                state,
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }
}

impl<T, S> Select<T, S> for Vec<usize>
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{
    fn select(&self, topology: T, state: S) -> Result<Selection<T, S>> {
        let index: Vec<usize> = self.iter().cloned().sorted().dedup().collect();
        if index.len() >0 {
            Ok(Selection {
                topology,
                state,
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }
}

//---------------------------------------
pub struct Selection<T, S>
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{
    topology: T,
    state: S,
    index: Vec<usize>,
}

impl<T, S> Selection<T, S>
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{
    /// Subselection from expression
    pub fn subsel_from_expr(&self, expr: &SelectionExpr) -> Result<Selection<T, S>> {
        let index = expr.apply_subset(&self.topology.read(), &self.state.read(), &self.index)?;
        if index.len() >0 {
            Ok(Selection {
                topology: self.topology.clone(),
                state: self.state.clone(),
                index,
            })    
        } else {
            bail!("Selection is empty")
        }
    }

    /// Subselection from string
    pub fn subsel_from_str(&self, sel_str: &str) -> Result<Selection<T, S>> {
        let expr = SelectionExpr::try_from(sel_str)?;
        self.subsel_from_expr(&expr)
    }

    /// Subselection from the range of local selection indexes
    pub fn subsel_from_local_range(
        &self,
        range: std::ops::Range<usize>,
    ) -> Result<Selection<T, S>> {
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
                topology: self.topology.clone(),
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
    ) -> Result<Selection<T, S>> {
        let index: Vec<usize> = iter.sorted().dedup().map(|i| self.index[i]).collect();
        if index.len() >0 {
            Ok(Selection {
                topology: self.topology.clone(),
                state: self.state.clone(),
                index,
            })    
        } else {
            bail!("Selection is empty")
        }
    }

    pub fn read<'a>(&'a self) -> SelectionReadGuard<'a, T, S> {
        SelectionReadGuard {
            topology_ref: self.topology.read(),
            state_ref: self.state.read(),
            index: &self.index,
        }
    }

    pub fn write<'a>(&'a self) -> SelectionWriteGuard<'a, T, S> {
        SelectionWriteGuard {
            topology_ref: self.topology.write(),
            state_ref: self.state.write(),
            index: &self.index,
        }
    }

}
//----------------------------------------------------

/// Scoped guard giving read-only access to selection
pub struct SelectionReadGuard<'a, T, S>
where
    T: UniRcLock<Topology> + 'a,
    S: UniRcLock<State> + 'a,
{
    topology_ref: <T as UniRcLock<Topology>>::OutRead<'a>,
    state_ref: <S as UniRcLock<State>>::OutRead<'a>,
    index: &'a Vec<usize>,
}

impl<T, S> SelectionReadGuard<'_, T, S>
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{
    fn iter(&self) -> impl ParticleIterator<'_> {
        ParticleIteratorAdaptor::new(
            &self.topology_ref.atoms,
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
pub struct SelectionWriteGuard<'a, T, S>
where
    T: UniRcLock<Topology> + 'a,
    S: UniRcLock<State> + 'a,
{
    topology_ref: <T as UniRcLock<Topology>>::OutWrite<'a>,
    state_ref: <S as UniRcLock<State>>::OutWrite<'a>,
    index: &'a Vec<usize>,
}

impl<T, S> SelectionWriteGuard<'_, T, S>
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{
    fn iter(&mut self) -> impl ParticleMutIterator<'_> {
        ParticleMutIteratorAdaptor::new(
            self.topology_ref.atoms.iter_mut(),
            self.state_ref.coords.iter_mut(),
            self.index.iter().cloned(),
        )
    }
}

//==================================================================
// Implement analysis traits

impl<T,S> BoxProvider for SelectionReadGuard<'_,T,S> 
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{
    fn get_box(&self) -> Result<&PeriodicBox> {
        let r = self.state_ref.box_
            .as_ref()
            .ok_or(anyhow!("No periodic box"))?;
        Ok(&r)
    }
}

impl<T,S> BoxProvider for SelectionWriteGuard<'_,T,S> 
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{
    fn get_box(&self) -> Result<&PeriodicBox> {
        let r = self.state_ref.box_
            .as_ref()
            .ok_or(anyhow!("No periodic box"))?;
        Ok(&r)
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
        core::{selection::Select, Topology},
        io::*,
    };
    use lazy_static::lazy_static;

    fn read_test_pdb() -> (Topology, State) {
        let mut h = FileHandler::new_reader("tests/triclinic.pdb").unwrap();
        let top = h.read_topology().unwrap();
        let state = h.read_next_state().unwrap().unwrap();
        (top, state)
    }

    // Read the test PDB file once and provide the content for tests
    lazy_static! {
        static ref SS: (Topology, State) = read_test_pdb();
    }

    #[test]
    fn test_sel1() {
        let t = SS.0.clone().to_rc();
        let s = SS.1.clone().to_rc();
        let sel = "name CA".select(t, s).unwrap();

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
