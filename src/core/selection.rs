use std::{rc::Rc, cell::RefCell, sync::{RwLock, Arc}};
use crate::io::{IoIndexAndTopologyProvider, IoIndexAndStateProvider};
use super::{particle::*, selection_parser::SelectionExpr, State, Topology, Pos, MeasureBox, PeriodicBox, MeasureMasses, ModifyParticles, MeasurePeriodic, ModifyPeriodic, IndexIterator, ModifyRandomAccess, MeasurePos, MeasureAtoms};
use anyhow::{bail, Result};
use itertools::Itertools;
use uni_rc_lock::UniRcLock;

//-----------------------------------------------------------------

/// Trait whic provides select method operating with generic smart pointers
pub trait Select<T, S>
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

pub type SelectionRc = Selection< Rc<RefCell<Topology>>, Rc<RefCell<State>> >;
pub type SelectionArc = Selection< Arc<RwLock<Topology>>, Arc<RwLock<State>> >;

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

    pub fn query<'a>(&'a self) -> SelectionQueryGuard<'a, T, S> {
        SelectionQueryGuard {
            topology_ref: self.topology.read(),
            state_ref: self.state.read(),
            index: &self.index,
        }
    }

    pub fn modify<'a>(&'a self) -> SelectionModifyGuard<'a, T, S> {
        SelectionModifyGuard {
            topology_ref: self.topology.write(),
            state_ref: self.state.write(),
            index: &self.index,
        }
    }


}
//---------------,-------------------------------------

/// Scoped guard giving read-only access to selection
pub struct SelectionQueryGuard<'a, T, S>
where
    T: UniRcLock<Topology> + 'a,
    S: UniRcLock<State> + 'a,
{
    topology_ref: <T as UniRcLock<Topology>>::OutRead<'a>,
    state_ref: <S as UniRcLock<State>>::OutRead<'a>,
    index: &'a Vec<usize>,
}

impl<'a,T,S> IoIndexAndTopologyProvider for SelectionQueryGuard<'a, T, S>
where
    T: UniRcLock<Topology> + 'a,
    S: UniRcLock<State> + 'a,
{
    fn get_index_and_topology(&self) -> (impl IndexIterator, &Topology) {
        (self.index.iter().cloned(), &self.topology_ref)
    }
}

impl<'a,T,S> IoIndexAndStateProvider for SelectionQueryGuard<'a, T, S>
where
    T: UniRcLock<Topology> + 'a,
    S: UniRcLock<State> + 'a,
{
    fn get_index_and_state(&self) -> (impl IndexIterator, &State) {
        (self.index.iter().cloned(), &self.state_ref)
    }
}

/// Scoped guard giving read-write access to selection
pub struct SelectionModifyGuard<'a, T, S>
where
    T: UniRcLock<Topology> + 'a,
    S: UniRcLock<State> + 'a,
{
    topology_ref: <T as UniRcLock<Topology>>::OutWrite<'a>,
    state_ref: <S as UniRcLock<State>>::OutWrite<'a>,
    index: &'a Vec<usize>,
}

//==================================================================
// Implement analysis traits

impl<T,S> MeasureBox for SelectionQueryGuard<'_,T,S> 
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{
    fn get_box(&self) -> Result<&PeriodicBox> {
        self.state_ref.get_box()
    }
}

impl<T, S> MeasurePos for SelectionQueryGuard<'_, T, S>
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{
    fn iter_pos(&self) -> impl super::PosIterator<'_> {
        self.index.iter().map(|i| &self.state_ref.coords[*i])
    }
}

impl<T, S> MeasureAtoms for SelectionQueryGuard<'_, T, S>
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{
    fn iter_atoms(&self) -> impl super::AtomIterator<'_> {
        self.index.iter().map(|i| &self.topology_ref.atoms[*i])
    }    
}


impl<T, S> MeasureMasses for SelectionQueryGuard<'_, T, S>
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{
    /*
    fn iter_particles(&self) -> impl ParticleIterator<'_> {
        ParticleIteratorAdaptor::new(
            &self.topology_ref.atoms,
            &self.state_ref.coords,
            &self.index,
        )
    }
    */

    fn iter_masses(&self) -> impl ExactSizeIterator<Item = f32> {
        self.index.iter().map(|i| self.topology_ref.atoms[*i].mass)
    }
}

impl<T, S> MeasurePeriodic for SelectionQueryGuard<'_, T, S>
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{}

impl<T,S> MeasureBox for SelectionModifyGuard<'_,T,S> 
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{
    fn get_box(&self) -> Result<&PeriodicBox> {
        self.state_ref.get_box()
    }
}

impl<T, S> ModifyParticles for SelectionModifyGuard<'_, T, S>
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{
    fn iter_particles_mut(&mut self) -> impl ParticleMutIterator<'_> {
        ParticleMutIteratorAdaptor::new(
            &mut self.topology_ref.atoms,
            &mut self.state_ref.coords,
            &self.index,
        )
    }
}

impl<T, S> ModifyPeriodic for SelectionModifyGuard<'_, T, S>
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{}

impl<T, S> ModifyRandomAccess for SelectionModifyGuard<'_, T, S>
where
    T: UniRcLock<Topology>,
    S: UniRcLock<State>,
{
    fn nth_particle_mut(&mut self, i: usize) -> ParticleMut {
        ParticleMut{
            id: i,
            pos: &mut self.state_ref.coords[self.index[i]],
            atom: &mut self.topology_ref.atoms[self.index[i]],
        }
    }

    fn nth_pos_mut(&mut self, i: usize) -> &mut Pos {
        &mut self.state_ref.coords[self.index[i]]
    }
}

//############################################################
//#  Tests
//############################################################

#[cfg(test)]
mod tests {
    use crate::{
        core::State,
        core::{selection::Select, Topology, MeasureMasses, MeasurePos, ModifyParticles, Vector3f, ModifyRandomAccess, fit_transform},
        io::*,
    };
    use lazy_static::lazy_static;
    use nalgebra::Unit;

    use super::{SelectionRc, SelectionAll};

    fn read_test_pdb() -> (Topology, State) {
        let mut h = FileHandler::new_reader("tests/no_ATP.pdb").unwrap();
        let top = h.read_topology().unwrap();
        let state = h.read_next_state().unwrap().unwrap();
        (top, state)
    }

    // Read the test PDB file once and provide the content for tests
    lazy_static! {
        static ref SS: (Topology, State) = read_test_pdb();
    }

    fn make_sel() -> anyhow::Result<SelectionRc> {
        let t = SS.0.clone().to_rc();
        let s = SS.1.clone().to_rc();
        let sel = SelectionAll{}.select(t, s)?;
        Ok(sel)
    }

    fn make_sel_prot() -> anyhow::Result<SelectionRc> {
        let t = SS.0.clone().to_rc();
        let s = SS.1.clone().to_rc();
        let sel = "not resname TIP3 POT CLA".select(t, s)?;
        Ok(sel)
    }

    #[test]
    fn test_make_sel() ->anyhow::Result<()>{
        make_sel().map(|_| ())
    }

    #[test]
    fn test_measure() -> anyhow::Result<()>{
        let sel = make_sel()?;

        println!("before {}",sel.query().iter_pos().next().unwrap());

        let (minv,maxv) = sel.query().min_max();
        println!("{minv}:{maxv}");

        //sel.translate(&Vector3f::new(10.0,10.0,10.0));
        println!("after {}",sel.query().iter_pos().next().unwrap());
        println!("{:?}",sel.query().min_max());
        Ok(())
    }
    
    #[test]
    fn test_measure_pbc() -> anyhow::Result<()>{
        let sel = make_sel()?;

        let cm = sel.query().center_of_mass()?;
        println!("{cm}");
        Ok(())
    }

    #[test]
    fn test_translate() {
        let sel = make_sel().unwrap();
        let mut w = sel.modify();

        println!("before {}",w.iter_pos_mut().next().unwrap());
        w.translate(Vector3f::new(10.0,10.0,10.0));
        println!("after {}",w.iter_pos_mut().next().unwrap());
    }

    #[test]
    fn test_write_to_file() -> anyhow::Result<()> {
        let sel = make_sel()?;
        let q = sel.query();

        let mut h = FileHandler::new_writer("f.pdb")?;
        h.write_topology(&q)?;
        h.write_next_state(&q)?;
        Ok(())
    }

    #[test]
    fn test_unwrap_connectivity_1() -> anyhow::Result<()> {
        let sel = make_sel_prot()?;
        sel.modify().unwrap_connectivity_dim(0.2, &[true,true,true])?;
        
        let mut h = FileHandler::new_writer("unwrapped.pdb")?;
        let q = sel.query();
        h.write_topology(&q)?;
        h.write_next_state(&q)?;
        Ok(())
    }

    #[test]
    fn eigen_test() -> anyhow::Result<()> {
        let sel1 = make_sel_prot().unwrap();
        let sel2 = make_sel_prot().unwrap();
        //let cm1 = sel1.query().center_of_mass()?.coords;
        //let cm2 = sel2.query().center_of_mass()?.coords;
        //sel1.modify().translate(-cm1);
        //sel2.modify().translate(-cm2);
        
        sel2.modify().rotate(&Unit::new_normalize(Vector3f::x()), 80.0_f32.to_radians());
        
        let mut h = FileHandler::new_writer("sel2.pdb")?;
        let q = sel2.query();
        h.write_topology(&q)?;
        h.write_next_state(&q)?;
        drop(q);

        let mut h = FileHandler::new_writer("sel1_before.pdb")?;
        let q = sel1.query();
        h.write_topology(&q)?;
        h.write_next_state(&q)?;
        drop(q);


        let m = fit_transform(sel1.query(), sel2.query())?;
        println!("{m}");
        sel1.modify().apply_transform(&m);


        let mut h = FileHandler::new_writer("sel1_after.pdb")?;
        let q = sel1.query();
        h.write_topology(&q)?;
        h.write_next_state(&q)?;

        Ok(())
    }
}
