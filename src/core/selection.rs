use std::{cell::{Ref, RefCell, RefMut}, rc::Rc, sync::{Arc, RwLock}};
use crate::io::{IoIndexProvider, IoTopologyProvider, IoStateProvider};
use super::{IndexIterator, MeasureAtoms, MeasureBox, MeasureMasses, MeasurePeriodic, MeasurePos, ModifyPeriodic, ModifyPos, ModifyRandomAccess, PeriodicBox, Pos, PosIterator, PosMutIterator, State, StateRc, Topology, TopologyRc};
use anyhow::{bail, Result};
use itertools::Itertools;
use uni_rc_lock::UniRcLock;

pub use super::selection_parser::SelectionExpr;
//-----------------------------------------------------------------

/// Trait whic provides select method operating with generic smart pointers
pub trait Select {
    fn select(&self, topology: &TopologyRc, state: &StateRc) -> Result<Selection>;
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

pub struct SelectionAll {}

impl SelectionAll {
    pub fn new() -> Self {
        Self{}
    }
}

impl Select for SelectionAll {
    fn select(&self, topology: &TopologyRc, state: &StateRc) -> Result<Selection> {
        check_sizes(&topology.read(), &state.read())?;
        let index: Vec<usize> = (0..state.read().coords.len()).collect();
        if index.len() >0 {
            Ok(Selection {
                topology: topology.clone(),
                state: state.clone(),
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }
}

impl Select for SelectionExpr {
    fn select(&self, topology: &TopologyRc, state: &StateRc) -> Result<Selection> {
        check_sizes(&topology.read(), &state.read())?;
        let index = self.apply_whole(&topology.read(), &state.read()).unwrap();
        if index.len() >0 {
            Ok(Selection {
                topology: topology.clone(),
                state: state.clone(),
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }
}

impl Select for str {
    fn select(&self, topology: &TopologyRc, state: &StateRc) -> Result<Selection> {
        check_sizes(&topology.read(), &state.read())?;
        let index = SelectionExpr::try_from(self)
            .unwrap()
            .apply_whole(&topology.read(), &state.read())
            .unwrap();
        if index.len() >0 {
            Ok(Selection {
                topology: topology.clone(),
                state: state.clone(),
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }
}

impl Select for std::ops::Range<usize> {
    fn select(&self, topology: &TopologyRc, state: &StateRc) -> Result<Selection> {
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
                topology: topology.clone(),
                state: state.clone(),
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }
}

impl Select for Vec<usize> {
    fn select(&self, topology: &TopologyRc, state: &StateRc) -> Result<Selection> {
        let index: Vec<usize> = self.iter().cloned().sorted().dedup().collect();
        if index.len() >0 {
            Ok(Selection {
                topology: topology.clone(),
                state: state.clone(),
                index,
            })
        } else {
            bail!("Selection is empty")
        }
    }
}

//---------------------------------------
pub struct Selection {
    topology: TopologyRc,
    state: StateRc,
    index: Vec<usize>,
}

impl Selection {
    /// Subselection from expression
    pub fn subsel_from_expr(&self, expr: &SelectionExpr) -> Result<Selection> {
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
    pub fn subsel_from_str(&self, sel_str: &str) -> Result<Selection> {
        let expr = SelectionExpr::try_from(sel_str)?;
        self.subsel_from_expr(&expr)
    }

    /// Subselection from the range of local selection indexes
    pub fn subsel_from_local_range(
        &self,
        range: std::ops::Range<usize>,
    ) -> Result<Selection> {
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
    ) -> Result<Selection> {
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

    pub fn query<'a>(&'a self) -> SelectionQueryGuard<'a> {
        SelectionQueryGuard {
            topology_ref: self.topology.read(),
            state_ref: self.state.read(),
            index: &self.index,
        }
    }

    pub fn modify<'a>(&'a self) -> SelectionModifyGuard<'a> {
        SelectionModifyGuard {
            topology_ref: self.topology.write(),
            state_ref: self.state.write(),
            index: &self.index,
        }
    }


}
//---------------,-------------------------------------

/// Scoped guard giving read-only access to selection
pub struct SelectionQueryGuard<'a> {
    topology_ref: Ref<'a,Topology>,
    state_ref: Ref<'a,State>,
    index: &'a Vec<usize>,
}

impl<'a> IoIndexProvider for SelectionQueryGuard<'a> {
    fn get_index(&self) -> impl IndexIterator {
        self.index.iter().cloned()
    }
}

impl<'a> IoTopologyProvider for SelectionQueryGuard<'a> {
    #[allow(refining_impl_trait)]
    fn get_topology(&self) -> &Topology {
        &self.topology_ref
    }
}

impl<'a> IoStateProvider for SelectionQueryGuard<'a> {
    #[allow(refining_impl_trait)]
    fn get_state(&self) -> &State {
        &self.state_ref
    }
}


/// Scoped guard giving read-write access to selection
pub struct SelectionModifyGuard<'a> {
    topology_ref: RefMut<'a,Topology>,
    state_ref: RefMut<'a,State>,
    index: &'a Vec<usize>,
}

//==================================================================
// Implement analysis traits

impl MeasureBox for SelectionQueryGuard<'_> {
    fn get_box(&self) -> Result<&PeriodicBox> {
        self.state_ref.get_box()
    }
}

impl MeasurePos for SelectionQueryGuard<'_> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.index.iter().map(|i| &self.state_ref.coords[*i])
    }
}

impl MeasureAtoms for SelectionQueryGuard<'_> {
    fn iter_atoms(&self) -> impl super::AtomIterator<'_> {
        self.index.iter().map(|i| &self.topology_ref.atoms[*i])
    }    
}

impl MeasureMasses for SelectionQueryGuard<'_> {
    fn iter_masses(&self) -> impl ExactSizeIterator<Item = f32> {
        self.index.iter().map(|i| self.topology_ref.atoms[*i].mass)
    }
}

impl MeasurePeriodic for SelectionQueryGuard<'_> {}

impl MeasureBox for SelectionModifyGuard<'_> {
    fn get_box(&self) -> Result<&PeriodicBox> {
        self.state_ref.get_box()
    }
}

//-------------------------------------------------------

impl MeasurePos for SelectionModifyGuard<'_> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.index.iter().map(|i| &self.state_ref.coords[*i])
    }
}

impl ModifyPos for SelectionModifyGuard<'_> {
    fn iter_pos_mut(&mut self) -> impl PosMutIterator<'_> {
        self.index.iter().map(|i|
            unsafe {
                &mut *self.state_ref.coords.as_mut_ptr().add(*i)
            }
        )
    }
}


/*
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
*/

impl ModifyPeriodic for SelectionModifyGuard<'_> {}

impl ModifyRandomAccess for SelectionModifyGuard<'_> {
    /*
    fn nth_particle_mut(&mut self, i: usize) -> ParticleMut {
        ParticleMut{
            id: i,
            pos: &mut self.state_ref.coords[self.index[i]],
            atom: &mut self.topology_ref.atoms[self.index[i]],
        }
    }
    */

    #[inline(always)]
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
        core::{selection::Select, Topology, MeasureMasses, MeasurePos, Vector3f, ModifyRandomAccess, fit_transform, PBC_FULL, ModifyPos},
        io::*,
    };
    use lazy_static::lazy_static;
    use nalgebra::Unit;

    use super::{Selection, SelectionAll};

    fn read_test_pdb() -> (Topology, State) {
        let mut h = FileHandler::open("tests/no_ATP.pdb").unwrap();
        //let top = h.read_topology().unwrap();
        //let state = h.read_state().unwrap().unwrap();
        //(top, state)
        h.read().unwrap()
    }

    // Read the test PDB file once and provide the content for tests
    lazy_static! {
        static ref SS: (Topology, State) = read_test_pdb();
    }

    fn make_sel() -> anyhow::Result<Selection> {
        let t = SS.0.clone().to_rc();
        let s = SS.1.clone().to_rc();
        let sel = SelectionAll{}.select(&t, &s)?;
        Ok(sel)
    }

    fn make_sel_prot() -> anyhow::Result<Selection> {
        let t = SS.0.clone().to_rc();
        let s = SS.1.clone().to_rc();
        let sel = "not resname TIP3 POT CLA".select(&t, &s)?;
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

        let mut h = FileHandler::create("f.pdb")?;
        h.write_topology(&q)?;
        h.write_state(&q)?;
        Ok(())
    }

    #[test]
    fn test_unwrap_connectivity_1() -> anyhow::Result<()> {
        let sel = make_sel_prot()?;
        sel.modify().unwrap_connectivity_dim(0.2, &PBC_FULL)?;
        
        let mut h = FileHandler::create("unwrapped.pdb")?;
        let q = sel.query();
        h.write_topology(&q)?;
        h.write_state(&q)?;
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
        
        let mut h = FileHandler::create("sel2.pdb")?;
        let q = sel2.query();
        h.write_topology(&q)?;
        h.write_state(&q)?;
        drop(q);

        let mut h = FileHandler::create("sel1_before.pdb")?;
        let q = sel1.query();
        h.write_topology(&q)?;
        h.write_state(&q)?;
        drop(q);


        let m = fit_transform(sel1.query(), sel2.query())?;
        println!("{m}");
        sel1.modify().apply_transform(&m);


        let mut h = FileHandler::create("sel1_after.pdb")?;
        let q = sel1.query();
        h.write_topology(&q)?;
        h.write_state(&q)?;

        Ok(())
    }
}
