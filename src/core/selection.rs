use std::cell::{Ref, RefMut};
use crate::io::{IoIndexProvider, IoTopologyProvider, IoStateProvider};
use super::{AtomIterator, AtomsProvider, BoxProvider, IndexIterator, MassesProvider, GuardedQuery, MeasureMasses, MeasurePeriodic, MeasurePos, GuardedModify, ModifyPos, ModifyRandomAccess, PeriodicBox, Pos, PosIterator, PosMutIterator, PosMutProvider, PosProvider, RandomPosMutProvider, State, StateRc, Topology, TopologyRc};
use anyhow::{bail, Result};
use itertools::Itertools;

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
        let st = state.borrow();
        check_sizes(&topology.borrow(), &st)?;
        let index: Vec<usize> = (0..st.coords.len()).collect();
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
        check_sizes(&topology.borrow(), &state.borrow())?;
        let index = self.apply_whole(&topology.borrow(), &state.borrow()).unwrap();
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
        let top = topology.borrow();
        let st = state.borrow();
        check_sizes(&top, &st)?;
        let index = SelectionExpr::try_from(self)
            .unwrap()
            .apply_whole(&top, &st)
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
        let top = topology.borrow();
        let st = state.borrow();
        check_sizes(&top, &st)?;
        let n = top.atoms.len();
        if self.start > n - 1 || self.end > n - 1 {
            bail!(
                "Index range {}:{} is invalid, 0:{} is allowed.",
                self.start,
                self.end,
                top.atoms.len()
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
pub struct Selection 
{
    topology: TopologyRc,
    state: StateRc,
    index: Vec<usize>,
}

impl Selection {
    /// Subselection from expression
    pub fn subsel_from_expr(&self, expr: &SelectionExpr) -> Result<Selection> {
        let index = expr.apply_subset(
            &self.topology.borrow(),
            &self.state.borrow(),
            self.index.iter().cloned()
        )?;
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
}
//---------------,-------------------------------------

/// Scoped guard giving read-only access to selection
pub struct SelectionQueryGuard<'a> {
    topology_ref: Ref<'a,Topology>,
    state_ref: Ref<'a,State>,
    index: &'a Vec<usize>,
}

/// Scoped guard giving read-write access to selection
#[allow(dead_code)]
pub struct SelectionModifyGuard<'a> {
    topology_ref: RefMut<'a,Topology>,
    state_ref: RefMut<'a,State>,
    index: &'a Vec<usize>,
}

//---------------------------------------------
// Implement traits for IO

impl<'a> IoIndexProvider for SelectionQueryGuard<'a> {
    fn get_index(&self) -> impl IndexIterator {
        self.index.iter().cloned()
    }
}

impl<'a> IoTopologyProvider for SelectionQueryGuard<'a> {
    fn get_topology(&self) -> &Topology {
        &self.topology_ref
    }
}

impl<'a> IoStateProvider for SelectionQueryGuard<'a> {
    fn get_state(&self) -> &State {
        &self.state_ref
    }
}

//==================================================================
// Implement analysis traits

impl GuardedQuery for Selection {
    type Guard<'a> = SelectionQueryGuard<'a>;
    fn guard<'a>(&'a self) -> Self::Guard<'a> {
        SelectionQueryGuard {
            topology_ref: self.topology.borrow(),
            state_ref: self.state.borrow(),
            index: &self.index,
        }
    }
}

impl BoxProvider for SelectionQueryGuard<'_> {
    fn get_box(&self) -> Result<&PeriodicBox> {
        self.state_ref.get_box()
    }
}

impl MeasurePeriodic for Selection {}

impl PosProvider for SelectionQueryGuard<'_> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.index.iter().map(|i| &self.state_ref.coords[*i])
    }
}

impl MeasurePos for Selection {}

impl AtomsProvider for SelectionQueryGuard<'_> {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        self.index.iter().map(|i| &self.topology_ref.atoms[*i])
    }    
}

impl MassesProvider for SelectionQueryGuard<'_> {
    fn iter_masses(&self) -> impl ExactSizeIterator<Item = f32> {
        self.index.iter().map(|i| self.topology_ref.atoms[*i].mass)
    }
}

impl MeasureMasses for Selection {}

//-------------------------------------------------------
impl GuardedModify for Selection {
    type GuardMut<'a> = SelectionModifyGuard<'a>;
    fn guard_mut<'a>(&'a self) -> Self::GuardMut<'a> {
        SelectionModifyGuard {
            topology_ref: self.topology.borrow_mut(),
            state_ref: self.state.borrow_mut(),
            index: &self.index,
        }
    }
}

impl BoxProvider for SelectionModifyGuard<'_> {
    fn get_box(&self) -> Result<&PeriodicBox> {
        self.state_ref.get_box()
    }
}

impl PosMutProvider for SelectionModifyGuard<'_> {
    fn iter_pos_mut(&mut self) -> impl PosMutIterator<'_> {
        self.index.iter().map(|i|
            unsafe {
                &mut *self.state_ref.coords.as_mut_ptr().add(*i)
            }
        )
    }
}

impl PosProvider for SelectionModifyGuard<'_> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.index.iter().map(|i| &self.state_ref.coords[*i])
    }
}

impl RandomPosMutProvider for SelectionModifyGuard<'_> {
    #[inline(always)]
    fn nth_pos_mut(&mut self, i: usize) -> &mut Pos {
        &mut self.state_ref.coords[self.index[i]]
    }    
}

impl ModifyPos for Selection {}
impl ModifyRandomAccess for Selection {}

//############################################################
//#  Tests
//############################################################

#[cfg(test)]
mod tests {
    use crate::{
        core::{algorithms::*, selection::Select, State, Topology, Vector3f, PBC_FULL},
        io::*,
    };
    use lazy_static::lazy_static;

    use super::{Selection, SelectionAll};

    fn read_test_pdb() -> (Topology, State) {
        let mut h = FileHandler::open("tests/no_ATP.pdb").unwrap();
        //let top = h.read_topology().unwrap();
        //let state = h.read_state().unwrap().unwrap();
        //(top, state)
        let (top,st) = h.read_raw().unwrap();
        (top,st)
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
        println!("before {}",sel.guard().iter_pos().next().unwrap());

        let (minv,maxv) = sel.min_max();
        println!("{minv}:{maxv}");

        //sel.translate(&Vector3f::new(10.0,10.0,10.0));
        println!("after {}",sel.guard().iter_pos().next().unwrap());
        println!("{:?}",sel.min_max());
        Ok(())
    }
    
    #[test]
    fn test_measure_pbc() -> anyhow::Result<()>{
        let sel = make_sel()?;

        let cm = sel.center_of_mass()?;
        println!("{cm}");
        Ok(())
    }

    #[test]
    fn test_translate() -> anyhow::Result<()> {
        let sel = make_sel()?;

        println!("before {}",sel.guard_mut().iter_pos_mut().next().unwrap());
        sel.translate(Vector3f::new(10.0,10.0,10.0));
        println!("after {}",sel.guard_mut().iter_pos_mut().next().unwrap());
        Ok(())
    }

    #[test]
    fn test_write_to_file() -> anyhow::Result<()> {
        let sel = make_sel()?;
        
        let mut h = FileHandler::create("f.pdb")?;
        h.write_topology(&sel)?;
        h.write_state(&sel)?;
        Ok(())
    }

    #[test]
    fn test_unwrap_connectivity_1() -> anyhow::Result<()> {
        let sel = make_sel_prot()?;
        sel.unwrap_connectivity_dim(0.2, &PBC_FULL)?;
        
        let mut h = FileHandler::create("unwrapped.pdb")?;
        h.write_topology(&sel)?;
        h.write_state(&sel)?;
        Ok(())
    }

    #[test]
    fn eigen_test() -> anyhow::Result<()> {
        let sel1 = make_sel_prot()?;
        let sel2 = make_sel_prot()?;
        
        sel2.rotate(&Vector3f::x_axis(), 80.0_f32.to_radians());        

        let mut h = FileHandler::create("sel2.pdb")?;
        h.write_topology(&sel2)?;
        h.write_state(&sel2)?;

        let mut h = FileHandler::create("sel1_before.pdb")?;
        h.write_topology(&sel1)?;
        h.write_state(&sel1)?;
        
        let m = Selection::fit_transform(&sel1,&sel2)?;
        println!("{m}");

        sel1.apply_transform(&m);

        let mut h = FileHandler::create("sel1_after.pdb")?;
        h.write_topology(&sel1)?;
        h.write_state(&sel1)?;

        Ok(())
    }
}
