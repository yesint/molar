use std::{cell::{Ref, RefCell, RefMut}, rc::Rc};

use anyhow::Result;
use crate::prelude::*;

use super::utils::*;


//-----------------------------------------------------


pub struct BuilderSel {
    system: Rc<RefCell<System>>,
    index: Vec<usize>,
}

pub struct BuilderSelGuard<'a>( Ref<'a,System> );
pub struct BuilderSelGuardMut<'a>( RefMut<'a,System> );

impl BuilderSel {
    fn check_index(&self,n: usize) {
        if self.index[0] >= n || self.index[self.index.len()-1] >= n {
            panic!(
                "BuilderSel indexes [{}:{}] are out of allowed range [0:{}]",
                self.index[0],self.index[self.index.len()-1],n
            );
        }
    }

    pub fn guard(&self) -> BuilderSelGuard<'_> {
        let g = self.system.borrow();        
        let n = g.topology.num_atoms();
        self.check_index(n);
        BuilderSelGuard(g)
    }

    pub fn guard_mut(&self) -> BuilderSelGuardMut<'_> {
        let g = self.system.borrow_mut();
        let n = g.topology.num_atoms();
        self.check_index(n);
        BuilderSelGuardMut(g)
    }
}

//---------------------------------------------
impl PosProvider for BuilderSelGuard<'_> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.0.state.iter_pos()
    }
}

impl MeasurePos for BuilderSelGuard<'_> {}

impl AtomsProvider for BuilderSelGuard<'_> {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        self.0.topology.iter_atoms()
    }
}

impl TopologyProvider for BuilderSelGuard<'_> {
    fn num_atoms(&self) -> usize {
        self.0.topology.num_atoms()
    }
}

impl BoxProvider for BuilderSelGuard<'_> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.0.state.get_box()
    }
}

impl StateProvider for BuilderSelGuard<'_> {
    fn get_time(&self) -> f32 {
        self.0.state.get_time()
    }

    fn num_coords(&self) -> usize {
        self.0.state.num_coords()
    }
}
//---------------------------------------------
impl PosProvider for BuilderSelGuardMut<'_> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.0.state.iter_pos()
    }
}

impl MeasurePos for BuilderSelGuardMut<'_> {}

impl AtomsProvider for BuilderSelGuardMut<'_> {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        self.0.topology.iter_atoms()
    }
}

impl TopologyProvider for BuilderSelGuardMut<'_> {
    fn num_atoms(&self) -> usize {
        self.0.topology.num_atoms()
    }
}

impl BoxProvider for BuilderSelGuardMut<'_> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.0.state.get_box()
    }
}

impl StateProvider for BuilderSelGuardMut<'_> {
    fn get_time(&self) -> f32 {
        self.0.state.get_time()
    }

    fn num_coords(&self) -> usize {
        self.0.state.num_coords()
    }
}
//----------------------------------------------

pub struct Builder {
    system: Rc<RefCell<System>>,
}

impl Builder {
    pub fn new(topology: Topology,state: State) -> Result<Self> {
        Ok(Self{
            system: Rc::new(RefCell::new(System{topology,state}))
        })
    }

    pub fn select_all(&self) -> Result<BuilderSel> {
        let guard = self.system.borrow();
        let vec = index_from_all(guard.topology.num_atoms());
        Ok(
            BuilderSel {
                system: Rc::clone(&self.system),
                index: vec,
            }
        )
    }

    pub fn select_str(&self, selstr: &str) -> Result<BuilderSel> {
        let guard = self.system.borrow();
        let vec = index_from_str(selstr, &guard.topology, &guard.state)?;
        Ok(
            BuilderSel {
                system: Rc::clone(&self.system),
                index: vec,
            }
        )
    }

    pub fn select_expr(&self, expr: &SelectionExpr) -> Result<BuilderSel> {
        let guard = self.system.borrow();
        let vec = index_from_expr(expr, &guard.topology, &guard.state)?;
        Ok(
            BuilderSel {
                system: Rc::clone(&self.system),
                index: vec,
            }
        )
    }

    pub fn select_range(&self, range: &std::ops::Range<usize>) -> Result<BuilderSel> {
        let guard = self.system.borrow();
        let vec = index_from_range(range, guard.topology.num_atoms())?;
        Ok(
            BuilderSel {
                system: Rc::clone(&self.system),
                index: vec,
            }
        )
    }

    pub fn append(&mut self, data: &(impl TopologyProvider + StateProvider)) {
        let guard = self.system.borrow_mut();
        guard.topology.get_storage_mut().add_atoms(data.iter_atoms());
        guard.state.get_storage_mut().add_coords(data.iter_pos());
    }

    pub fn remove(&mut self, to_remove: &impl IndexProvider) -> Result<()> {
        let guard = self.system.borrow_mut();
        guard.topology.get_storage_mut().remove_atoms(to_remove.iter_index())?;
        guard.state.get_storage_mut().remove_coords(to_remove.iter_index())?;
        Ok(())
    }


}

//-----------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;

   
    #[test]
    fn test_new_sel() -> Result<()> {
        let (top,st) = FileHandler::open("tests/protein.pdb")?.read()?;
        
        let mut bu = Builder::new(top,st)?;
        let sel1 = bu.select_str("index 1:20")?;
        let cm = sel1.guard().center_of_geometry();
        //bu.append(&sel1);
        Ok(())
    }
}