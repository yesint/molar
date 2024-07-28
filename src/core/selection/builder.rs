use std::{cell::{Ref, RefCell, RefMut}, rc::Rc};

use anyhow::Result;
use crate::prelude::*;

use super::utils::*;


//-----------------------------------------------------


pub struct BuilderSel<'a> {
    system: &'a RefMut<'a,System>,
    index: Vec<usize>,
}

impl BuilderSel<'_> {
    fn check_index(&self,n: usize) {
        if self.index[0] >= n || self.index[self.index.len()-1] >= n {
            panic!(
                "BuilderSel indexes [{}:{}] are out of allowed range [0:{}]",
                self.index[0],self.index[self.index.len()-1],n
            );
        }
    }
}

//---------------------------------------------
impl PosProvider for BuilderSel<'_> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        unsafe { self.index.iter().map(|i| self.system.state.nth_pos_unchecked(*i)) }
    }
}

impl MeasurePos for BuilderSel<'_> {}

impl AtomsProvider for BuilderSel<'_> {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        unsafe {
            self.index
                .iter()
                .map(|i| self.system.topology.nth_atom_unchecked(*i))
        }
    }
}

impl TopologyProvider for BuilderSel<'_> {
    fn num_atoms(&self) -> usize {
        self.system.topology.num_atoms()
    }
}

impl BoxProvider for BuilderSel<'_> {
    fn get_box(&self) -> Option<&PeriodicBox> {
        self.system.state.get_box()
    }
}

impl StateProvider for BuilderSel<'_> {
    fn get_time(&self) -> f32 {
        self.system.state.get_time()
    }

    fn num_coords(&self) -> usize {
        self.system.state.num_coords()
    }
}
//----------------------------------------------

pub struct Builder {
    system: Rc<RefCell<System>>,
}
pub struct BuilderGuard<'a>{
    sys_ref: Ref<'a,System>,
    sys_rc: &'a Rc<RefCell<System>>,
}
pub struct BuilderGuardMut<'a>{
    sys_ref: RefMut<'a,System>,
    sys_rc: &'a Rc<RefCell<System>>,
}


impl Builder {
    pub fn new(topology: Topology,state: State) -> Result<Self> {
        Ok(Self{
            system: Rc::new(RefCell::new(System{topology,state}))
        })
    }

    pub fn guard(&self) -> BuilderGuard<'_> {
        let g = self.system.borrow();
        BuilderGuard{sys_ref: g, sys_rc: &self.system}
    }

    pub fn guard_mut(&self) -> BuilderGuardMut<'_> {
        let g = self.system.borrow_mut();
        BuilderGuardMut{sys_ref: g, sys_rc: &self.system}
    }
}

impl BuilderGuardMut<'_> {
    pub fn select_all(&self) -> Result<BuilderSel> {
        let vec = index_from_all(self.sys_ref.topology.num_atoms());
        Ok(
            BuilderSel {
                system: &self.sys_ref,
                index: vec,
            }
        )
    }

    pub fn select_str(&self, selstr: &str) -> Result<BuilderSel> {
        let vec = index_from_str(selstr, &self.sys_ref.topology, &self.sys_ref.state)?;
        Ok(
            BuilderSel {
                system: &self.sys_ref,
                index: vec,
            }
        )
    }

    pub fn select_expr(&self, expr: &SelectionExpr) -> Result<BuilderSel> {
        let vec = index_from_expr(expr, &self.sys_ref.topology, &self.sys_ref.state)?;
        Ok(
            BuilderSel {
                system: &self.sys_ref,
                index: vec,
            }
        )
    }

    pub fn select_range(&self, range: &std::ops::Range<usize>) -> Result<BuilderSel> {
        let vec = index_from_range(range, self.sys_ref.topology.num_atoms())?;
        Ok(
            BuilderSel {
                system: &self.sys_ref,
                index: vec,
            }
        )
    }

    pub fn append(&self, data: &(impl TopologyProvider + StateProvider)) {
        self.sys_ref.topology.get_storage_mut().add_atoms(data.iter_atoms());
        self.sys_ref.state.get_storage_mut().add_coords(data.iter_pos());
    }

    pub fn remove(&self, to_remove: &impl IndexProvider) -> Result<()> {
        self.sys_ref.topology.get_storage_mut().remove_atoms(to_remove.iter_index())?;
        self.sys_ref.state.get_storage_mut().remove_coords(to_remove.iter_index())?;
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
        let g = bu.guard_mut();
        let sel1 = g.select_str("index 1:20")?;
        let cm = sel1.center_of_geometry();
        g.append(&sel1);
        Ok(())
    }
}