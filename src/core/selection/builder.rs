use std::{cell::{Ref, RefCell, RefMut}, rc::Rc};

use anyhow::Result;
use crate::prelude::*;

use super::utils::index_from_str;


//-----------------------------------------------------


struct SelSerial {
    system: Rc<RefCell<System>>,
    index: Vec<usize>,
}

struct SelGuard<'a>( Ref<'a,System> );
struct SelGuardMut<'a>( RefMut<'a,System> );

impl SelSerial {
    fn guard(&self) -> SelGuard<'_> {
        SelGuard(self.system.borrow())
    }

    fn guard_mut(&self) -> SelGuardMut<'_> {
        SelGuardMut(self.system.borrow_mut())
    }
}

impl PosProvider for SelGuard<'_> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.0.state.iter_pos()
    }
}

impl MeasurePos for SelGuard<'_> {}


struct SourceSerial {
    system: Rc<RefCell<System>>,
}

impl SourceSerial {
    fn new(topology: Topology,state: State) -> Result<Self> {
        Ok(Self{
            system: Rc::new(RefCell::new(System{topology,state}))
        })
    }

    pub fn select_str(&self, selstr: &str) -> Result<SelSerial> {
        let guard = self.system.borrow();
        let vec = index_from_str(selstr, &guard.topology, &guard.state)?;
        Ok(
            SelSerial {
                system: Rc::clone(&self.system),
                index: vec,
            }
        )
    }
}



//-----------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;

   
    #[test]
    fn test_new_sel() -> Result<()> {
        let (top,st) = FileHandler::open("tests/protein.pdb")?.read()?;
        
        let src = SourceSerial::new(top,st)?;
        let sel1 = src.select_str("index 1:20")?;
        let sel2 = src.select_str("index 30:120")?;
        let cm = sel1.guard().center_of_geometry();
        Ok(())
    }
}