use std::cell::{RefCell, RefMut, Ref};

use anyhow::{Result,bail};

use crate::prelude::*;

pub struct SelChecked {
    sel: Sel<MutableSerial>,
}

impl SelChecked {
    fn check_validity(sel: &Sel) -> Result<()> {        
        let b = sel.first_index();
        let e = sel.last_index();
        if sel.num_atoms() != sel.num_coords() {
            bail!("Topology has {} atoms, while state has {} coordinates!",sel.num_atoms(),sel.num_coords());
        }        
        if b > sel.num_atoms() || e > sel.num_atoms() {
            bail!("Indexes [{}:{}] are out of allowed range [0:{}]",b,e,sel.len());
        }
        Ok(())
    }

    pub fn get(&self) -> Result<&Sel<MutableSerial>>{
        // Check validity of selection        
        Self::check_validity(&self.sel)?;
        Ok(&self.sel)
    }
}

pub struct Builder {
    topology: TopologyStorage,
    state: StateStorage,
}

impl Builder {
    pub fn clone_as_source(&self) -> Source {
        Source::new(
            self.topology.clone().into(),
            self.state.clone().into()
        ).unwrap()
    }

    
}