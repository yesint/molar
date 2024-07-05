use std::cell::{RefCell, RefMut, Ref};

use anyhow::{Result,bail};

use crate::prelude::*;

pub struct SelChecked {
    sel: RefCell<Sel<MutableSerial>>,
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

    pub fn get<'a>(&'a self) -> Result<Ref<'a, Sel<MutableSerial>>>{
        // Check validity of selection
        let ret = self.sel.borrow();
        Self::check_validity(&ret)?;
        Ok(ret)
    }

    pub fn get_mut<'a>(&'a self) -> Result<RefMut<'a, Sel<MutableSerial>>>{
        // Check validity of selection
        let ret = self.sel.borrow_mut();
        Self::check_validity(&ret)?;
        Ok(ret)
    }
}

struct Builder {
    src: RefCell<Source>,
}