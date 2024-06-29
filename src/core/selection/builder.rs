use std::cell::{RefCell, RefMut};

use crate::prelude::*;

struct SelChecked {
    sel: RefCell<Sel>,
}

impl SelChecked {
    fn get<'a>(&'a self) -> RefMut<'a, Sel>{
        // Check validity here
        self.sel.borrow_mut()
    }
}

struct Builder {
    src: RefCell<Source>,
}