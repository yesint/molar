use std::{cell::{Ref, RefCell}, ops::Deref, rc::Rc};

use anyhow::{anyhow, bail, Result};
use crate::prelude::*;

pub struct SelCheckedGuard<'a> {
    _src_ref: Ref<'a,Source>,
    sel_ref: &'a Sel<MutableSerial>,
}

impl Deref for SelCheckedGuard<'_> {
    type Target = Sel<MutableSerial>;
    fn deref(&self) -> &Self::Target {
        self.sel_ref
    }
}

pub struct SelChecked{
    sel: Sel<MutableSerial>,
    src: Rc<RefCell<Source>>,
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

    pub fn get(&self) -> Result<SelCheckedGuard<'_>>{
        // Check validity of selection        
        Self::check_validity(&self.sel)?;
        Ok(
            SelCheckedGuard {
                sel_ref: &self.sel,
                _src_ref: self.src.try_borrow()?,
            }
        )
    }
}

impl IndexProvider for SelChecked {
    fn iter_index(&self) -> impl Iterator<Item = usize> {
        self.sel.iter_index()
    }
}

impl AtomsProvider for SelChecked {
    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        self.sel.iter_atoms()
    }
}

impl PosProvider for SelChecked {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.sel.iter_pos()
    }
}

#[derive(Default)]
pub struct Builder(
    Rc<RefCell<Source>>,
);

impl Builder {
    pub fn new_from_source(src: Source) -> Self {
        Self( Rc::new(RefCell::new(src)) )
    }

    pub fn new(topology: TopologyUArc,state: StateUArc) -> Result<Self> {
        Ok(Self::new_from_source(Source::new(topology,state)?))
    }

    pub fn into_source(self) -> Result<Source> {
        Ok(Rc::into_inner(self.0).ok_or_else(|| anyhow!("There are still selections pointing to Builder!"))?.into_inner())
    }

    pub fn append(&mut self, sel: &(impl AtomsProvider+PosProvider)) -> Result<SelChecked> {
        let guard = self.0.try_borrow_mut()?;
        let top = guard.get_topology().get_mut();
        let st = guard.get_state().get_mut();

        let sz = top.atoms.len(); // Size before appending
        let added = sel.iter_atoms().len();
        top.add_atoms(sel.iter_atoms());
        st.add_coords(sel.iter_pos());
        // Make output selection
        Ok(SelChecked{
            sel: guard.select_range(&(sz..sz+added))?,
            src: Rc::clone(&self.0),
        })
    }

    pub fn remove(&mut self, sel: &impl IndexProvider) -> Result<()> {
        let guard = self.0.try_borrow_mut()?;
        let top = guard.get_topology().get_mut();
        let st = guard.get_state().get_mut();

        top.remove_atoms(sel.iter_index())?;
        st.remove_coords(sel.iter_index())?;
        Ok(())
    }

    pub fn len(&self) -> usize {
        self.0.borrow().get_topology().num_atoms()
    }
    //pub fn rearrange(&mut self, begin: &[&dyn IndexProvider], end: &[&dyn IndexIterator]) -> Result<()> {

    //}

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_builder() -> Result<()> {
        // Empty builder
        let mut b = Builder::default();
        // Read file and make source
        let (top,st) = FileHandler::open("tests/protein.pdb")?.read()?;
        let src = Source::new(top, st)?;
        // Create two overlapping selections
        let sel1 = src.select_from_iter(0..10)?;
        //let sel2 = src.select_from_iter(5..15)?;
        // Add to builder
        let added = b.append(&sel1)?;
        let _cm = added.get()?.center_of_mass();
        //b.append(&sel2)?;

        b.append(&added)?;
        assert_eq!(b.len(), 20);
        Ok(())
    }

    #[test]
    #[should_panic]
    fn test_builder_guard_fail() {
        // Read file and make source
        let (top,st) = FileHandler::open("tests/protein.pdb").unwrap().read().unwrap();
        let src = Source::new(top, st).unwrap();
        
        let mut b = Builder::default();

        // Create two overlapping selections
        let sel1 = src.select_from_iter(0..10).unwrap();
        
        // Add to builder
        let added = b.append(&sel1).unwrap();
        let _guard = added.get().unwrap();

        // This have to fail
        b.append(&sel1).unwrap();
    }
}