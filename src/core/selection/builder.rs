use anyhow::{Result,bail};
use crate::prelude::*;

pub struct SelChecked(Sel<MutableSerial>);

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
        Self::check_validity(&self.0)?;
        Ok(&self.0)
    }
}

pub struct Builder(Source);

impl Builder {
    pub fn new_from_source(src: Source) -> Self {
        Self(src)
    }

    pub fn new(topology: TopologyUArc,state: StateUArc) -> Result<Self> {
        Ok(Self::new_from_source(Source::new(topology,state)?))
    }

    #[inline(always)]
    fn get_storage_mut(&self) -> (&mut TopologyStorage, &mut StateStorage) {
        (self.0.get_topology().get_mut(), self.0.get_state().get_mut())
    }

    pub fn into_source(self) -> Source {
        self.0
    }

    pub fn append(&mut self, sel: &(impl AtomsProvider+PosProvider)) -> Result<SelChecked> {
        let (top,st) = self.get_storage_mut();
        let sz = top.atoms.len(); // Size before appending
        let added = sel.iter_atoms().len();
        top.add_atoms(sel.iter_atoms());
        st.add_coords(sel.iter_pos());
        // Make output selection
        Ok(SelChecked(self.0.select_range(&(sz+1..sz+added))?))
    }

    pub fn remove(&mut self, sel: &impl IndexProvider) -> Result<()> {
        let (top,st) = self.get_storage_mut();
        top.remove_atoms(sel.iter_index())?;
        st.remove_coords(sel.iter_index())?;
        Ok(())
    }

    //pub fn rearrange(&mut self, begin: &[&dyn IndexProvider], end: &[&dyn IndexIterator]) -> Result<()> {

    //}

}