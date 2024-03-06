use super::{AtomIterator, AtomMutIterator, PeriodicBox, Pos, PosIterator, PosMutIterator};

//--------------------------------------------------------------
// Immutable providers
//--------------------------------------------------------------
pub trait PosProvider {
    fn iter_pos(&self) -> impl PosIterator<'_>;
}

pub trait MassesProvider {
    fn iter_masses(&self) -> impl ExactSizeIterator<Item = f32>;
}

pub trait AtomsProvider {
    fn iter_atoms(&self) -> impl AtomIterator<'_>;
}

pub trait BoxProvider {
    fn get_box(&self) -> Option<&PeriodicBox>;
}

//--------------------------------------------------------------
// Mutable providers
//--------------------------------------------------------------

pub trait PosMutProvider {
    fn iter_pos_mut(&mut self) -> impl PosMutIterator<'_>;
}

pub trait RandomPosMutProvider {
    unsafe fn nth_pos_unchecked_mut(&mut self, i: usize) -> &mut Pos;

    unsafe fn nth_pos_unchecked(&mut self, i: usize) -> &Pos {
        self.nth_pos_unchecked_mut(i)
    }
}

pub trait AtomsMutProvider {
    fn iter_atoms_mut(&mut self) -> impl AtomMutIterator<'_>;
}
