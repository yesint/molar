use crate::prelude::*;

//--------------------------------------------------------------
// Immutable providers
//--------------------------------------------------------------
pub trait LenProvider {
    fn len(&self) -> usize;
}

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

pub trait ParticleProvider: IndexProvider {
    fn iter_particle(&self) -> impl ExactSizeIterator<Item = Particle<'_>>;
}

pub trait RandomPos: LenProvider {       
    unsafe fn nth_pos_unchecked(&self, i: usize) -> &Pos;
    
    fn nth_pos(&self, i: usize) -> Option<&Pos> {
        if i<self.len() {
            Some(unsafe{self.nth_pos_unchecked(i)})
        } else {
            None
        }
    }
}

pub trait RandomAtom: LenProvider {
    unsafe fn nth_atom_unchecked(&self, i: usize) -> &Atom;

    fn nth_atom(&self, i: usize) -> Option<&Atom> {
        if i<self.len() {
            Some(unsafe{self.nth_atom_unchecked(i)})
        } else {
            None
        }
    }
}

pub trait BondsProvider: LenProvider {

    fn num_bonds(&self) -> usize;

    
    fn iter_bonds(&self) -> impl Iterator<Item = &[usize; 2]>;

    unsafe fn nth_bond_unchecked(&self, i: usize) -> &[usize;2];
    
    fn nth_bond(&self, i: usize) -> Option<&[usize; 2]> {
        if i<self.num_bonds() {
            Some(unsafe{self.nth_bond_unchecked(i)})
        } else {
            None
        }
    }
}

pub trait MoleculesProvider: LenProvider {
    fn num_molecules(&self) -> usize;
    
    fn iter_molecules(&self) -> impl Iterator<Item = &[usize; 2]>;

    unsafe fn nth_molecule_unchecked(&self, i: usize) -> &[usize;2];
    
    fn nth_molecule(&self, i: usize) -> Option<&[usize; 2]> {
        if i<self.num_molecules() {
            Some(unsafe{self.nth_molecule_unchecked(i)})
        } else {
            None
        }
    }
}


//--------------------------------------------------------------
// Mutable providers
//--------------------------------------------------------------

pub trait PosMutProvider: PosProvider {
    fn iter_pos_mut(&self) -> impl PosMutIterator<'_>;
}

pub trait RandomPosMut: RandomPos + PosMutProvider {
    fn nth_pos_mut(&self, i: usize) -> Option<&mut Pos>;
    unsafe fn nth_pos_mut_unchecked(&self, i: usize) -> &mut Pos;
}

pub trait AtomsMutProvider: AtomsProvider {
    fn iter_atoms_mut(&self) -> impl AtomMutIterator<'_>;
}

pub trait ParticleMutProvider: IndexProvider {
    fn iter_particle_mut(&self) -> impl ExactSizeIterator<Item = ParticleMut<'_>>;
}

pub trait RandomAtomMut {
    fn nth_atom_mut(&self, i: usize) -> Option<&mut Atom>;
    unsafe fn nth_atom_mut_unchecked(&self, i: usize) -> &mut Atom;
}

pub trait BoxMutProvider {
    fn get_box_mut(&self) -> Option<&mut PeriodicBox>;
}

//----------------------------------------------------
impl PosProvider for Vec<Pos> {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.iter()
    }
}

impl PosProvider for Pos {
    fn iter_pos(&self) -> impl PosIterator<'_> {
        std::iter::once(self)
    }
}
