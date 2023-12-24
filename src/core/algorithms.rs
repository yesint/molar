use anyhow::{bail, Result};
use num_traits::Zero;
use num_traits::Bounded;

use super::Atom;
use super::AtomIterator;
use super::AtomMutIterator;
use super::ParticleIteratorAdaptor;
use super::ParticleMutIterator;
use super::PeriodicBox;
use super::PosMutIterator;
use super::{
    ParticleIterator,
    Pos, Vector3f, PosIterator,
};
 
// Trait that provides periodic box information
pub trait BoxProvider {
    fn get_box(&self) -> Result<&PeriodicBox>;
}

/// Trait for measuring various properties that requires only
/// the iterator of particles. User types should 
/// implement `iter`
pub trait Measure {
    fn iter(&self) -> impl ParticleIterator<'_>;
    
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.iter().map(|p| p.pos)
    }

    fn iter_atoms(&self) -> impl AtomIterator<'_> {
        self.iter().map(|p| p.atom)
    }

    fn min_max(&self) -> (Pos,Pos) {
        let mut lower = Pos::max_value();
        let mut upper = Pos::min_value();
        for p in self.iter_pos() {
            for d in 0..3 {
                if p[d] < lower[d] { lower[d] = p[d] }
                if p[d] > upper[d] { upper[d] = p[d] }
            }
        }
        (lower,upper)
    }

    fn center_of_geometry(&self) -> Pos {
        let iter = self.iter_pos();
        let n = iter.len();
        let c = iter.fold(
            Pos::new(0.0, 0.0, 0.0),
            |acc, el| acc + el.coords
        );
        c / n as f32
    }
    

}

/// The trait for measuring properties that requires
/// a periodic box information.
pub trait MeasurePeriodic: Measure + BoxProvider {
    fn center_of_mass(&self) -> Result<Pos> {
        let mut c = self.iter().fold(
            (Vector3f::zero(),0.0), 
            |acc, p| {
                (acc.0+p.pos.coords*p.atom.mass, acc.1+p.atom.mass)
            }
        );
        
        if c.1==0.0 {
            bail!("Zero mass in COM!")
        } else {
            Ok(Pos::from(c.0/c.1))
        }
    }
}

/// The trait for modifying the particles. User types should
/// implement `iter_mut`.
pub trait Modify {
    fn iter(&mut self) -> impl ParticleMutIterator<'_>;

    fn iter_pos(&mut self) -> impl PosMutIterator<'_> {
        self.iter().map(|p| p.pos)
    }

    fn iter_atoms(&mut self) -> impl AtomMutIterator<'_> {
        self.iter().map(|p| p.atom)
    }

    fn translate(&mut self, shift: Vector3f) {
        for el in self.iter_pos() {
            *el += shift;
        }
    }
}

/// The trait for modifying the particles that requires
/// the periodic box.
pub trait ModifyPeriodic: Modify + BoxProvider {
    
}
