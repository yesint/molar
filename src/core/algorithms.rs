use anyhow::{bail, Result};
use num_traits::Zero;
use num_traits::Bounded;

use super::ParticleIteratorAdaptor;
use super::{
    ParticleIterator,
    Pos, Vector3f, PosIterator,
};


trait Measure {
    /* TODO: Waiting for Rust 1.74 to stabilize this feature
    fn iter(&self) -> impl ParticleIterator<'a>;
    
    fn iter_pos(&self) -> impl PosIterator<'_> {
        self.iter().map(|p| p.pos)
    }

    fn iter_atoms(&self) -> impl ExactSizeIterator<Item = &'_ Atom> {
        self.iter().map(|p| p.atom)
    }
    */
}

pub fn min_max<'a>(coords: impl PosIterator<'a>) -> (Pos,Pos) {
    let mut lower = Pos::max_value();
    let mut upper = Pos::min_value();
    for p in coords {
        for d in 0..3 {
            if p[d] < lower[d] { lower[d] = p[d] }
            if p[d] > upper[d] { upper[d] = p[d] }
        }
    }
    (lower,upper)
}

pub fn center_of_geometry<'a>(particles: impl PosIterator<'a>) -> Pos {
    let n = particles.len();
    let c = particles.fold(
        Pos::new(0.0, 0.0, 0.0),
        |acc, el| acc + el.coords
    );
    c / n as f32
}

pub fn center_of_mass<'a>(particles: impl ParticleIterator<'a>) -> Result<Pos> {
    let mut c = particles.fold(
        (Vector3f::zero(),0.0), 
        |acc, p| {
            (acc.0+p.pos.coords * p.atom.mass, acc.1+p.atom.mass)
        }
    );
    
    if c.1==0.0 {
        bail!("Zero mass in COM!")
    } else {
        c.0 /= c.1;
        Ok(Pos::from(c.0))
    }
}
