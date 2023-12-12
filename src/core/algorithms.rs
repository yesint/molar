use anyhow::{bail, Result};
use num_traits::Zero;
use num_traits::Bounded;

use super::{
    ParticleIterator, ParticleMutIterator,
    Pos, Vector3f, PbcDims, PeriodicBox, IdPosIterator, PosIterator,
};

fn min_max<'a>(coords: impl PosIterator<'a>) -> (Pos,Pos) {
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

fn center_of_geometry<'a>(particles: impl PosIterator<'a>) -> Result<Pos> {
    let n = particles.len();
    if n > 0 {
        let c = particles.fold(Pos::new(0.0, 0.0, 0.0), |acc, el| acc + el.coords);
        Ok(c / n as f32)
    } else {
        bail!("Selection is empty")
    }
}

fn center_of_mass<'a>(particles: impl ParticleIterator<'a>) -> Result<Pos> {
    let n = particles.len();
    if n > 0 {
        let mut c = particles.fold(
            (Vector3f::zero(),0.0), 
            |acc, el| {
                (acc.0+el.pos.coords * el.atom.mass, acc.1+el.atom.mass)
            });
        c.0 /= c.1;
        Ok(Pos::from(c.0))
    } else {
        bail!("Selection is empty")
    }
}
