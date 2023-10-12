use anyhow::{bail, Result};
use num_traits::Zero;

use super::{
    selection::{ParticleIterator, ParticleMutIterator},
    Pos, Vector3f, PbcDims, PeriodicBox,
};

fn center_of_geometry<'a>(particles: impl ParticleIterator<'a>) -> Result<Pos> {
    let n = particles.len();
    if n > 0 {
        let c = particles.fold(Pos::new(0.0, 0.0, 0.0), |acc, el| acc + el.pos.coords);
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

