use std::iter::zip;

use super::{Matrix3f, PBC_FULL};
use super::{
    AtomIterator, PbcDims,
    PeriodicBox, Pos, PosIterator, PosMutIterator, Vector3f,
};
use crate::distance_search::search::{DistanceSearcherSingle, SearchConnectivity};
use anyhow::{bail, Result};
use itertools::izip;
use nalgebra::Rotation3;
use nalgebra::Unit;
use nalgebra::SVD;
use num_traits::Bounded;
use num_traits::Zero;

// Trait that provides periodic box information
pub trait MeasureBox {
    fn get_box(&self) -> Result<&PeriodicBox>;
}

/// Trait for measuring various properties that requires only
/// the iterator of positions.
pub trait MeasurePos {
    fn iter_pos(&self) -> impl PosIterator<'_>;

    fn min_max(&self) -> (Pos, Pos) {
        let mut lower = Pos::max_value();
        let mut upper = Pos::min_value();
        for p in self.iter_pos() {
            for d in 0..3 {
                if p[d] < lower[d] {
                    lower[d] = p[d]
                }
                if p[d] > upper[d] {
                    upper[d] = p[d]
                }
            }
        }
        (lower, upper)
    }

    fn center_of_geometry(&self) -> Pos {
        let iter = self.iter_pos();
        let n = iter.len();
        let mut cog = Vector3f::zero();
        for c in iter {
            cog += c.coords;
        }
        Pos::from(cog / n as f32)
    }
}

pub trait MeasureAtoms {
    fn iter_atoms(&self) -> impl AtomIterator<'_>;
}

/// Trait for measuring various properties that requires only
/// the iterator of particles. User types should
/// implement `iter`
pub trait MeasureMasses: MeasurePos {
    fn iter_masses(&self) -> impl ExactSizeIterator<Item = f32>;

    fn center_of_mass(&self) -> Result<Pos> {
        let mut cm = Vector3f::zero();
        let mut mass = 0.0;
        for (c, m) in zip(self.iter_pos(), self.iter_masses()) {
            cm += c.coords * m;
            mass += m;
        }

        if mass == 0.0 {
            bail!("Zero mass in COM!")
        } else {
            Ok(Pos::from(cm / mass))
        }
    }
}

/// The trait for measuring properties that requires
/// a periodic box information.
pub trait MeasurePeriodic: MeasureMasses + MeasureBox {
    fn center_of_mass_pbc(&self) -> Result<Pos> {
        let b = self.get_box()?;
        let mut pos_iter = self.iter_pos();
        let mut mass_iter = self.iter_masses();
        
        let mut mass = mass_iter.next().unwrap();
        let p0 = pos_iter.next().unwrap();
        let mut cm = p0.coords;

        for (c, m) in zip(pos_iter, mass_iter) {
            let im = b.closest_image(c, p0).coords;
            cm += im * m;
            mass += m;
        }

        if mass == 0.0 {
            bail!("Zero mass in COM!")
        } else {
            Ok(Pos::from(cm / mass))
        }
    }
}

pub trait ModifyPos: MeasurePos {
    fn iter_pos_mut(&mut self) -> impl PosMutIterator<'_>;

    fn translate(&mut self, shift: Vector3f) {
        for el in self.iter_pos_mut() {
            *el += shift;
        }
    }

    fn rotate(&mut self, ax: &Unit<Vector3f>, ang: f32) {
        let tr = Rotation3::<f32>::from_axis_angle(ax, ang);
        for p in self.iter_pos_mut() {
            p.coords = tr * p.coords;
        }
    }

    fn apply_transform(&mut self, tr: &nalgebra::IsometryMatrix3<f32>) {
        for p in self.iter_pos_mut() {
            *p = tr * (*p);
        }
    }
}

/// The trait for modifying the particles. User types should
/// implement `iter_mut`.
/*
pub trait ModifyParticles: ModifyPos {
    fn iter_particles_mut(&mut self) -> impl ParticleMutIterator<'_>;

    fn iter_particles(&mut self) -> impl ParticleIterator<'_> {
        self.iter_particles_mut().map(|p| p.into())
    }

    fn iter_atoms_mut(&mut self) -> impl AtomMutIterator<'_> {
        self.iter_particles_mut().map(|p| p.atom)
    }
}
*/

/// The trait for modifying the particles that requires
/// the periodic box.
pub trait ModifyPeriodic: ModifyPos + MeasureBox {
    fn unwrap_simple_dim(&mut self, dims: PbcDims) -> Result<()> {
        let b = self.get_box()?.clone();
        let mut iter = self.iter_pos_mut();
        if iter.len() > 0 {
            let p0 = iter.next().unwrap();
            for p in iter {
                *p = b.closest_image_dims(p, p0, &dims);
            }
        }
        Ok(())
    }

    fn unwrap_simple(&mut self) -> Result<()> {
        self.unwrap_simple_dim(PBC_FULL)
    }
}

pub trait ModifyRandomAccess: ModifyPeriodic {
    //fn nth_particle_mut(&mut self, i: usize) -> ParticleMut;
    fn nth_pos_mut(&mut self, i: usize) -> &mut Pos;

    fn nth_pos(&mut self, i: usize) -> &Pos {
        self.nth_pos_mut(i)
    }

    fn unwrap_connectivity(&mut self, cutoff: f32) -> Result<()> {
        self.unwrap_connectivity_dim(cutoff, &PBC_FULL)
    }

    fn unwrap_connectivity_dim(&mut self, cutoff: f32, dims: &PbcDims) -> Result<()> {
        let b = self.get_box()?.to_owned();
        let conn: SearchConnectivity = DistanceSearcherSingle::new_periodic(
            cutoff,
            self.iter_pos().enumerate(),
            &b,
            &dims,
        )
        .search();

        // used atoms
        let mut used = vec![false; conn.len()];
        // Centers to unwrap
        let mut todo = Vec::<usize>::with_capacity(conn.len() / 2);
        // Place first center to the stack
        todo.push(0);
        used[0] = true;

        // Loop while stack is not empty
        while let Some(c) = todo.pop() {
            // Central point
            let p0 = self.nth_pos(c).to_owned();
            // Iterate over connected points
            for ind in &conn[c] {
                // Unwrap this point if it is not used yet
                if !used[*ind] {
                    let p = self.nth_pos_mut(*ind);
                    *p = b.closest_image_dims(p, &p0, &dims);
                    // Add it to the stack
                    todo.push(*ind);
                    used[*ind] = true;
                }
            }
            //println!(">> {:?}",todo);
        }

        if used.len() != conn.len() {
            bail!("Selection is not compact for cutoff={}", cutoff)
        }

        Ok(())
    }
}

// Straighforward implementation of Kabsch algorithm
pub fn rot_transform<'a>(
    pos1: impl Iterator<Item = Vector3f>,
    pos2: impl Iterator<Item = Vector3f>,
    masses: impl Iterator<Item = f32>,
) -> Rotation3<f32> {
    //Calculate the covariance matrix
    let mut cov = Matrix3f::zeros();

    for (p1, p2, m) in izip!(pos1, pos2, masses) {
        cov += p2 * p1.transpose() * m;
    }

    // Perform Singular Value Decomposition (SVD) on the covariance matrix
    let svd = SVD::new(cov, true, true);
    let u = svd.u.unwrap();
    let v_t = svd.v_t.unwrap();

    // Determine if a reflection is necessary
    let d = if (u * v_t).determinant() < 0.0 {
        -1.0
    } else {
        1.0
    };

    // Create a diagonal matrix for correcting the reflection
    let mut d_matrix = Matrix3f::identity();
    d_matrix[(2, 2)] = d;

    // Compute the optimal rotation matrix
    Rotation3::from_matrix_unchecked(u * d_matrix * v_t)
}

pub fn fit_transform(
    sel1: impl MeasureMasses,
    sel2: impl MeasureMasses,
) -> Result<nalgebra::IsometryMatrix3<f32>> {
    let cm1 = sel1.center_of_mass()?;
    let cm2 = sel2.center_of_mass()?;

    //let rot = rot_transform_matrix(coords1.iter(), coords2.iter(), masses.iter());
    let rot = rot_transform(
        sel1.iter_pos().map(|p| *p-cm1),
        sel2.iter_pos().map(|p| *p-cm2),
        sel1.iter_masses()
    );

    Ok(nalgebra::Translation3::from(cm2) * rot * nalgebra::Translation3::from(-cm1))
}

// Version for selection with CM alredy at zero
pub fn fit_transform_at_origin(
    sel1: impl MeasureMasses,
    sel2: impl MeasureMasses,
) -> Result<nalgebra::IsometryMatrix3<f32>> {
    
    let rot = rot_transform(
        sel1.iter_pos().map(|p| p.coords),
        sel2.iter_pos().map(|p| p.coords),
        sel1.iter_masses()
    );

    Ok(nalgebra::convert(rot))
}

/// Mass-weighted RMSD
pub fn rmsd_mw(sel1: impl MeasureMasses,sel2: impl MeasureMasses) -> Result<f32> {
    let mut res = 0.0;
    let mut m_tot = 0.0;
    let iter1 = sel1.iter_pos();
    let iter2 = sel2.iter_pos();

    if iter1.len() != iter2.len() {
        bail!("Different sizes in rmsd_mw: {} and {}",iter1.len(),iter2.len());
    }

    for (p1,p2,m) in izip!(iter1,iter2,sel1.iter_masses()){
        res += (p2-p1).norm_squared()*m;
        m_tot += m;
    }

    if m_tot==0.0 {
        bail!("Zero mass in rmsd_mw")
    } else {
        Ok((res/m_tot).sqrt())
    }
}

/// RMSD
pub fn rmsd(sel1: impl MeasurePos,sel2: impl MeasurePos) -> Result<f32> {
    let mut res = 0.0;
    let iter1 = sel1.iter_pos();
    let iter2 = sel2.iter_pos();

    if iter1.len() != iter2.len() {
        bail!("Different sizes in rmsd: {} and {}",iter1.len(),iter2.len());
    }

    let n = iter1.len();
    if n==0 {
        bail!("No atoms in rmsd")
    }

    for (p1,p2) in zip(iter1,iter2){
        res += (p2-p1).norm_squared();
    }

    Ok((res/n as f32).sqrt())
}