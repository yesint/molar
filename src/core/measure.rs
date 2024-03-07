use std::iter::zip;

use super::algorithms::*;
use super::providers::*;
use super::Matrix3f;
use super::Pos;
use super::Vector3f;
use anyhow::{bail,Result,anyhow};
use itertools::izip;
use nalgebra::Rotation3;
use num_traits::Bounded;
//==============================================================
// Traits for measuring (immutable access)
//==============================================================

/// Trait for analysis requiring only positions
pub trait MeasurePos: PosProvider {
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
        let mut cog = Vector3f::zeros();
        for c in iter {
            cog += c.coords;
        }
        Pos::from(cog / n as f32)
    }

    fn rmsd(sel1: &Self, sel2: &Self) -> Result<f32> {
        let mut res = 0.0;
        let iter1 = sel1.iter_pos();
        let iter2 = sel2.iter_pos();

        if iter1.len() != iter2.len() {
            bail!(
                "Different sizes in rmsd: {} and {}",
                iter1.len(),
                iter2.len()
            );
        }

        let n = iter1.len();
        if n == 0 {
            bail!("No atoms in rmsd")
        }

        for (p1, p2) in std::iter::zip(iter1, iter2) {
            res += (p2 - p1).norm_squared();
        }

        Ok((res / n as f32).sqrt())
    }
}

/// Trait for analysis requiring positions and masses
pub trait MeasureMasses: PosProvider + MassesProvider {
    fn center_of_mass(&self) -> Result<Pos> {
        let mut cm = Vector3f::zeros();
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

    fn fit_transform(sel1: &Self, sel2: &Self) -> Result<nalgebra::IsometryMatrix3<f32>> {
        let cm1 = sel1.center_of_mass()?;
        let cm2 = sel2.center_of_mass()?;

        //let rot = rot_transform_matrix(coords1.iter(), coords2.iter(), masses.iter());
        let rot = rot_transform(
            sel1.iter_pos().map(|p| *p - cm1),
            sel2.iter_pos().map(|p| *p - cm2),
            sel1.iter_masses(),
        );

        Ok(nalgebra::Translation3::from(cm2) * rot * nalgebra::Translation3::from(-cm1))
    }

    fn fit_transform_at_origin(sel1: &Self, sel2: &Self) -> Result<nalgebra::IsometryMatrix3<f32>> {
        let rot = rot_transform(
            sel1.iter_pos().map(|p| p.coords),
            sel2.iter_pos().map(|p| p.coords),
            sel1.iter_masses(),
        );
        Ok(nalgebra::convert(rot))
    }

    fn rmsd_mw(sel1: &Self, sel2: &Self) -> Result<f32> {
        let mut res = 0.0;
        let mut m_tot = 0.0;
        let iter1 = sel1.iter_pos();
        let iter2 = sel2.iter_pos();

        if iter1.len() != iter2.len() {
            bail!(
                "Different sizes in rmsd_mw: {} and {}",
                iter1.len(),
                iter2.len()
            );
        }

        for (p1, p2, m) in izip!(iter1, iter2, sel1.iter_masses()) {
            res += (p2 - p1).norm_squared() * m;
            m_tot += m;
        }

        if m_tot == 0.0 {
            bail!("Zero mass in rmsd_mw")
        } else {
            Ok((res / m_tot).sqrt())
        }
    }
}

/// Trait for analysis requiring positions, masses and pbc
pub trait MeasurePeriodic: PosProvider + MassesProvider + BoxProvider {
    fn center_of_mass_pbc(&self) -> Result<Pos> {
        let b = self.get_box().ok_or_else(|| anyhow!("No periodicity!"))?;
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

// Straighforward implementation of Kabsch algorithm
pub fn rot_transform(
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
    let svd = nalgebra::SVD::new(cov, true, true);
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
