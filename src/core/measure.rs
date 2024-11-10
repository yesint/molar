use super::providers::*;
use super::Matrix3f;
use super::PbcDims;
use super::Pos;
use super::Vector3f;
use itertools::izip;
use nalgebra::SymmetricEigen;
use nalgebra::{IsometryMatrix3, Rotation3, Translation3};
use num_traits::Bounded;
use std::iter::zip;
use thiserror::Error;
//==============================================================
// Traits for measuring (immutable access)
//==============================================================

#[derive(Error, Debug)]
pub enum MeasureError {
    #[error("invalid data sizes: {0} and {1}")]
    Sizes(usize, usize),
    #[error("zero mass")]
    ZeroMass,
    #[error("SVD failed")]
    Svd,
    #[error("no periodic box")]
    NoPbc,
    #[error("can't unwrap disjoint pieces")]
    Disjoint,
}

/// Trait for analysis requiring only positions
pub trait MeasurePos: PosProvider + LenProvider {
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
        let n = self.len();
        let mut cog = Vector3f::zeros();
        for c in iter {
            cog += c.coords;
        }
        Pos::from(cog / n as f32)
    }

    fn rmsd(sel1: &Self, sel2: &Self) -> Result<f32, MeasureError> {
        let mut res = 0.0;
        let iter1 = sel1.iter_pos();
        let iter2 = sel2.iter_pos();

        if sel1.len() != sel2.len() {
            return Err(MeasureError::Sizes(sel1.len(), sel2.len()));
        }

        let n = sel1.len();
        if n == 0 {
            return Err(MeasureError::Sizes(sel1.len(), sel2.len()));
        }

        for (p1, p2) in std::iter::zip(iter1, iter2) {
            res += (p2 - p1).norm_squared();
        }

        Ok((res / n as f32).sqrt())
    }
}

/// Trait for analysis requiring positions and masses
pub trait MeasureMasses: PosProvider + MassesProvider + LenProvider {
    fn center_of_mass(&self) -> Result<Pos, MeasureError> {
        let mut cm = Vector3f::zeros();
        let mut mass = 0.0;
        for (c, m) in zip(self.iter_pos(), self.iter_masses()) {
            cm += c.coords * m;
            mass += m;
        }

        if mass == 0.0 {
            Err(MeasureError::ZeroMass)
        } else {
            Ok(Pos::from(cm / mass))
        }
    }

    fn fit_transform(
        sel1: &Self,
        sel2: &Self,
    ) -> Result<nalgebra::IsometryMatrix3<f32>, MeasureError> {
        let cm1 = sel1.center_of_mass()?;
        let cm2 = sel2.center_of_mass()?;

        //let rot = rot_transform_matrix(coords1.iter(), coords2.iter(), masses.iter());
        let rot = rot_transform(
            sel1.iter_pos().map(|p| *p - cm1),
            sel2.iter_pos().map(|p| *p - cm2),
            sel1.iter_masses(),
        )?;

        Ok(Translation3::from(cm2) * rot * Translation3::from(-cm1))
    }

    fn fit_transform_at_origin(
        sel1: &Self,
        sel2: &Self,
    ) -> Result<nalgebra::IsometryMatrix3<f32>, MeasureError> {
        let rot = rot_transform(
            sel1.iter_pos().map(|p| p.coords),
            sel2.iter_pos().map(|p| p.coords),
            sel1.iter_masses(),
        )?;
        Ok(nalgebra::convert(rot))
    }

    fn rmsd_mw(sel1: &Self, sel2: &Self) -> Result<f32, MeasureError> {
        let mut res = 0.0;
        let mut m_tot = 0.0;
        let iter1 = sel1.iter_pos();
        let iter2 = sel2.iter_pos();

        if sel1.len() != sel2.len() {
            return Err(MeasureError::Sizes(sel1.len(), sel2.len()));
        }

        for (p1, p2, m) in izip!(iter1, iter2, sel1.iter_masses()) {
            res += (p2 - p1).norm_squared() * m;
            m_tot += m;
        }

        if m_tot == 0.0 {
            Err(MeasureError::ZeroMass)
        } else {
            Ok((res / m_tot).sqrt())
        }
    }

    fn gyration(&self) -> Result<f32, MeasureError> {
        let c = self.center_of_mass()?;
        Ok(do_gyration(self.iter_pos().map(|pos| pos - c), self.iter_masses()))
    }

    fn inertia(&self) -> Result<(Vector3f, Matrix3f), MeasureError> {
        let c = self.center_of_mass()?;
        Ok(do_inertia(self.iter_pos().map(|pos| pos - c), self.iter_masses()))
    }

    fn principal_transform(&self) -> Result<IsometryMatrix3<f32>, MeasureError> {
        let c = self.center_of_mass()?;
        let (_, axes) = do_inertia(self.iter_pos().map(|pos| pos - c), self.iter_masses());
        Ok(do_principal_transform(axes, c.coords))
    }
}

/// Trait for analysis requiring positions, masses and pbc
pub trait MeasurePeriodic: PosProvider + MassesProvider + BoxProvider + LenProvider {
    fn center_of_mass_pbc(&self) -> Result<Pos, MeasureError> {
        let b = self.get_box().ok_or_else(|| MeasureError::NoPbc)?;
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
            Err(MeasureError::ZeroMass)
        } else {
            Ok(Pos::from(cm / mass))
        }
    }

    fn center_of_mass_pbc_dims(&self, dims: PbcDims) -> Result<Pos, MeasureError> {
        let b = self.get_box().ok_or_else(|| MeasureError::NoPbc)?;
        let mut pos_iter = self.iter_pos();
        let mut mass_iter = self.iter_masses();

        let mut mass = mass_iter.next().unwrap();
        let p0 = pos_iter.next().unwrap();
        let mut cm = p0.coords;

        for (c, m) in zip(pos_iter, mass_iter) {
            let im = b.closest_image_dims(c, p0, dims).coords;
            cm += im * m;
            mass += m;
        }

        if mass == 0.0 {
            Err(MeasureError::ZeroMass)
        } else {
            Ok(Pos::from(cm / mass))
        }
    }

    fn center_of_geometry_pbc(&self) -> Result<Pos, MeasureError> {
        let b = self.get_box().ok_or_else(|| MeasureError::NoPbc)?;
        let mut pos_iter = self.iter_pos();

        let p0 = pos_iter.next().unwrap();
        let mut cm = p0.coords;

        for c in pos_iter {
            cm += b.closest_image(c, p0).coords;
        }

        Ok(Pos::from(cm / self.len() as f32))
    }

    fn center_of_geometry_pbc_dims(&self, dims: PbcDims) -> Result<Pos, MeasureError> {
        let b = self.get_box().ok_or_else(|| MeasureError::NoPbc)?;
        let mut pos_iter = self.iter_pos();

        let p0 = pos_iter.next().unwrap();
        let mut cm = p0.coords;

        for c in pos_iter {
            cm += b.closest_image_dims(c, p0, dims).coords;
        }

        Ok(Pos::from(cm / self.len() as f32))
    }

    fn gyration_pbc(&self) -> Result<f32, MeasureError> {
        let b = self.get_box().ok_or_else(|| MeasureError::NoPbc)?;
        let c = self.center_of_mass_pbc()?;
        Ok(do_gyration(
            self.iter_pos().map(|pos| b.shortest_vector(&(pos - c))),
            self.iter_masses(),
        ))
    }

    fn inertia_pbc(&self) -> Result<(Vector3f, Matrix3f), MeasureError> {
        let b = self.get_box().ok_or_else(|| MeasureError::NoPbc)?;
        let c = self.center_of_mass_pbc()?;
        Ok(do_inertia(
            self.iter_pos().map(|pos| b.shortest_vector(&(pos - c))),
            self.iter_masses(),
        ))
    }

    fn principal_transform_pbc(&self) -> Result<IsometryMatrix3<f32>, MeasureError> {
        let b = self.get_box().ok_or_else(|| MeasureError::NoPbc)?;
        let c = self.center_of_mass_pbc()?;
        let (_, axes) = do_inertia(
            self.iter_pos().map(|pos| b.shortest_vector(&(pos - c))),
            self.iter_masses(),
        );
        Ok(do_principal_transform(axes, c.coords))
    }
}

fn do_gyration(
    dists: impl Iterator<Item = Vector3f>,
    masses: impl Iterator<Item = f32>,
) -> f32 {
    let mut sd = 0.0;
    let mut sm = 0.0;
    for (d, m) in zip(dists, masses) {
        sd += d.norm_squared() * m;
        sm += m;
    }
    // Don't need to test for zero mass since this is already done in center_of_mass()
    (sd / sm).sqrt()
}

fn do_inertia(
    dists: impl Iterator<Item = Vector3f>,
    masses: impl Iterator<Item = f32>,
) -> (Vector3f, Matrix3f) {
    let mut tens = Matrix3f::zeros();

    for (d, m) in zip(dists, masses) {
        tens[(0, 0)] += m * (d.y * d.y + d.z * d.z);
        tens[(1, 1)] += m * (d.x * d.x + d.z * d.z);
        tens[(2, 2)] += m * (d.x * d.x + d.y * d.y);
        tens[(0, 1)] -= m * d.x * d.y;
        tens[(0, 2)] -= m * d.x * d.z;
        tens[(1, 2)] -= m * d.y * d.z;
    }
    tens[(1, 0)] = tens[(0, 1)];
    tens[(2, 0)] = tens[(0, 2)];
    tens[(2, 1)] = tens[(1, 2)];

    // Compute axes and moments
    let eig = SymmetricEigen::new(tens);
    // Eigenvectors are NOT sorted, so sort them manually
    let mut s = eig
        .eigenvalues
        .into_iter()
        .enumerate()
        .collect::<Vec<(_, _)>>();
    s.sort_unstable_by(|a, b| a.1.partial_cmp(b.1).unwrap());
    // Create reordered moments and axes
    let moments = Vector3f::new(*s[0].1, *s[1].1, *s[2].1);

    let col0 = eig.eigenvectors.column(s[0].0).normalize();
    let col1 = eig.eigenvectors.column(s[1].0).normalize();
    let col2 = col0.cross(&col1); // Ensure right-handed system

    let axes = Matrix3f::from_columns(&[col0, col1, col2]);

    (moments, axes)
}

// Straighforward implementation of Kabsch algorithm
pub fn rot_transform(
    pos1: impl Iterator<Item = Vector3f>,
    pos2: impl Iterator<Item = Vector3f>,
    masses: impl Iterator<Item = f32>,
) -> Result<Rotation3<f32>, MeasureError> {
    //Calculate the covariance matrix
    let mut cov = Matrix3f::zeros();

    for (p1, p2, m) in izip!(pos1, pos2, masses) {
        cov += p2 * p1.transpose() * m;
    }

    // Perform Singular Value Decomposition (SVD) on the covariance matrix
    let svd = nalgebra::SVD::new(cov, true, true);
    let u = svd.u.ok_or_else(|| MeasureError::Svd)?;
    let v_t = svd.v_t.ok_or_else(|| MeasureError::Svd)?;

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
    Ok(Rotation3::from_matrix_unchecked(u * d_matrix * v_t))
}

// Principal transform algorithms
fn do_principal_transform(mut axes: Matrix3f, cm: Vector3f) -> IsometryMatrix3<f32> {
    axes.try_inverse_mut();
    Translation3::from(cm) * Rotation3::from_matrix_unchecked(axes) * Translation3::from(-cm)
}
