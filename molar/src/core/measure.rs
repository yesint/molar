use super::providers::*;
use super::{Matrix3f, PbcDims, Pos, Vector3f};
use itertools::izip;
use nalgebra::{IsometryMatrix3, Rotation3, SymmetricEigen, Translation3};
use num_traits::Bounded;
use std::f32::consts::PI;
use std::iter::zip;
use thiserror::Error;

//==============================================================
// Traits for measuring (immutable access)
//==============================================================

/// Errors that can occur during measurements
#[derive(Error, Debug)]
pub enum MeasureError {
    /// Mismatch in sizes between two selections
    #[error("invalid data sizes: {0} and {1}")]
    Sizes(usize, usize),

    /// Total mass of the selection is zero
    #[error("zero mass")]
    ZeroMass,

    /// Singular Value Decomposition failed
    #[error("SVD failed")]
    Svd,

    /// Operation requires periodic boundary conditions but none are defined
    #[error("no periodic box")]
    NoPbc,

    /// Cannot unwrap coordinates due to disjoint pieces
    #[error("can't unwrap disjoint pieces")]
    Disjoint,

    #[error("lipid order error")]
    LipidOrder(#[from] LipidOrderError)
}

/// Trait for analysis requiring only positions
pub trait MeasurePos: PosProvider + LenProvider {
    /// Returns the minimum and maximum coordinates across all dimensions
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

    /// Calculates the geometric center (centroid) of all positions
    fn center_of_geometry(&self) -> Pos {
        let iter = self.iter_pos();
        let n = self.len();
        let mut cog = Vector3f::zeros();
        for c in iter {
            cog += c.coords;
        }
        Pos::from(cog / n as f32)
    }

    /// Calculates the Root Mean Square Deviation between two selections
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
    /// Calculates the center of mass of the selection
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

    /// Computes the transformation that best fits sel1 onto sel2
    fn fit_transform(
        sel1: &Self,
        sel2: &impl MeasureMasses,
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

    /// Like fit_transform but assumes both selections are centered at origin
    fn fit_transform_at_origin(
        sel1: &Self,
        sel2: &impl MeasureMasses,
    ) -> Result<nalgebra::IsometryMatrix3<f32>, MeasureError> {
        let rot = rot_transform(
            sel1.iter_pos().map(|p| p.coords),
            sel2.iter_pos().map(|p| p.coords),
            sel1.iter_masses(),
        )?;
        Ok(nalgebra::convert(rot))
    }

    /// Calculates the mass-weighted Root Mean Square Deviation between two selections
    fn rmsd_mw(sel1: &Self, sel2: &impl MeasureMasses) -> Result<f32, MeasureError> {
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

    /// Calculates the radius of gyration of the selection
    fn gyration(&self) -> Result<f32, MeasureError> {
        let c = self.center_of_mass()?;
        Ok(do_gyration(
            self.iter_pos().map(|pos| pos - c),
            self.iter_masses(),
        ))
    }

    /// Calculates the moments and axes of inertia
    fn inertia(&self) -> Result<(Vector3f, Matrix3f), MeasureError> {
        let c = self.center_of_mass()?;
        Ok(do_inertia(
            self.iter_pos().map(|pos| pos - c),
            self.iter_masses(),
        ))
    }

    /// Computes transformation to principal axes of inertia
    fn principal_transform(&self) -> Result<IsometryMatrix3<f32>, MeasureError> {
        let c = self.center_of_mass()?;
        let (_, axes) = do_inertia(self.iter_pos().map(|pos| pos - c), self.iter_masses());
        Ok(do_principal_transform(axes, c.coords))
    }
}

/// Trait for analysis requiring positions, masses and periodic boundary conditions
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

/// Helper function to calculate radius of gyration from a set of distances and masses
fn do_gyration(dists: impl Iterator<Item = Vector3f>, masses: impl Iterator<Item = f32>) -> f32 {
    let mut sd = 0.0;
    let mut sm = 0.0;
    for (d, m) in zip(dists, masses) {
        sd += d.norm_squared() * m;
        sm += m;
    }
    // Don't need to test for zero mass since this is already done in center_of_mass()
    (sd / sm).sqrt()
}

/// Helper function to calculate inertia tensor from a set of distances and masses
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

/// Implements the Kabsch algorithm for finding optimal rotation between two sets of points
fn rot_transform(
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

/// Creates a transformation that aligns a structure with its principal axes
fn do_principal_transform(mut axes: Matrix3f, cm: Vector3f) -> IsometryMatrix3<f32> {
    axes.try_inverse_mut();
    Translation3::from(cm) * Rotation3::from_matrix_unchecked(axes) * Translation3::from(-cm)
}

/// Trait for modification requiring random access positions and pbc
pub trait MeasureRandomAccess: RandomPosProvider {
    /// Computes order parameter of the lipid tail. Each position in [Self] is
    /// supposed to represent a carbon atom of a single lipid tail. The size of the array
    /// of normals is either `1` or `N-2`, where `N` is the number of position in [Self].
    /// In the first case the same normal is used for all bonds, in the second case each
    /// bond in the range `1..N-1` has its own normal.
    /// `bonds_orders` should contain either 1 or 2 for all `N-1` bonds.
    /// Formulas are taken from:
    /// - https://pubs.acs.org/doi/10.1021/acs.jctc.7b00643
    /// - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3882000/
    fn lipid_tail_order(
        &self,
        order_type: OrderType,
        normals: &Vec<Vector3f>,
        bond_orders: &Vec<u8>,
    ) -> Result<Vec<f32>, LipidOrderError> {
        //atoms:  0 - 1 - 2 - 3 = 4 - 5 - 6
        //bonds:    0   1   2   3   4   5
        //normals:  0   1   2   3   4   5

        // Size check
        if self.num_coords() < 3 {
            return Err(LipidOrderError::TailTooShort(self.num_coords()));
        }
        
        if normals.len() != 1 && normals.len() != self.num_coords() - 2 {
            return Err(LipidOrderError::NormalsCount(self.num_coords(), self.num_coords()-2));
        }

        if bond_orders.len() != self.num_coords() - 1 {
            return Err(LipidOrderError::BondOrderCount(self.num_coords(), self.num_coords()-1));
        }

        let mut order = vec![0.0; self.num_coords() - 2];
        if order_type == OrderType::Sz {
            // Iterate over atoms
            for at in 1..self.num_coords() - 1 {
                // Vector from at+1 to at-1
                let v = unsafe { self.nth_pos_unchecked(at + 1) - self.nth_pos_unchecked(at - 1) };
                // Normal
                let normal = if normals.len() == 1 {
                    &normals[0]
                } else {
                    &normals[at - 1] // indexing starts from zero
                };
                // Angle
                let ang = v.angle(&normal);
                order[at - 1] = 1.5 * ang.cos().powi(2) - 0.5;
            }
        } else {
            // Compute deuterium order
            // We iterate over bonds and treat differently single and double bonds
            for i in 0..self.num_coords() - 2 {
                if bond_orders[i] == 1 {
                    // Single bond between atoms i:i+1
                    // If next bond is also single, compute order for atom i+1
                    if bond_orders[i + 1] == 1 {
                        // Compute single order for atom i+1
                        /*
                         * C(i)[1]
                         *  \
                         *   C(i+1)[2]---H1,2
                         *  /
                         * C(i+2)[3]
                         */

                        let p1 = unsafe { self.nth_pos_unchecked(i) };
                        let p2 = unsafe { self.nth_pos_unchecked(i + 1) };
                        let p3 = unsafe { self.nth_pos_unchecked(i + 2) };

                        let local_z = (p3 - p1).normalize();
                        let local_x = ((p1 - p2).cross(&(p3 - p2))).normalize();
                        let local_y = local_x.cross(&local_z);

                        // Normal
                        let n = if normals.len() == 1 {
                            &normals[0]
                        } else {
                            &normals[i] // indexing starts from zero, (i+1)-1 = i
                        };

                        let ang_x = local_x.angle(n);
                        let ang_y = local_y.angle(n);
                        let sxx = 0.5 * (3.0 * ang_x.cos().powi(2) - 1.0);
                        let syy = 0.5 * (3.0 * ang_y.cos().powi(2) - 1.0);
                        // Instantaneous order
                        // Index in order array is c_atom-1: (i+1)-1=i
                        order[i] = -(2.0 * sxx + syy) / 3.0;
                        // For single bonds there is no difference between ideal
                        // and corrected versions of Scd
                    }
                    // If next bond is double we do nothing and wait for next iteration
                } else {
                    // Double bond between atoms i:i+1
                    // Compute double order for atoms i and i+1
                    /*
                     * C(i-1)[1]
                     *  \
                     *   C(i)[2]----H1
                     *   ||
                     *   C(i+1)[3]--H2
                     *  /
                     * C(i+2)[4]
                     *
                     * a1 = 0.5*(pi-ang(1,2,3))
                     * a2 = 0.5*(pi-ang(2,3,4))
                     */

                    let c1 = i - 1;
                    let c2 = i;
                    let c3 = i + 1;
                    let c4 = i + 2;

                    let p1 = unsafe { self.nth_pos_unchecked(c1) };
                    let p2 = unsafe { self.nth_pos_unchecked(c2) };
                    let p3 = unsafe { self.nth_pos_unchecked(c3) };
                    let p4 = unsafe { self.nth_pos_unchecked(c4) };

                    let a1 = 0.5 * (PI - (p1 - p2).angle(&(p3 - p2)));
                    let a2 = 0.5 * (PI - (p2 - p3).angle(&(p4 - p3)));

                    // For atom i
                    let local_z = (p3 - p2).normalize();
                    let local_x = ((p1 - p2).cross(&local_z)).normalize();
                    let local_y = local_x.cross(&local_z);

                    let n1 = if normals.len() == 1 {
                        &normals[0]
                    } else {
                        &normals[i] // indexing starts from zero, (i+1)-1 = i
                    };

                    let ang_y = local_y.angle(&n1);
                    let ang_z = local_z.angle(&n1);
                    let szz = 0.5 * (3.0 * ang_z.cos().powi(2) - 1.0);
                    let syy = 0.5 * (3.0 * ang_y.cos().powi(2) - 1.0);
                    let syz = 1.5 * ang_y.cos() * ang_z.cos();
                    // Index in order array is c_atom-1: i-1
                    if order_type == OrderType::ScdCorr {
                        order[i - 1] = -(a1.cos().powi(2) * syy + a1.sin().powi(2) * szz
                            - 2.0 * a1.cos() * a1.sin() * syz);
                    } else {
                        // SCD
                        order[i - 1] = -(szz / 4.0 + 3.0 * syy / 4.0 - 3.0_f32.sqrt() * syz / 2.0);
                    }

                    // For atom i+1
                    // same local_z is used
                    let local_x = ((p3 - p4).cross(&local_z)).normalize();
                    let local_y = local_x.cross(&local_z);

                    let n2 = if normals.len() == 1 {
                        &normals[0]
                    } else {
                        &normals[i + 1] // indexing starts from zero, (i+2)-1 = i+1
                    };

                    let ang_y = local_y.angle(n2);
                    let ang_z = local_z.angle(n2);
                    let szz = 0.5 * (3.0 * ang_z.cos().powi(2) - 1.0);
                    let syy = 0.5 * (3.0 * ang_y.cos().powi(2) - 1.0);
                    let syz = 1.5 * ang_y.cos() * ang_z.cos();
                    // Index in order array is c_atom-1: (i+1)-1=i
                    if order_type == OrderType::ScdCorr {
                        order[i] = -(a2.cos().powi(2) * syy
                            + a2.sin().powi(2) * szz
                            + 2.0 * a2.cos() * a2.sin() * syz);
                    } else {
                        // SCD
                        order[i] = -(szz / 4.0 + 3.0 * syy / 4.0 + 3.0_f32.sqrt() * syz / 2.0);
                    }
                } // if single/double
            } // for bonds
        }
        Ok(order)
    }
}

/// Type of order parameter calculation
#[derive(PartialEq, Debug)]
pub enum OrderType {
    // Sz order parameter identical to gromacs -szonly option
    Sz,
    // Deuterium order parameter computed for ideal H positions for double bonds
    Scd,
    // Deuterium order parameter computed for corrected H positions for double bonds
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3882000/
    ScdCorr,
}

#[derive(Error,Debug)]
pub enum LipidOrderError {
    #[error("for {0} tail carbons # of normals should be 1 or {1}")]
    NormalsCount(usize,usize),
    
    #[error("for {0} tail carbons # of bond orders should be {1}")]
    BondOrderCount(usize,usize),

    #[error("tail should have at least 3 carbons, not {0}")]
    TailTooShort(usize),
}