use core::f32::consts::SQRT_2;

use super::Matrix3f;
use super::Selection;
use super::State;
use super::Topology;
use super::{
    AtomIterator, AtomMutIterator, ParticleIterator, ParticleMut, ParticleMutIterator, PbcDims,
    PeriodicBox, Pos, PosIterator, PosMutIterator, Vector3f,
};
use crate::distance_search::search::{DistanceSearcherSingle, SearchConnectivity};
use anyhow::{bail, Result};
use nalgebra::Matrix3;
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
    fn iter_masses(&self) -> impl ExactSizeIterator<Item = &f32>;

    fn center_of_mass(&self) -> Result<Pos> {
        let mut cm = Vector3f::zero();
        let mut mass = 0.0;
        for (c, m) in std::iter::zip(self.iter_pos(), self.iter_masses()) {
            cm += c.coords * (*m);
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
        
        let mut mass = *mass_iter.next().unwrap();
        let p0 = pos_iter.next().unwrap();
        let mut cm = p0.coords;

        for (c, m) in std::iter::zip(pos_iter, mass_iter) {
            let im = b.closest_image(c, p0).coords;
            cm += im * (*m);
            mass += *m;
        }

        if mass == 0.0 {
            bail!("Zero mass in COM!")
        } else {
            Ok(Pos::from(cm / mass))
        }
    }
}

/// The trait for modifying the particles. User types should
/// implement `iter_mut`.
pub trait ModifyParticles {
    fn iter_particles_mut(&mut self) -> impl ParticleMutIterator<'_>;

    fn iter_particles(&mut self) -> impl ParticleIterator<'_> {
        self.iter_particles_mut().map(|p| p.into())
    }

    fn iter_pos_mut(&mut self) -> impl PosMutIterator<'_> {
        self.iter_particles_mut().map(|p| p.pos)
    }

    fn iter_atoms_mut(&mut self) -> impl AtomMutIterator<'_> {
        self.iter_particles_mut().map(|p| p.atom)
    }

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

/// The trait for modifying the particles that requires
/// the periodic box.
pub trait ModifyPeriodic: ModifyParticles + MeasureBox {
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
        self.unwrap_simple_dim([true, true, true])
    }
}

pub trait ModifyRandomAccess: ModifyPeriodic {
    fn nth_particle_mut(&mut self, i: usize) -> ParticleMut;
    fn nth_pos_mut(&mut self, i: usize) -> &mut Pos;

    fn nth_pos(&mut self, i: usize) -> &Pos {
        self.nth_pos_mut(i)
    }

    fn unwrap_connectivity(&mut self, cutoff: f32) -> Result<()> {
        self.unwrap_connectivity_dim(cutoff, &[true, true, true])
    }

    fn unwrap_connectivity_dim(&mut self, cutoff: f32, dims: &PbcDims) -> Result<()> {
        let b = self.get_box()?.to_owned();
        let conn: SearchConnectivity = DistanceSearcherSingle::new_periodic(
            cutoff,
            self.iter_particles().map(|p| (p.id, p.pos)),
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

// Rotational fitting using quaternions
// as described here: https://arxiv.org/pdf/physics/0506177.pdf, Page 3.
// Positions are assumed to be at the center of masseses
// This is SLOWER than Kabsch method!
#[allow(non_snake_case)]
pub fn rot_transform_quat<'a>(
    pos1: impl Iterator<Item = &'a Vector3f>,
    pos2: impl Iterator<Item = &'a Vector3f>,
    masses: impl Iterator<Item = &'a f32>,
) -> Unit<nalgebra::Quaternion<f32>> {
    // a = pos2+pos1
    // b = pos2-pos1
    // A =  0   -b1 -b2 -b3
    //      b1   0  -a3  a2
    //      b2   a3  0  -a1
    //      b3  -a2  a1  0
    let mut B = nalgebra::Matrix4::<f32>::zeros();

    for (p1, p2, m) in itertools::izip!(pos1, pos2, masses) {
        let a = p2 + p1;
        let b = p2 - p1;
        let A = nalgebra::Matrix4::<f32>::new(
            0.0, -b[0], -b[1], -b[2], b[0], 0.0, -a[2], a[1], b[1], a[2], 0.0, a[0], b[2], -a[1],
            a[0], 0.0,
        );

        B += A.transpose() * A * (*m);
    }
    // We can skip normalizing by mass because resulting quaternion
    // will be normalized anyway
    // B /= M;

    let eig = nalgebra_lapack::SymmetricEigen::new(B);
    let om = eig.eigenvectors;

    //println!("om =\n {}",om);
    //println!("val =\n {}", eig.eigenvalues);

    let q = nalgebra::Quaternion::<f32>::from_parts(om[(0, 0)], om.fixed_view::<3, 1>(1, 0));
    Unit::new_normalize(q)
}

// Implementation of Kabsch algorithm from Gromacs with 6x6 matrix
// This is less intuitive than direct implementation with SVD and currently requires
// workarrounds because of nalgebra weird behavior with eigenvectors signs.
pub fn rot_transform_gmx<'a>(
    pos1: impl Iterator<Item = &'a Vector3f>,
    pos2: impl Iterator<Item = &'a Vector3f>,
    masses: impl Iterator<Item = &'a f32>,
) -> Matrix3<f32> {
    //Calculate the matrix U
    let mut u = Matrix3f::zeros();

    for (p1, p2, m) in itertools::izip!(pos1, pos2, masses) {
        u += p1 * p2.transpose() * (*m);
    }

    //println!("u =\n {}",u);

    //Construct omega
    /*
     u= 1 4 7
        2 5 8
        3 6 9
    omega =
    0 0 0 1 2 3
    0 0 0 4 5 6
    0 0 0 7 8 9
    1 4 7 0 0 0
    2 5 8 0 0 0
    3 6 9 0 0 0
    */
    let mut omega = nalgebra::Matrix6::<f32>::zeros();
    omega.fixed_view_mut::<3, 3>(0, 3).copy_from(&u.transpose());
    omega.fixed_view_mut::<3, 3>(3, 0).copy_from(&u);

    //println!("omega =\n {}",omega);

    //Finding eigenvalues of omega
    let eig = nalgebra_lapack::SymmetricEigen::new(omega);
    let om = eig.eigenvectors;

    //println!("om =\n {}",om);
    //println!("val =\n {}", eig.eigenvalues);

    /*  Copy only the first two eigenvectors
        The eigenvectors are already sorted ascending by their eigenvalues!
    */

    /*
     i0 1 2 3 4 5
    j*-----------
    0|0 0 0 0 1 7
    1|0 0 0 0 2 8
    2|0 0 0 0 3 9
    3|0 0 0 0 4 10
    4|0 0 0 0 5 11
    5|0 0 0 0 6 12

    vh:
    7 8 9
    1 2 3
    ? ? ?

    vk:
    10 11 12
     4  5  6
     ?  ?  ?
    */

    // Calculate the last eigenvector as the cross-product of the first two.
    // This insures that the conformation is not mirrored and
    // prevents problems with completely flat reference structures.
    let vh0 = om.fixed_view::<3, 1>(0, 5) * SQRT_2;
    // For unknown reason nalgevra produces second-last eigenvector with wrong sign!
    let vh1 = -om.fixed_view::<3, 1>(0, 4) * SQRT_2;
    let vh2 = vh0.cross(&vh1);
    let vh = Matrix3f::from_columns(&[vh0, vh1, vh2]).transpose();

    let vk0 = om.fixed_view::<3, 1>(3, 5) * SQRT_2;
    // For unknown reason nalgevra produces second-last eigenvector with wrong sign!
    let vk1 = -om.fixed_view::<3, 1>(3, 4) * SQRT_2;
    let vk2 = vk0.cross(&vk1);
    let vk = Matrix3f::from_columns(&[vk0, vk1, vk2]).transpose();

    //println!("vh =\n {}",vh);
    //println!("vk =\n {}",vk);

    /* Determine rotational part */
    let mut rot = Matrix3f::zeros();
    for r in 0..3 {
        for c in 0..3 {
            rot[(c, r)] =
                vk[(0, r)] * vh[(0, c)] + vk[(1, r)] * vh[(1, c)] + vk[(2, r)] * vh[(2, c)];
        }
    }
    //rot

    //let rot = vk * vh.transpose();
    //println!("rot =\n {}",rot);
    rot
}

// Straighforward implementation of Kabsch algorithm (written by AI).
pub fn rot_transform_kabsch<'a>(
    pos1: impl Iterator<Item =  Vector3f>,
    pos2: impl Iterator<Item =  Vector3f>,
    masses: impl Iterator<Item =  &'a f32>,
) -> Matrix3<f32> {
    //Calculate the covariance matrix
    let mut cov = Matrix3f::zeros();

    for (p1, p2, m) in itertools::izip!(pos1, pos2, masses) {
        cov += p2 * p1.transpose() * *m;
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
    u * d_matrix * v_t
}

pub fn fit_transform_gmx<'a>(
    sel1: impl ParticleIterator<'a>,
    sel2: impl ParticleIterator<'a>,
) -> Result<nalgebra::Isometry3<f32>> {
    let (mut coords1, masses): (Vec<Vector3f>, Vec<f32>) =
        sel1.map(|p| (p.pos.coords, p.atom.mass)).unzip();
    let mut coords2: Vec<Vector3f> = sel2.map(|p| p.pos.coords).collect();

    let cm1 = center_of_mass(coords1.iter(), masses.iter())?;
    let cm2 = center_of_mass(coords2.iter(), masses.iter())?;

    //println!("cm1={}",cm1);
    //println!("cm2={}",cm2);

    for c in coords1.iter_mut() {
        *c -= cm1;
    }
    for c in coords2.iter_mut() {
        *c -= cm2;
    }
    //let rot = rot_transform_matrix(coords1.iter(), coords2.iter(), masses.iter());
    let rot = rot_transform_quat(coords1.iter(), coords2.iter(), masses.iter());

    Ok(nalgebra::Translation3::from(cm2) * rot * nalgebra::Translation3::from(-cm1))
}

pub fn fit_transform_matrix<'a>(
    sel1: impl ParticleIterator<'a>,
    sel2: impl ParticleIterator<'a>,
) -> Result<nalgebra::IsometryMatrix3<f32>> {
    let (mut coords1, masses): (Vec<Vector3f>, Vec<f32>) =
        sel1.map(|p| (p.pos.coords, p.atom.mass)).unzip();
    let mut coords2: Vec<Vector3f> = sel2.map(|p| p.pos.coords).collect();

    let cm1 = center_of_mass(coords1.iter(), masses.iter())?;
    let cm2 = center_of_mass(coords2.iter(), masses.iter())?;

    //println!("cm1={}",cm1);
    //println!("cm2={}",cm2);

    for c in coords1.iter_mut() {
        *c -= cm1;
    }
    for c in coords2.iter_mut() {
        *c -= cm2;
    }
    //let rot = rot_transform_matrix(coords1.iter(), coords2.iter(), masses.iter());
    let rot = rot_transform_gmx(coords1.iter(), coords2.iter(), masses.iter());

    Ok(nalgebra::Translation3::from(cm2)
        * Rotation3::from_matrix_unchecked(rot)
        * nalgebra::Translation3::from(-cm1))
}

pub fn fit_transform(
    sel1: impl MeasureMasses,
    sel2: impl MeasureMasses,
) -> Result<nalgebra::IsometryMatrix3<f32>> {
    let cm1 = sel1.center_of_mass()?;
    let cm2 = sel2.center_of_mass()?;

    //let rot = rot_transform_matrix(coords1.iter(), coords2.iter(), masses.iter());
    let rot = rot_transform_kabsch(
        sel1.iter_pos().map(|p| *p-cm1),
        sel2.iter_pos().map(|p| *p-cm2),
        sel1.iter_masses()
    );

    Ok(nalgebra::Translation3::from(cm2)
        * Rotation3::from_matrix_unchecked(rot)
        * nalgebra::Translation3::from(-cm1))
}


fn center_of_mass<'a>(
    coords: impl Iterator<Item = &'a Vector3f>,
    masses: impl Iterator<Item = &'a f32>,
) -> Result<Vector3f> {
    
    let mut cm = Vector3f::zero();
    let mut mass = 0.0;
    for (c, m) in std::iter::zip(coords, masses) {
        cm += c * (*m);
        mass += m;
    }

    if mass == 0.0 {
        bail!("Zero mass in COM!")
    } else {
        Ok(cm / mass)
    }
}
