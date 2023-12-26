use core::f32::consts::SQRT_2;

use anyhow::{bail, Result};
use num_traits::Zero;
use num_traits::Bounded;
use crate::distance_search::search::{SearchConnectivity, DistanceSearcherSingle};
use super::Matrix3f;
use super::{AtomIterator, AtomMutIterator, ParticleMut, ParticleMutIterator, PbcDims, PeriodicBox, PosMutIterator, ParticleIterator, Pos, Vector3f, PosIterator};
 
// Trait that provides periodic box information
pub trait MeasureBox {
    fn get_box(&self) -> Result<&PeriodicBox>;
}

/// Trait for measuring various properties that requires only
/// the iterator of positions.
pub trait MeasurePos {
    fn iter_pos(&self) -> impl PosIterator<'_>;

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

pub trait MeasureAtoms {
    fn iter_atoms(&self) -> impl AtomIterator<'_>;
}
/// Trait for measuring various properties that requires only
/// the iterator of particles. User types should 
/// implement `iter`
pub trait MeasureParticles: MeasurePos {
    fn iter_particles(&self) -> impl ParticleIterator<'_>;

    fn center_of_mass(&self) -> Result<Pos> {
        let c = self.iter_particles().fold(
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

/// The trait for measuring properties that requires
/// a periodic box information.
pub trait MeasurePeriodic: MeasureParticles + MeasureBox {
    fn center_of_mass_pbc(&self) -> Result<Pos> {
        let b = self.get_box()?;
        let mut iter = self.iter_particles();
        let p0 = iter.next().unwrap().pos;
        let c = iter.fold(
            (Vector3f::zero(),0.0), 
            |acc, p| {
                let im = b.closest_image(p.pos,p0).coords;
                (acc.0+im*p.atom.mass, acc.1+p.atom.mass)
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
pub trait ModifyParticles {
    fn iter_particles_mut(&mut self) -> impl ParticleMutIterator<'_>;
    
    fn iter_particles(&mut self) -> impl ParticleIterator<'_>{
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
}

/// The trait for modifying the particles that requires
/// the periodic box.
pub trait ModifyPeriodic: ModifyParticles + MeasureBox {
    fn unwrap_simple_dim(&mut self, dims: PbcDims) -> Result<()> {
        let b = self.get_box()?.clone();
        let mut iter = self.iter_pos_mut();
        if iter.len()>0 {
            let p0 = iter.next().unwrap();
            for p in iter {
                *p = b.closest_image_dims(p, p0, &dims);
            }
        }
        Ok(())
    }

    fn unwrap_simple(&mut self) -> Result<()> {
        self.unwrap_simple_dim([true,true,true])
    }
}

pub trait ModifyRandomAccess: ModifyPeriodic {
    fn nth_particle_mut(&mut self, i: usize) -> ParticleMut;
    fn nth_pos_mut(&mut self, i: usize) -> &mut Pos;

    fn nth_pos(&mut self, i: usize) -> &Pos {
        self.nth_pos_mut(i)
    }

    fn unwrap_connectivity(&mut self, cutoff: f32) -> Result<()> {
        self.unwrap_connectivity_dim(cutoff, &[true,true,true])
    }
    
    fn unwrap_connectivity_dim(&mut self, cutoff: f32, dims: &PbcDims) -> Result<()> {
        let b = self.get_box()?.to_owned();
        let conn: SearchConnectivity = DistanceSearcherSingle::new_periodic(
            cutoff,
            self.iter_particles().map(|p| (p.id, p.pos)),
            &b,
            &dims
        ).search();

        // used atoms
        let mut used = vec![false;conn.len()];
        // Centers to unwrap
        let mut todo = Vec::<usize>::with_capacity(conn.len()/2);
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
            bail!("Selection is not compact for cutoff={}",cutoff)
        }

        Ok(())
    }
}

/// Computes a rotational part of the fit transform 
pub fn rot_transform_matrix<'a>(
    sel1: impl ParticleIterator<'a>, 
    sel2: impl ParticleIterator<'a>
) -> nalgebra::Matrix3<f32> 
{
    let n = sel1.len();

    //Calculate the matrix U
    let mut u = Matrix3f::zeros();

    for (p1,p2) in std::iter::zip(sel1,sel2) {
        u += p1.pos.coords * p2.pos.coords.transpose() * p1.atom.mass;
    }

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
    omega.fixed_view_mut::<3,3>(0,3).copy_from(&u.transpose());
    omega.fixed_view_mut::<3,3>(3,0).copy_from(&u);

    //Finding eigenvalues of omega
    // Lapack binding produces sorted eigenvectors, while native nalgebra do not!
    let eig = nalgebra_lapack::SymmetricEigen::new(omega);
    let om = eig.eigenvectors;
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
    let vh0 = om.fixed_view::<3,1>(0, 5).transpose() * SQRT_2;
    let vh1 = om.fixed_view::<3,1>(0, 4).transpose() * SQRT_2;
    let vh2 = vh0.cross(&vh1);
    let vh = Matrix3f::from_rows(&[vh0,vh1,vh2]);

    let vk0 = om.fixed_view::<3,1>(3, 5).transpose() * SQRT_2;
    let vk1 = om.fixed_view::<3,1>(3, 4).transpose() * SQRT_2;
    let vk2 = vk0.cross(&vk1);
    let vk = Matrix3f::from_rows(&[vk0,vk1,vk2]);

    /* Determine rotational part */
    let mut rot = Matrix3f::zeros();
    for r in 0..3 {
        for c in 0..3 {
            rot[(c,r)] = vk[(0,r)]*vh[(0,c)] + vk[(1,r)]*vh[(1,c)] + vk[(2,r)]*vh[(2,c)];
        }
    }
    
    rot
}
