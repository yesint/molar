use std::iter::zip;

use super::Matrix3f;
use super::{ PbcDims, Pos, Vector3f};
use crate::distance_search::search::{DistanceSearcherSingle, SearchConnectivity};
use anyhow::{bail, Result};
use itertools::izip;
use nalgebra::Rotation3;
use nalgebra::Unit;
use nalgebra::SVD;
use num_traits::Bounded;
use num_traits::Zero;
use super::providers::*;

//---------------------------------------------------
// Free functions for computing properties that
// acceps a needed type of data provider as argument
//---------------------------------------------------
pub fn min_max(dp: &impl PosProvider) -> (Pos, Pos) {
    let mut lower = Pos::max_value();
    let mut upper = Pos::min_value();
    for p in dp.iter_pos() {
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

pub fn center_of_geometry(dp: &impl PosProvider) -> Pos {
    let iter = dp.iter_pos();
    let n = iter.len();
    let mut cog = Vector3f::zero();
    for c in iter {
        cog += c.coords;
    }
    Pos::from(cog / n as f32)
}

pub fn center_of_mass(dp: &(impl PosProvider+MassesProvider)) -> Result<Pos> {
    let mut cm = Vector3f::zero();
    let mut mass = 0.0;
    for (c, m) in zip(dp.iter_pos(), dp.iter_masses()) {
        cm += c.coords * m;
        mass += m;
    }

    if mass == 0.0 {
        bail!("Zero mass in COM!")
    } else {
        Ok(Pos::from(cm / mass))
    }
}

pub fn center_of_mass_pbc(dp: &(impl PosProvider+MassesProvider+BoxProvider)) -> Result<Pos> {
    let b = dp.get_box()?;
    let mut pos_iter = dp.iter_pos();
    let mut mass_iter = dp.iter_masses();
    
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

//-------------------------------------------------------
// Free functions for modifying atoms
//-------------------------------------------------------

pub fn assign_resindex(dp: &mut impl AtomsMutProvider) {
    let mut resindex = 0usize;
    let mut at_iter = dp.iter_atoms_mut();
    if at_iter.len()>1 {
        let at0 = at_iter.next().unwrap();
        let mut cur_resid = at0.resid;
        at0.resindex = resindex;
        for at in at_iter {
            if at.resid != cur_resid {
                cur_resid = at.resid;
                resindex += 1;
            }
            at.resindex = resindex;
        }
    }
}

//-------------------------------------------------------
// Free functions for modifying that
// acceps a needed type of mut data provider as argument
//-------------------------------------------------------
pub fn translate(dp: &mut impl PosMutProvider, shift: Vector3f) {
    for el in dp.iter_pos_mut() {
        *el += shift;
    }
}

pub fn rotate(dp: &mut impl PosMutProvider, ax: &Unit<Vector3f>, ang: f32) {
    let tr = Rotation3::<f32>::from_axis_angle(ax, ang);
    for p in dp.iter_pos_mut() {
        p.coords = tr * p.coords;
    }
}

pub fn apply_transform(dp: &mut impl PosMutProvider, tr: &nalgebra::IsometryMatrix3<f32>) {
    for p in dp.iter_pos_mut() {
        *p = tr * (*p);
    }
}

pub fn unwrap_simple_dim(dp: &mut (impl PosMutProvider+BoxProvider), dims: PbcDims) -> Result<()> {
    let b = dp.get_box()?.to_owned();
    let mut iter = dp.iter_pos_mut();
    if iter.len() > 0 {
        let p0 = iter.next().unwrap();
        for p in iter {
            *p = b.closest_image_dims(p, p0, &dims);
        }
    }
    Ok(())
}

pub fn unwrap_connectivity_dim(dp: &mut (impl PosMutProvider + PosProvider + BoxProvider + RandomPosMutProvider), cutoff: f32, dims: &PbcDims) -> Result<()> {
    let b = dp.get_box()?.to_owned();
    let conn: SearchConnectivity = DistanceSearcherSingle::new_periodic(
        cutoff,
        dp.iter_pos().enumerate(),
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
        let p0 = dp.nth_pos(c).to_owned();
        // Iterate over connected points
        for ind in &conn[c] {
            // Unwrap this point if it is not used yet
            if !used[*ind] {
                let p = dp.nth_pos_mut(*ind);
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

//---------------------------------------------------
// Free functions for RMSD and fitting
//---------------------------------------------------
// Straighforward implementation of Kabsch algorithm
pub fn rot_transform (
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

pub fn fit_transform (
    dp1: &(impl PosProvider+MassesProvider),
    dp2: &(impl PosProvider+MassesProvider),
) -> Result<nalgebra::IsometryMatrix3<f32>> 
{
    let cm1 = center_of_mass(dp1)?;
    let cm2 = center_of_mass(dp2)?;

    //let rot = rot_transform_matrix(coords1.iter(), coords2.iter(), masses.iter());
    let rot = rot_transform(
        dp1.iter_pos().map(|p| *p-cm1),
        dp2.iter_pos().map(|p| *p-cm2),
        dp1.iter_masses()
    );

    Ok(nalgebra::Translation3::from(cm2) * rot * nalgebra::Translation3::from(-cm1))
}

// Version for selection with CM alredy at zero
pub fn fit_transform_at_origin (
    dp1: &(impl PosProvider+MassesProvider),
    dp2: &(impl PosProvider+MassesProvider),
) -> Result<nalgebra::IsometryMatrix3<f32>> 
{
    let rot = rot_transform(
        dp1.iter_pos().map(|p| p.coords),
        dp2.iter_pos().map(|p| p.coords),
        dp1.iter_masses()
    );
    Ok(nalgebra::convert(rot))
}

/// Mass-weighted RMSD
pub fn rmsd_mw (
    dp1: &(impl PosProvider+MassesProvider),
    dp2: &(impl PosProvider+MassesProvider),
) -> Result<f32> 
{
    let mut res = 0.0;
    let mut m_tot = 0.0;
    let iter1 = dp1.iter_pos();
    let iter2 = dp2.iter_pos();

    if iter1.len() != iter2.len() {
        bail!("Different sizes in rmsd_mw: {} and {}",iter1.len(),iter2.len());
    }

    for (p1,p2,m) in izip!(iter1,iter2,dp1.iter_masses()){
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
pub fn rmsd (
    dp1: &impl PosProvider,
    dp2: &impl PosProvider,
) -> Result<f32> 
{
    let mut res = 0.0;
    let iter1 = dp1.iter_pos();
    let iter2 = dp2.iter_pos();

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

