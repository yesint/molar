use crate::distance_search::DistanceSearcherSingle;
use crate::distance_search::SearchConnectivity;

use super::providers::*;
use super::PbcDims;
use super::Vector3f;
use super::PBC_FULL;
use anyhow::{anyhow, bail, Result};
use nalgebra::Rotation3;
use nalgebra::Unit;

//==============================================================
// Traits for modification (mutable access)
//==============================================================

/// Trait for modification requiring only positions
pub trait ModifyPos: PosMutProvider + PosProvider {
    fn translate(&self, shift: &Vector3f) {
        for el in self.iter_pos_mut() {
            *el += shift;
        }
    }

    fn rotate(&self, ax: &Unit<Vector3f>, ang: f32) {
        let tr = Rotation3::<f32>::from_axis_angle(ax, ang);
        for p in self.iter_pos_mut() {
            p.coords = tr * p.coords;
        }
    }

    fn apply_transform(&self, tr: &nalgebra::IsometryMatrix3<f32>) {
        for p in self.iter_pos_mut() {
            *p = tr * (*p);
        }
    }
}

/// Trait for modification requiring positions and pbc
pub trait ModifyPeriodic: PosMutProvider + BoxProvider {
    fn unwrap_simple_dim(&self, dims: PbcDims) -> Result<()> {
        let b = self
            .get_box()
            .ok_or_else(|| anyhow!("No periodicity!"))?
            .to_owned();
        let mut iter = self.iter_pos_mut();
        if iter.len() > 0 {
            let p0 = iter.next().unwrap();
            for p in iter {
                *p = b.closest_image_dims(p, p0, &dims);
            }
        }
        Ok(())
    }

    fn unwrap_simple(&self) -> Result<()> {
        self.unwrap_simple_dim(PBC_FULL)
    }
}

/// Trait for modification requiring random access positions and pbc
pub trait ModifyRandomAccess:
    PosMutProvider + PosProvider + BoxProvider + RandomPosMutProvider
{
    fn unwrap_connectivity(&self, cutoff: f32) -> Result<()> {
        self.unwrap_connectivity_dim(cutoff, &PBC_FULL)
    }

    fn unwrap_connectivity_dim(&self, cutoff: f32, dims: &PbcDims) -> Result<()> {
        let b = self
            .get_box()
            .ok_or_else(|| anyhow!("No periodicity!"))?
            .to_owned();
        let conn: SearchConnectivity =
            DistanceSearcherSingle::new_periodic(
                cutoff,
                self.iter_pos().cloned().enumerate(),
                &b,
                &dims
            ).search();

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
            let p0 = unsafe { self.nth_pos_unchecked(c) }.to_owned();
            // Iterate over connected points
            for ind in &conn[c] {
                // Unwrap this point if it is not used yet
                if !used[*ind] {
                    let p = unsafe { self.nth_pos_unchecked_mut(*ind) };
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

/// Trait for modification requiring atoms
pub trait ModifyAtoms: AtomsMutProvider {
    fn assign_resindex(&self) {
        let mut resindex = 0usize;
        let mut at_iter = self.iter_atoms_mut();
        if at_iter.len() > 1 {
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
}
