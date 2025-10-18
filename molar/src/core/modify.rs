use super::providers::*;
use super::MeasureError;
use super::PbcDims;
use super::Vector3f;
use super::PBC_FULL;
use crate::prelude::*;
use nalgebra::Const;
use nalgebra::Rotation3;
use nalgebra::Unit;

//==============================================================
// Traits for modification (mutable access)
//==============================================================

/// Trait for modification requiring only positions
pub trait ModifyPos: PosIterMutProvider {
    //pub fn from_matrix<S>(matrix: nalgebra::Matrix<f32,Const<3>,Const<3>,S>) -> Result<Self, PeriodicBoxError>
    //where S: nalgebra::storage::Storage<f32, Const<3>, Const<3>>

    fn translate<S>(&mut self, shift: &nalgebra::Matrix<f32, Const<3>, Const<1>, S>)
    where
        S: nalgebra::storage::Storage<f32, Const<3>, Const<1>>,
    {
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

/// Trait for modification requiring positions and pbc
pub trait ModifyPeriodic: PosIterMutProvider + BoxMutProvider + LenProvider + BoxProvider {
    fn unwrap_simple_dim(&mut self, dims: PbcDims) -> Result<(), MeasureError> {
        let n = self.len();
        let b = self.require_box()?.to_owned();
        let mut iter = self.iter_pos_mut();
        if n > 0 {
            let p0 = iter.next().unwrap();
            for p in iter {
                *p = b.closest_image_dims(p, p0, dims);
            }
        }
        Ok(())
    }

    fn unwrap_simple(&mut self) -> Result<(), MeasureError> {
        self.unwrap_simple_dim(PBC_FULL)
    }
}

/// Trait for modification requiring random access positions and pbc
pub trait ModifyRandomAccess: PosIterProvider + PosIterMutProvider + BoxMutProvider + BoxProvider + RandomPosMutProvider + Sized {
    fn unwrap_connectivity(&mut self, cutoff: f32) -> Result<(), MeasureError> {
        self.unwrap_connectivity_dim(cutoff, PBC_FULL)
    }

    fn unwrap_connectivity_dim(&mut self, cutoff: f32, dims: PbcDims) -> Result<(), MeasureError> {
        let b = self.require_box()?.to_owned();
        let conn: SearchConnectivity =
            distance_search_single_pbc(cutoff, self.iter_pos(), 0..self.len(), &b, dims);

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
            let p0 = unsafe { self.get_pos_mut_unchecked(c) }.to_owned();
            // Iterate over connected points
            for ind in &conn[c] {
                // Unwrap this point if it is not used yet
                if !used[*ind] {
                    let p = unsafe { self.get_pos_mut_unchecked(*ind) };
                    *p = b.closest_image_dims(p, &p0, dims);
                    // Add it to the stack
                    todo.push(*ind);
                    used[*ind] = true;
                }
            }
            //println!(">> {:?}",todo);
        }

        if used.len() != conn.len() {
            Err(MeasureError::Disjoint)
        } else {
            Ok(())
        }
    }
}

/// Trait for modification requiring atoms
pub trait ModifyAtoms: AtomIterMutProvider + LenProvider {
    fn assign_resindex(&mut self) {
        let n = self.len();
        let mut resindex = 0usize;
        let mut at_iter = self.iter_atoms_mut();
        if n > 1 {
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
