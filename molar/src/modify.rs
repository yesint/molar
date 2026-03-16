use crate::prelude::*;
use itertools::Itertools;
use nalgebra::Const;
use nalgebra::Rotation3;
use nalgebra::Unit;

//==============================================================
// Unified Modify trait
//==============================================================

/// Unified trait for all modification operations on selections.
///
/// The base supertrait is [`PosMutProvider`]. Methods that require additional
/// capabilities carry `where Self: …` bounds.
pub trait Modify: PosMutProvider {
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

    // ---- BoxProvider methods ----

    fn unwrap_simple_dim(&mut self, dims: PbcDims) -> Result<(), MeasureError>
    where
        Self: BoxProvider,
    {
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

    fn unwrap_simple(&mut self) -> Result<(), MeasureError>
    where
        Self: BoxProvider,
    {
        self.unwrap_simple_dim(PBC_FULL)
    }

    // ---- BoxProvider + Sized + Selectable methods ----

    fn unwrap_connectivity(&mut self, cutoff: f32) -> Result<Vec<Sel>, MeasureError>
    where
        Self: BoxProvider + Sized + Selectable,
    {
        self.unwrap_connectivity_dim(cutoff, PBC_FULL)
    }

    fn unwrap_connectivity_dim(&mut self, cutoff: f32, dims: PbcDims) -> Result<Vec<Sel>, MeasureError>
    where
        Self: BoxProvider + Sized + Selectable,
    {
        let b = self.require_box()?.to_owned();
        let conn: SearchConnectivity =
            distance_search_single_pbc(cutoff, self.iter_pos(), 0..self.len(), &b, PBC_FULL);

        // used atoms
        let mut used = vec![false; self.len()];
        // Centers to unwrap
        let mut todo = Vec::<usize>::with_capacity(conn.len() / 2);
        // Place first center to the stack
        todo.push(0);
        used[0] = true;
        let mut sel_vec = vec![];
        let mut res_sels = vec![];
        loop {
            // Loop while stack is not empty
            while let Some(c) = todo.pop() {
                // Central point
                let p0 = unsafe { self.get_pos_mut_unchecked(c) }.to_owned();
                // Iterate over connected points
                if let Some(v) = conn.get(c) {
                    for ind in v {
                        // Unwrap this point if it is not used yet
                        if !used[*ind] {
                            let p = unsafe { self.get_pos_mut_unchecked(*ind) };
                            *p = b.closest_image_dims(p, &p0, dims);
                            // Add it to the stack
                            todo.push(*ind);
                            used[*ind] = true;
                            // Add to current sel
                            sel_vec.push(*ind);
                        }
                    }
                }
            }

            // Check if any unused points remained
            if let Some((i, _)) = used.iter().find_position(|el| **el == false) {
                // If any found, add it to used and go on
                todo.push(i);
                used[i] = true;
                // Create output selection
                if !sel_vec.is_empty() {
                    res_sels.push(self.select(&sel_vec).unwrap());
                }
                sel_vec.clear();
            } else {
                // Add remaining indices to the last selection
                if !sel_vec.is_empty() {
                    res_sels.push(self.select(&sel_vec).unwrap());
                }
                break;
            }
        }

        Ok(res_sels)
    }

    // ---- AtomMutProvider method ----

    fn assign_resindex(&mut self)
    where
        Self: AtomMutProvider,
    {
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
