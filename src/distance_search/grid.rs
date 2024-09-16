use num_traits::clamp_min;

use crate::core::{PbcDims, PeriodicBox, Pos, Vector3f};

//====================================================================
// Cell location in the grid
pub type CellLoc = nalgebra::Vector3<usize>;

/// Trait for grid items
/// Objects implementing this trait will work for distance search.
pub trait GridItem: Send + Sync { 
    fn get_pos(&self) -> &Pos;
    fn get_id(&self) -> usize;
}

impl GridItem for (usize,Pos) {
    fn get_pos(&self) -> &Pos {
        &self.1
    }
    
    fn get_id(&self) -> usize {
        self.0
    }
}

impl GridItem for (usize,Pos,f32) {    
    fn get_pos(&self) -> &Pos {
        &self.1
    }

    fn get_id(&self) -> usize {
        self.0
    }
}

//======================================================================
// The pair of cells for searching
#[derive(Debug, Clone)]
pub struct CellPair {
    pub c1: CellLoc,
    pub c2: CellLoc,
    pub wrapped: PbcDims,
}

impl CellPair {
    pub fn swaped(&self) -> CellPair {
        let mut ret = self.clone();
        std::mem::swap(&mut ret.c1, &mut ret.c2);
        ret
    }
}

// Periodicity stuff for searching
#[derive(Debug, Clone)]
pub struct GridPbc {
    pub pbox: PeriodicBox,
    pub dims: PbcDims,
}

//=============================================================

// Grid
#[derive(Debug, Clone)]
pub struct Grid<T> {
    pub data: ndarray::Array3<T>,
    pub pbc: Option<GridPbc>,
    pub n_items: usize, // Number of items in the grid
}

impl<T> Grid<T>
where
    T: Default,
{
    pub fn new(sz: [usize; 3]) -> Self {
        Self {
            data: ndarray::Array3::<T>::from_shape_simple_fn(sz, || T::default()),
            pbc: None,
            n_items: 0,
        }
    }

    pub fn from_cutoff_and_extents(cutoff: f32, extents: &Vector3f) -> Self {
        let mut sz = [0, 0, 0];
        // Cell size should be >= cutoff for all dimentions
        for d in 0..3 {
            sz[d] = clamp_min((extents[d] / cutoff).floor() as usize, 1);
        }
        Self::new(sz)
    }

    pub fn from_cutoff_and_min_max(cutoff: f32, min: &Vector3f, max: &Vector3f) -> Self {
        Self::from_cutoff_and_extents(cutoff, &(max - min))
    }

    pub fn from_cutoff_and_box(cutoff: f32, box_: &PeriodicBox) -> Self {
        Self::from_cutoff_and_extents(cutoff, &box_.get_extents())
    }

    pub fn dim(&self) -> [usize; 3] {
        let d = self.data.dim();
        [d.0, d.1, d.2]
    }
}

impl<I: GridItem> Grid<Vec<I>> 
{
    pub fn populate(
        &mut self,
        iter: impl Iterator<Item = I>,
        lower: &Vector3f,
        upper: &Vector3f,
    ) {
        let dim = self.dim();
        let dim_sz = upper - lower;
        'outer: for el in iter {
            let mut ind = [0usize, 0, 0];
            for d in 0..3 {
                let n = (dim[d] as f32 * (el.get_pos()[d] - lower[d]) / dim_sz[d]).floor() as isize;
                if n < 0 || n >= dim[d] as isize {
                    continue 'outer;
                } else {
                    ind[d] = n as usize;
                }
            }
            self.data[ind].push(el);
            self.n_items += 1;
        }
    }

    pub fn populate_periodic(
        &mut self,
        iter: impl Iterator<Item = I>,
        box_: &PeriodicBox,
        pbc_dims: &PbcDims,
    ) {
        let dim = &self.dim();

        'outer: for el in iter {
            // Relative coordinates
            let rel = box_.to_box_coords(&el.get_pos().coords);
            let mut ind = [0usize, 0, 0];
            for d in 0..3 {
                // If dimension in not periodic and
                // out of bounds - skip the point
                if !pbc_dims[d] && (rel[d] >= 1.0 || rel[d] < 0.0) {
                    continue 'outer;
                }
                ind[d] = (rel[d] * dim[d] as f32).floor().rem_euclid(dim[d] as f32) as usize;
            }
            //println!("{:?}",ind);
            self.data[ind].push(el);
            self.n_items += 1;
        }
        self.pbc = Some(GridPbc {
            pbox: box_.clone(),
            dims: pbc_dims.clone(),
        });
    }
}

// Implement indexing by nalgebra Vector3 for generic Grid
impl<T> std::ops::Index<&CellLoc> for Grid<T> {
    type Output = T;

    fn index(&self, index: &CellLoc) -> &Self::Output {
        &self.data[(index.x, index.y, index.z)]
    }
}

impl<T> std::ops::IndexMut<&CellLoc> for Grid<T> {
    fn index_mut(&mut self, index: &CellLoc) -> &mut Self::Output {
        &mut self.data[(index.x, index.y, index.z)]
    }
}

//============================================================================
#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    #[test]
    fn test_grid() {
        use crate::core::PBC_FULL;
        use crate::io::*;
        use std::iter::zip;

        let mut r = FileHandler::open("tests/protein.pdb").unwrap();
        let st = r.read_state().unwrap().unwrap();

        let mut gr = Grid::new([10, 10, 10]);
        gr.populate_periodic(
            zip(0..st.num_coords(), st.iter_pos().cloned()),
            &st.get_box().unwrap(),
            &PBC_FULL,
        );
    }
}
