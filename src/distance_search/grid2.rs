use num_traits::clamp_min;

use crate::{
    core::{BoxProvider, MeasurePos, PbcDims, PeriodicBox, Pos, PosProvider, RandomAtom, RandomPos, RandomPosMut, Vector3f},
    io::IndexProvider,
};
use super::{grid::GridItem, SearchOutputType};

pub trait DistanceSearchePosProvider {
    fn nth_pos_unchecked(&self, n: usize) -> &Pos;
}

pub trait DistanceSearchePosVdwProvider {
    fn nth_pos_vdw_unchecked(&self, n: usize) -> (&Pos,f32);
}

struct Grid3d {
    cells: Vec<Vec<usize>>,
    dims: [usize; 3],
}

static MASK: [([usize; 3], [usize; 3]); 14] = [
    // Center
    ([0, 0, 0], [0, 0, 0]),
    // Edges
    ([0, 0, 0], [1, 0, 0]), //X
    ([0, 0, 0], [0, 1, 0]), //Y
    ([0, 0, 0], [0, 0, 1]), //Z
    // Face angles
    ([0, 0, 0], [1, 1, 0]), //XY
    ([0, 0, 0], [1, 0, 1]), //XZ
    ([0, 0, 0], [0, 1, 1]), //YZ
    // Far angls
    ([0, 0, 0], [1, 1, 1]), //XYZ
    // Face-diagonals
    ([1, 0, 0], [0, 1, 0]), // XY
    ([1, 0, 0], [0, 0, 1]), // XZ
    ([0, 1, 0], [0, 0, 1]), // YZ
    // Cross-diagonals
    ([1, 1, 0], [0, 0, 1]), // XY-Z
    ([1, 0, 1], [0, 1, 0]), // XZ-Y
    ([0, 1, 1], [1, 0, 0]), // YZ-X
];

impl Grid3d {
    fn new_internal(dims: [usize; 3]) -> Self {
        Self {
            cells: vec![vec![]; dims[0] * dims[1] * dims[2]],
            dims,
        }
    }

    #[inline(always)]
    fn loc_to_ind(&self, loc: [usize; 3]) -> usize {
        loc[0] + loc[1] * self.dims[0] + loc[2] * self.dims[0] * self.dims[1]
    }

    pub fn get(&self, loc: [usize; 3]) -> &Vec<usize> {
        &self.cells[self.loc_to_ind(loc)]
    }

    pub fn get_mut(&mut self, loc: [usize; 3]) -> &mut Vec<usize> {
        let i = self.loc_to_ind(loc);
        &mut self.cells[i]
    }

    fn dims_from_cutoff_and_extents(cutoff: f32, extents: &Vector3f) -> [usize; 3] {
        let mut sz = [0, 0, 0];
        // Cell size should be >= cutoff for all dimentions
        for d in 0..3 {
            sz[d] = clamp_min((extents[d] / cutoff).floor() as usize, 1);
        }
        sz
    }

    pub fn new(cutoff: f32, data: &impl MeasurePos) -> Self {
        // Find extents
        let (lower, upper) = data.min_max();
        let dim_sz = upper - lower;
        // Compute grid dimensions
        let dims = Self::dims_from_cutoff_and_extents(cutoff, &dim_sz);
        let mut slf = Self::new_internal(dims);

        // Data points are numbered sequentially from zero
        // So grid always stores the local index within the data
        'outer: for (id, pos) in data.iter_pos().enumerate() {
            let mut loc = [0usize, 0, 0];
            for d in 0..3 {
                let n = (slf.dims[d] as f32 * (pos[d] - lower[d]) / dim_sz[d]).floor() as isize;
                if n < 0 || n >= slf.dims[d] as isize {
                    continue 'outer;
                } else {
                    loc[d] = n as usize;
                }
            }
            slf.get_mut(loc).push(id);
        }
        slf
    }

    pub fn new_pbc(
        cutoff: f32,
        data: &impl PosProvider,
        box_: &PeriodicBox,
        pbc_dims: &PbcDims,
    ) -> Self {
        // Compute grid dimensionsdim_sz
        let dims = Self::dims_from_cutoff_and_extents(cutoff, &box_.get_extents());
        let mut slf = Self::new_internal(dims);

        'outer: for (id, pos) in data.iter_pos().enumerate() {
            // Relative coordinates
            let rel = box_.to_box_coords(&pos.coords);
            let mut loc = [0usize, 0, 0];
            for d in 0..3 {
                // If dimension in not periodic and
                // out of bounds - skip the point
                if !pbc_dims[d] && (rel[d] >= 1.0 || rel[d] < 0.0) {
                    continue 'outer;
                }
                loc[d] = (rel[d] * slf.dims[d] as f32)
                    .floor()
                    .rem_euclid(slf.dims[d] as f32) as usize;
            }
            //println!("{:?}",ind);
            slf.get_mut(loc).push(id);
        }
        slf
    }

    pub fn search_plan_single(grid: &Grid3d, periodic: bool) -> Vec<(usize, usize, bool)> {
        let mut plan = Vec::with_capacity(14 * grid.dims[0] * grid.dims[1] * grid.dims[2]);
        // Cycle over whole grid
        for x in 0..grid.dims[0] {
            for y in 0..grid.dims[1] {
                for z in 0..grid.dims[2] {
                    // go over possible pairs
                    for (v1, v2) in MASK {
                        let mut c1 = [x + v1[0], y + v1[1], z + v1[2]];
                        let mut c2 = [x + v2[0], y + v2[1], z + v2[2]];
                        // we only go to the right, so need to check the right edge
                        let mut wrapped = false;
                        for d in 0..3 {
                            if c1[d] == grid.dims[d] {
                                if periodic {
                                    c1[d] = 0;
                                    wrapped = true;
                                } else {
                                    continue;
                                }
                            }
                            if c2[d] == grid.dims[d] {
                                if periodic {
                                    c2[d] = 0;
                                    wrapped = true;
                                } else {
                                    continue;
                                }
                            }
                        }
                        // If we are here we need to add te cell pair to then plan
                        // Check if there are points in both cells
                        let i1 = grid.loc_to_ind(c1);
                        let i2 = grid.loc_to_ind(c2);
                        if grid.cells[i1].len()>0 && grid.cells[i2].len()>0 {
                            plan.push((i1, i2, wrapped));
                        }
                    }
                }
            }
        }
        plan
    }

    pub fn search_plan_double(
        grid1: &Grid3d,
        grid2: &Grid3d,
        periodic: bool,
    ) -> Vec<(usize, usize, bool)> {
        let mut plan = Vec::with_capacity(14 * grid1.dims[0] * grid1.dims[1] * grid1.dims[2]);
        // Cycle over whole grid
        for x in 0..grid1.dims[0] {
            for y in 0..grid1.dims[1] {
                for z in 0..grid1.dims[2] {
                    // go over possible pairs
                    for (v1, v2) in MASK {
                        let mut c1 = [x + v1[0], y + v1[1], z + v1[2]];
                        let mut c2 = [x + v2[0], y + v2[1], z + v2[2]];
                        // we only go to the right, so need to check the right edge
                        let mut wrapped = false;
                        for d in 0..3 {
                            if c1[d] == grid1.dims[d] {
                                if periodic {
                                    c1[d] = 0;
                                    wrapped = true;
                                } else {
                                    continue;
                                }
                            }
                            if c2[d] == grid1.dims[d] {
                                if periodic {
                                    c2[d] = 0;
                                    wrapped = true;
                                } else {
                                    continue;
                                }
                            }
                        }
                        // If we are here we need to add te cell pair to then plan
                        // Check if there are points in both cells
                        // Check both combinations
                        let i1 = grid1.loc_to_ind(c1);
                        let i2 = grid1.loc_to_ind(c2);
                        if (grid1.cells[i1].len()>0 && grid2.cells[i2].len()>0)
                        || (grid2.cells[i1].len()>0 && grid1.cells[i2].len()>0)
                        {
                            plan.push((i1, i2, wrapped));
                        }
                    }
                }
            }
        }
        plan
    }
}


fn search_cell_pair<T: SearchOutputType>(
    grid: &Grid3d,
    pair: (usize,usize,bool),    
    coords: impl RandomPos,    
    accept_func: fn(&Pos, &Pos, bool) -> Option<f32>,
) -> Vec<T> {
    let n1 = grid.cells[pair.0].len();
    let n2 = grid.cells[pair.1].len();

    let mut found = Vec::<T>::new();

    if pair.0 == pair.1 {
        // Same cell
        for i in 0..n1 - 1 {
            let ind1 = grid.cells[pair.0][i];
            let p1 = unsafe{ coords.nth_pos_unchecked(ind1) };
            for j in i + 1..n1 {
                let ind2 = grid.cells[pair.1][j];
                let p2 = unsafe{ coords.nth_pos_unchecked(ind2) };
                if let Some(d) = accept_func(p1, p2, pair.2) {
                    found.push(T::from_search_results(i,j,d));
                }
            }
        }
    } else {
        // Different cells
        for i in 0..n1 {
            let ind1 = grid.cells[pair.0][i];
            let p1 = unsafe{ coords.nth_pos_unchecked(ind1) };
            for j in 0..n2 {
                let ind2 = grid.cells[pair.1][j];
                let p2 = unsafe{ coords.nth_pos_unchecked(ind2) };
                if let Some(d) = accept_func(p1, p2, pair.2) {
                    found.push(T::from_search_results(i,j,d));
                }
            }
        }
    }

    found
}

// pub fn distance_search_single<I: GridItem>(cutoff: f32, data: impl std::ops::Index<usize, Output = I>) {
//     let grid = Grid3d::new(cutoff, data);
//     let plan = Grid3d::search_plan_single(&grid, false);
//     // Cycle over search plan and perform seacrh for each cell pair
//     for pair in plan {
//         search_cell_pair(&grid, pair, data, accept_func)
//     }
// }
