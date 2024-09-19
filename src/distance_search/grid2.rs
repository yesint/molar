use num_traits::clamp_min;
use rayon::iter::{FromParallelIterator, IndexedParallelIterator, IntoParallelIterator};

use super::SearchOutputType;
use crate::prelude::*;

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
    ([0, 0, 0], [1, 1, 0]), //XYWO #1
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
    pub fn debug(&self) {
        for x in 0..self.dims[0] {
            for y in 0..self.dims[1] {
                for z in 0..self.dims[2] {
                    let n = self.cells[self.loc_to_ind([x,y,z])].len();
                    if n>0 {
                        println!("{},{},{} {}",x,y,z,n);
                    }
                }
            }
        }
    }

    fn new_with_dims(dims: [usize; 3]) -> Self {
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

    fn from_cutoff_and_extents(cutoff: f32, extents: &Vector3f) -> Self {
        let mut sz = [0, 0, 0];
        // Cell size should be >= cutoff for all dimentions
        for d in 0..3 {
            sz[d] = clamp_min((extents[d] / cutoff).floor() as usize, 1);
        }
        Self::new_with_dims(sz)
    }

    pub fn from_cutoff_and_min_max(cutoff: f32, min: &Vector3f, max: &Vector3f) -> Self {
        Self::from_cutoff_and_extents(cutoff, &(max - min))
    }

    pub fn from_cutoff_and_box(cutoff: f32, box_: &PeriodicBox) -> Self {
        Self::from_cutoff_and_extents(cutoff, &box_.get_extents())
    }

    pub fn populate(&mut self, data: &impl RandomPos, lower: &Vector3f, upper: &Vector3f) {
        // Data points are numbered sequentially from zero
        // So grid always stores the local index within the data
        let dim_sz = upper - lower;
        'outer: for id in 0..data.len() {
            let pos = unsafe {data.nth_pos_unchecked(id)};
            let mut loc = [0usize, 0, 0];
            for d in 0..3 {
                let n = (self.dims[d] as f32 * (pos[d] - lower[d]) / dim_sz[d]).floor() as isize;
                if n < 0 || n >= self.dims[d] as isize {
                    continue 'outer;
                } else {
                    loc[d] = n as usize;
                }
            }
            self.get_mut(loc).push(id);
        }
    }

    pub fn populate_pbc(
        &mut self,
        data: &impl PosProvider,
        box_: &PeriodicBox,
        pbc_dims: &PbcDims,
    ) {
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
                loc[d] = (rel[d] * self.dims[d] as f32)
                    .floor()
                    .rem_euclid(self.dims[d] as f32) as usize;
            }
            //println!("{:?}",ind);
            self.get_mut(loc).push(id);
        }
    }

    pub fn search_plan(
        grid1: &Grid3d,
        grid2: Option<&Grid3d>,
        periodic: bool,
    ) -> Vec<(usize, usize, bool)> {
        let mut plan = Vec::with_capacity(14 * grid1.dims[0] * grid1.dims[1] * grid1.dims[2]);
        // Cycle over whole grid
        for x in 0..grid1.dims[0] {
            for y in 0..grid1.dims[1] {
                for z in 0..grid1.dims[2] {
                    // go over possible pairs
                    'mask: for (v1, v2) in MASK {
                        let mut c = [
                            [x + v1[0], y + v1[1], z + v1[2]],
                            [x + v2[0], y + v2[1], z + v2[2]],
                        ];
                        // we only go to the right, so need to check the right edge
                        let mut wrapped = false;
                        for i in 0..=1 {
                            for d in 0..3 {
                                if c[i][d] == grid1.dims[d] {
                                    if periodic {
                                        c[i][d] = 0;
                                        wrapped = true;
                                    } else {
                                        continue 'mask;
                                    }
                                }
                            }
                        }
                        // If we are here we need to add te cell pair to then plan
                        let i1 = grid1.loc_to_ind(c[0]);
                        let i2 = grid1.loc_to_ind(c[1]);

                        // Check if there are points in both cells
                        // This is different in single and double grid cases
                        if let Some(grid2) = grid2 {
                            if (grid1.cells[i1].len() > 0 && grid2.cells[i2].len() > 0)
                                || (grid2.cells[i1].len() > 0 && grid1.cells[i2].len() > 0)
                            {
                                plan.push((i1, i2, wrapped));
                                //plan.push((i2, i1, wrapped));
                            }
                        } else if grid1.cells[i1].len() > 0 && grid1.cells[i2].len() > 0 {
                            plan.push((i1, i2, wrapped));
                        }
                    }
                }
            }
        }

        // Assess the work load for each pair
        // let mut n = 0;
        // if let Some(grid2) = grid2 {
        //     for p in plan.iter() {
        //         let mut nl = 0;
        //         nl += grid1.cells[p.0].len() * grid2.cells[p.1].len();
        //         nl += grid1.cells[p.1].len() * grid2.cells[p.0].len();
        //         n+=nl;
        //     println!("{}",nl);
        //     }
        // } else {
        //     for p in plan.iter() {
        //         n += grid1.cells[p.0].len() * grid1.cells[p.1].len();
        //     }
        // }
        // println!("pairs: {}, total N: {}",plan.len(),n);

        plan
    }
    
    pub fn get_dims(&self) -> [usize; 3] {
        self.dims
    }
}

fn search_cell_within(
    cutoff2: f32,
    grid1: &Grid3d,
    grid2: &Grid3d,
    pair: (usize, usize, bool),
    coords1: &impl RandomPos,
    coords2: &impl RandomPos,
    found: &mut Vec<usize>,
) {
    let n1 = grid1.cells[pair.0].len();
    let n2 = grid2.cells[pair.1].len();

    for i in 0..n1 {
        let ind1 = grid1.cells[pair.0][i];
        let pos1 = unsafe { *coords1.nth_pos_unchecked(ind1) };
        for j in 0..n2 {
            let ind2 = grid2.cells[pair.1][j];
            //println!("{} {}",ind1,ind2);
            let pos2 = unsafe { coords2.nth_pos_unchecked(ind2) };

            let d2 = (pos2 - pos1).norm_squared();
            if d2 <= cutoff2 {
                found.push(ind1);
                break;
            }
        }
    }
}

// fn search_cell_pair_pbc(
//     cutoff: f32,
//     grid1: &Grid3d,
//     grid2: &Grid3d,
//     pair: (usize, usize, bool),
//     coords1: &impl RandomPos,
//     coords2: &impl RandomPos,
//     box_: &PeriodicBox,
// ) -> Vec<usize> {
//     let n1 = grid1.cells[pair.0].len();
//     let n2 = grid2.cells[pair.1].len();

//     let mut found = Vec::new();

//     for i in 0..n1 {
//         let ind1 = grid1.cells[pair.0][i];
//         let pos1 = unsafe { coords1.nth_pos_unchecked(ind1) };
//         for j in 0..n2 {
//             let ind2 = grid2.cells[pair.0][j];
//             let pos2 = unsafe { coords2.nth_pos_unchecked(ind2) };

//             let d = box_.distance(pos1, pos2, &PBC_FULL);
//             if d <= cutoff && d>10.0*f32::EPSILON {
//                 found.push(ind1);
//             }
//         }
//     }

//     found
// }

pub(crate) fn distance_search_within<C>(
    cutoff: f32,
    data1: &(impl RandomPos + Send + Sync),
    data2: &(impl RandomPos + Send + Sync),
    lower: &Vector3f,
    upper: &Vector3f,
) -> C 
where    
    C: FromIterator<usize> + FromParallelIterator<usize>,
{
    let mut grid1 = Grid3d::from_cutoff_and_min_max(cutoff, lower,upper);
    let mut grid2 = Grid3d::new_with_dims(grid1.get_dims());

    grid1.populate(data1, lower, upper);
    grid2.populate(data2, lower, upper);

    //grid1.debug();

    let plan = Grid3d::search_plan(&grid1, Some(&grid2), false);
    
    // Cycle over search plan and perform search for each cell pair
    plan.into_par_iter().with_min_len(3).map(|pair| {
            let mut found = Vec::new();
            search_cell_within(cutoff*cutoff, &grid1, &grid2, pair, data1, data2,&mut found);
            search_cell_within(cutoff*cutoff, &grid1, &grid2, (pair.1,pair.0,pair.2), data1, data2,&mut found);
            found
        }
    ).flatten().collect()

    // Serial
    // let mut found = Vec::with_capacity(10000); 
    // for pair in plan {
    //     search_cell_pair::<T>(cutoff*cutoff, &grid1, &grid2, pair, data1, data2,&mut found);
    //     search_cell_pair::<T>(cutoff*cutoff, &grid1, &grid2, (pair.1,pair.0,pair.2), data1, data2,&mut found);
    // }
    // found.into_iter().collect()


    // Parallel: SubsetType
    // plan.into_par_iter().map(|pair|
    //     search_cell_pair::<usize>(cutoff, &grid1, &grid2, pair, &data1, &data2)
    // ).flatten().collect()
}


#[cfg(test)]
mod tests {
    use crate::core::Source;
    
    #[test]
    fn within_plan_test() -> anyhow::Result<()> {
        let mut src = Source::serial_from_file("tests/albumin.pdb")?;
        let t = std::time::Instant::now();
        src.select_str("within 4.0 of resid 10:300")?;
        println!("elapsed {}", t.elapsed().as_secs_f32());
        Ok(())
    }
}
