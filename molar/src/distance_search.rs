use crate::prelude::*;
use num_traits::clamp_min;
use rayon::iter::{FromParallelIterator, IndexedParallelIterator, IntoParallelIterator};

pub trait SearchOutputType {
    fn from_ijd(i: usize, j: usize, d: f32) -> Self;
}

impl SearchOutputType for usize {
    fn from_ijd(i: usize, _j: usize, _d: f32) -> Self {
        i
    }
}

impl SearchOutputType for (usize, usize) {
    fn from_ijd(i: usize, j: usize, _d: f32) -> Self {
        (i, j)
    }
}

impl SearchOutputType for (usize, usize, f32) {
    fn from_ijd(i: usize, j: usize, d: f32) -> Self {
        (i, j, d)
    }
}

//---------------------------------------------------------------
#[derive(Debug, Default)]
pub struct SearchConnectivity(rustc_hash::FxHashMap<usize, Vec<usize>>);

impl SearchConnectivity {
    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn iter(&self) -> SearchConnectivityIter {
        SearchConnectivityIter(self.0.iter())
    }

    fn from_iter(iter: impl IntoIterator<Item = (usize, usize)>) -> Self {
        let mut res = Self::default();
        for (i, j) in iter {
            if res.0.contains_key(&i) {
                res.0.get_mut(&i).unwrap().push(j);
            } else {
                res.0.insert(i, vec![j]);
            }

            if res.0.contains_key(&j) {
                res.0.get_mut(&j).unwrap().push(i);
            } else {
                res.0.insert(j, vec![i]);
            }
        }
        res
    }
}

impl FromIterator<(usize, usize)> for SearchConnectivity {
    fn from_iter<T: IntoIterator<Item = (usize, usize)>>(iter: T) -> Self {
        Self::from_iter(iter)
    }
}

impl FromParallelIterator<(usize, usize)> for SearchConnectivity {
    fn from_par_iter<I>(par_iter: I) -> Self
    where
        I: IntoParallelIterator<Item = (usize, usize)>,
    {
        let v = par_iter.into_par_iter().collect::<Vec<_>>();
        Self::from_iter(v.iter().cloned())
    }
}

pub struct SearchConnectivityIter<'a>(std::collections::hash_map::Iter<'a, usize, Vec<usize>>);

impl<'a> Iterator for SearchConnectivityIter<'a> {
    type Item = (&'a usize, &'a Vec<usize>);
    fn next(&mut self) -> Option<Self::Item> {
        self.0.next()
    }
}

impl IntoIterator for SearchConnectivity {
    type Item = (usize, Vec<usize>);
    type IntoIter = std::collections::hash_map::IntoIter<usize, Vec<usize>>;
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl std::ops::Index<usize> for SearchConnectivity {
    type Output = Vec<usize>;
    fn index(&self, i: usize) -> &Self::Output {
        &self.0[&i]
    }
}

//--------------------------------------------------------------------------------

struct Grid<'a> {
    cells: Vec<Vec<(usize, &'a Pos)>>,
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

impl<'a> Grid<'a> {
    // pub fn debug(&self) {
    //     for x in 0..self.dims[0] {
    //         for y in 0..self.dims[1] {
    //             for z in 0..self.dims[2] {
    //                 let n = self.cells[self.loc_to_ind([x,y,z])].len();
    //                 if n>0 {
    //                     println!("{},{},{} {}",x,y,z,n);
    //                 }
    //             }
    //         }
    //     }
    // }

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

    // pub fn get_loc_mut(&mut self, loc: [usize; 3]) -> &mut Vec<(usize,&'a Pos)> {
    //     let i = self.loc_to_ind(loc);
    //     &mut self.cells[i]
    // }

    fn push_loc(&mut self, loc: [usize; 3], data: (usize, &'a Pos)) {
        let i = self.loc_to_ind(loc);
        self.cells[i].push(data);
    }

    fn push_ind(&mut self, ind: usize, data: (usize, &'a Pos)) {
        self.cells[ind].push(data);
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
        Self::from_cutoff_and_extents(cutoff, &box_.get_lab_extents())
    }

    pub fn populate(
        &mut self,
        data: impl PosIterator<'a>,
        ids: impl Iterator<Item = usize>,
        lower: &Vector3f,
        upper: &Vector3f,
    ) {
        // Data points are numbered sequentially from zero
        // So grid always stores the local index within the data
        let dim_sz = upper - lower;
        'outer: for (id, pos) in ids.zip(data) {
            let mut loc = [0usize, 0, 0];
            for d in 0..3 {
                let n = (self.dims[d] as f32 * (pos[d] - lower[d]) / dim_sz[d]).floor() as isize;
                if n < 0 || n >= self.dims[d] as isize {
                    continue 'outer;
                } else {
                    loc[d] = n as usize;
                }
            }
            self.push_loc(loc, (id, pos));
        }
    }

    pub fn populate_pbc(
        &mut self,
        data: impl PosIterator<'a>,
        ids: impl Iterator<Item = usize>,
        box_: &PeriodicBox,
        pbc_dims: PbcDims,
        wrapped_pos: &'a mut Vec<Pos>,
    ) {
        let mut wrapped_ind = vec![];
        wrapped_pos.clear();

        'outer: for (id, pos) in ids.zip(data) {
            // Relative coordinates
            let mut rel = box_.to_box_coords(&pos.coords);
            let mut loc = [0usize, 0, 0];

            // Check if point is correctly wrapped
            let mut correct = true;
            for d in 0..3 {
                if rel[d] < 0.0 || rel[d] >= 1.0 {
                    if !pbc_dims.get_dim(d) {
                        continue 'outer;
                    } else {
                        correct = false;
                        break;
                    }
                }
            }

            if correct {
                // Wrapped correctly
                for d in 0..3 {
                    loc[d] = ((rel[d] * self.dims[d] as f32).floor() as usize)
                        .clamp(0, self.dims[d] - 1);
                    // Accounts for float point errors when loc[d] could be 1.00001
                }
                self.push_loc(loc, (id, pos));
            } else {
                // Need to wrap the point
                for d in 0..3 {
                    if pbc_dims.get_dim(d) {
                        // For each periodic dim
                        rel[d] = rel[d].fract();
                        if rel[d] < 0.0 {
                            rel[d] = 1.0 - rel[d]
                        }
                    }
                    loc[d] = ((rel[d] * self.dims[d] as f32).floor() as usize)
                        .clamp(0, self.dims[d] - 1);
                    // Accounts for float point errors when loc[d] could be 1.00001
                }
                let wp = Pos::from(box_.to_lab_coords(&rel));
                wrapped_pos.push(wp);
                wrapped_ind.push((self.loc_to_ind(loc), id));
            }
        }

        // Add wrapped points to the grid if any
        for i in 0..wrapped_ind.len() {
            self.push_ind(wrapped_ind[i].0, (wrapped_ind[i].1, &wrapped_pos[i]));
        }
    }

    pub fn get_dims(&self) -> [usize; 3] {
        self.dims
    }
}

fn search_plan(
    grid1: &Grid,
    grid2: Option<&Grid>,
    pbc_dims: PbcDims,
) -> Vec<(usize, usize, PbcDims)> {
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
                    let mut wrapped = PBC_NONE;
                    for i in 0..=1 {
                        for d in 0..3 {
                            if c[i][d] == grid1.dims[d] {
                                if pbc_dims.get_dim(d) {
                                    c[i][d] = 0;
                                    wrapped.set_dim(d, true);
                                } else {
                                    // Drop point for non-periodic dimension
                                    continue 'mask;
                                }
                            }
                        }
                    }
                    // If we are here we need to add the cell pair to then plan
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
    plan
}

fn search_cell_pair_within(
    cutoff2: f32,
    grid1: &Grid,
    grid2: &Grid,
    pair: (usize, usize, PbcDims),
    found: &mut Vec<usize>,
) {
    let n1 = grid1.cells[pair.0].len();
    let n2 = grid2.cells[pair.1].len();

    for i in 0..n1 {
        let (ind1, pos1) = grid1.cells[pair.0][i];
        for j in 0..n2 {
            let (_, pos2) = grid2.cells[pair.1][j];

            let d2 = (pos2 - pos1).norm_squared();
            if d2 <= cutoff2 {
                found.push(ind1);
                break;
            }
        }
    }
}

fn search_cell_pair_within_pbc(
    cutoff2: f32,
    grid1: &Grid,
    grid2: &Grid,
    pair: (usize, usize, PbcDims),
    pbox: &PeriodicBox,
    found: &mut Vec<usize>,
) {
    let n1 = grid1.cells[pair.0].len();
    let n2 = grid2.cells[pair.1].len();

    for i in 0..n1 {
        let (ind1, pos1) = grid1.cells[pair.0][i];
        for j in 0..n2 {
            let (_, pos2) = grid2.cells[pair.1][j];

            let d2 = if pair.2.any() {
                pbox.distance_squared(&pos1, &pos2, pair.2)
            } else {
                (pos2 - pos1).norm_squared()
            };
            if d2 <= cutoff2 {
                found.push(ind1);
                break;
            }
        }
    }
}

fn search_cell_pair_double<T: SearchOutputType>(
    cutoff2: f32,
    grid1: &Grid,
    grid2: &Grid,
    pair: (usize, usize, PbcDims),
    found: &mut Vec<T>,
) {
    let n1 = grid1.cells[pair.0].len();
    let n2 = grid2.cells[pair.1].len();

    for i in 0..n1 {
        let (ind1, pos1) = grid1.cells[pair.0][i];
        for j in 0..n2 {
            let (ind2, pos2) = grid2.cells[pair.1][j];

            let d2 = (pos2 - pos1).norm_squared();
            if d2 <= cutoff2 {
                found.push(T::from_ijd(ind1, ind2, d2.sqrt()));
            }
        }
    }
}

fn search_cell_pair_double_pbc<T: SearchOutputType>(
    cutoff2: f32,
    grid1: &Grid,
    grid2: &Grid,
    pair: (usize, usize, PbcDims),
    pbox: &PeriodicBox,
    found: &mut Vec<T>,
) {
    let n1 = grid1.cells[pair.0].len();
    let n2 = grid2.cells[pair.1].len();

    for i in 0..n1 {
        let (ind1, pos1) = grid1.cells[pair.0][i];
        for j in 0..n2 {
            let (ind2, pos2) = grid2.cells[pair.1][j];

            let d2 = if pair.2.any() {
                pbox.distance_squared(&pos1, &pos2, pair.2)
            } else {
                (pos2 - pos1).norm_squared()
            };
            if d2 <= cutoff2 {
                found.push(T::from_ijd(ind1, ind2, d2.sqrt()));
            }
        }
    }
}

fn search_cell_pair_double_vdw<T: SearchOutputType>(
    grid1: &Grid,
    grid2: &Grid,
    pair: (usize, usize, PbcDims),
    vdw1: &Vec<f32>,
    vdw2: &Vec<f32>,
    found: &mut Vec<T>,
) {
    let n1 = grid1.cells[pair.0].len();
    let n2 = grid2.cells[pair.1].len();

    for i in 0..n1 {
        let (ind1, pos1) = grid1.cells[pair.0][i];
        for j in 0..n2 {
            let (ind2, pos2) = grid2.cells[pair.1][j];

            let d2 = (pos2 - pos1).norm_squared();
            let cutoff = vdw1[ind1] + vdw2[ind2] + f32::EPSILON;
            if d2 <= cutoff * cutoff {
                found.push(T::from_ijd(ind1, ind2, d2.sqrt()));
            }
        }
    }
}

fn search_cell_pair_double_vdw_pbc<T: SearchOutputType>(
    grid1: &Grid,
    grid2: &Grid,
    pair: (usize, usize, PbcDims),
    vdw1: &Vec<f32>,
    vdw2: &Vec<f32>,
    pbox: &PeriodicBox,
    found: &mut Vec<T>,
) {
    let n1 = grid1.cells[pair.0].len();
    let n2 = grid2.cells[pair.1].len();

    for i in 0..n1 {
        let (ind1, pos1) = grid1.cells[pair.0][i];
        for j in 0..n2 {
            let (ind2, pos2) = grid2.cells[pair.1][j];

            let d2 = if pair.2.any() {
                pbox.distance_squared(&pos1, &pos2, pair.2)
            } else {
                (pos2 - pos1).norm_squared()
            };

            let cutoff = vdw1[ind1] + vdw2[ind2] + f32::EPSILON;

            if d2 <= cutoff * cutoff {
                found.push(T::from_ijd(ind1, ind2, d2.sqrt()));
            }
        }
    }
}

fn search_cell_pair_single<T: SearchOutputType>(
    cutoff2: f32,
    grid: &Grid,
    pair: (usize, usize, PbcDims),
) -> Vec<T> {
    let mut found = Vec::<T>::new();

    if pair.0 == pair.1 {
        let n = grid.cells[pair.0].len();
        for i in 0..n - 1 {
            let (ind1, pos1) = grid.cells[pair.0][i];
            for j in i + 1..n {
                let (ind2, pos2) = grid.cells[pair.0][j];

                let d2 = (pos2 - pos1).norm_squared();
                if d2 <= cutoff2 {
                    found.push(T::from_ijd(ind1, ind2, d2.sqrt()));
                }
            }
        }
    } else {
        let n1 = grid.cells[pair.0].len();
        let n2 = grid.cells[pair.1].len();
        for i in 0..n1 {
            let (ind1, pos1) = grid.cells[pair.0][i];
            for j in 0..n2 {
                let (ind2, pos2) = grid.cells[pair.1][j];

                let d2 = (pos2 - pos1).norm_squared();
                if d2 <= cutoff2 {
                    found.push(T::from_ijd(ind1, ind2, d2.sqrt()));
                }
            }
        }
    }
    found
}

fn search_cell_pair_single_pbc<T: SearchOutputType>(
    cutoff2: f32,
    grid: &Grid,
    pair: (usize, usize, PbcDims),
    pbox: &PeriodicBox,
) -> Vec<T> {
    let mut found = Vec::<T>::new();

    if pair.0 == pair.1 {
        let n = grid.cells[pair.0].len();
        for i in 0..n - 1 {
            let (ind1, pos1) = grid.cells[pair.0][i];
            for j in i + 1..n {
                let (ind2, pos2) = grid.cells[pair.0][j];

                let d2 = if pair.2.any() {
                    pbox.distance_squared(&pos1, &pos2, pair.2)
                } else {
                    (pos2 - pos1).norm_squared()
                };

                if d2 <= cutoff2 {
                    found.push(T::from_ijd(ind1, ind2, d2.sqrt()));
                }
            }
        }
    } else {
        let n1 = grid.cells[pair.0].len();
        let n2 = grid.cells[pair.1].len();
        for i in 0..n1 {
            let (ind1, pos1) = grid.cells[pair.0][i];
            for j in 0..n2 {
                let (ind2, pos2) = grid.cells[pair.1][j];

                let d2 = if pair.2.any() {
                    pbox.distance_squared(&pos1, &pos2, pair.2)
                } else {
                    (pos2 - pos1).norm_squared()
                };

                if d2 <= cutoff2 {
                    found.push(T::from_ijd(ind1, ind2, d2.sqrt()));
                }
            }
        }
    }
    found
}

pub(crate) fn distance_search_within<'a, C>(
    cutoff: f32,
    data1: &impl PosProvider,
    data2: &impl PosProvider,
    ids1: impl Iterator<Item = usize>,
    ids2: impl Iterator<Item = usize>,
    lower: &Vector3f,
    upper: &Vector3f,
) -> C
where
    C: FromIterator<usize> + FromParallelIterator<usize>,
{
    let mut grid1 = Grid::from_cutoff_and_min_max(cutoff, lower, upper);
    let mut grid2 = Grid::new_with_dims(grid1.get_dims());

    grid1.populate(data1.iter_pos(), ids1, lower, upper);
    grid2.populate(data2.iter_pos(), ids2, lower, upper);

    //grid1.debug();

    let plan = search_plan(&grid1, Some(&grid2), PBC_NONE);

    // Cycle over search plan and perform search for each cell pair
    plan.into_par_iter()
        .with_min_len(3)
        .map(|pair| {
            let mut found = Vec::new();
            search_cell_pair_within(
                cutoff * cutoff,
                &grid1,
                &grid2,
                pair,
                &mut found
            );
            search_cell_pair_within(
                cutoff * cutoff,
                &grid1,
                &grid2,
                (pair.1, pair.0, pair.2),
                &mut found,
            );
            found
        })
        .flatten()
        .collect()
}

pub(crate) fn distance_search_within_pbc<C>(
    cutoff: f32,
    data1: &impl PosProvider,
    data2: &impl PosProvider,
    ids1: impl Iterator<Item = usize>,
    ids2: impl Iterator<Item = usize>,
    pbox: &PeriodicBox,
    pbc_dims: PbcDims,
) -> C
where
    C: FromIterator<usize> + FromParallelIterator<usize>,
{
    let mut grid1 = Grid::from_cutoff_and_box(cutoff, pbox);
    let mut grid2 = Grid::new_with_dims(grid1.get_dims());

    let mut buf1 = vec![];
    let mut buf2 = vec![];

    grid1.populate_pbc(data1.iter_pos(), ids1, pbox, pbc_dims, &mut buf1);
    grid2.populate_pbc(data2.iter_pos(), ids2, pbox, pbc_dims, &mut buf2);

    let plan = search_plan(&grid1, Some(&grid2), pbc_dims);

    // Cycle over search plan and perform search for each cell pair
    plan.into_par_iter()
        .with_min_len(3)
        .map(|pair| {
            let mut found = Vec::new();
            search_cell_pair_within_pbc(
                cutoff * cutoff,
                &grid1,
                &grid2,
                pair,
                pbox,
                &mut found
            );
            search_cell_pair_within_pbc(
                cutoff * cutoff,
                &grid1,
                &grid2,
                (pair.1, pair.0, pair.2),
                pbox,
                &mut found,
            );
            found
        })
        .flatten()
        .collect()
}

//-------------------------------------------------------------------------

fn compute_min_max(data: &impl PosProvider) -> (Vector3f, Vector3f) {
    let mut lower = Vector3f::zeros();
    let mut upper = Vector3f::zeros();
    for p in data.iter_pos() {
        for d in 0..3 {
            if p[d] < lower[d] {
                lower[d] = p[d];
            }
            if p[d] > upper[d] {
                upper[d] = p[d];
            }
        }
    }
    (lower, upper)
}

fn compute_bounding_box_double(
    cutoff: f32,
    data1: &impl PosProvider,
    data2: &impl PosProvider,
) -> (Vector3f, Vector3f) {
    let (l1, u1) = compute_min_max(data1);
    let (l2, u2) = compute_min_max(data2);

    let mut l = Vector3f::zeros();
    let mut u = Vector3f::zeros();
    for d in 0..3 {
        l[d] = l1[d].min(l2[d]);
        u[d] = u1[d].max(u2[d]);
    }

    l.add_scalar_mut(-cutoff - f32::EPSILON);
    u.add_scalar_mut(cutoff + f32::EPSILON);
    (l, u)
}

fn compute_bounding_box_single(cutoff: f32, data: &impl PosProvider) -> (Vector3f, Vector3f) {
    let (mut l, mut u) = compute_min_max(data);
    l.add_scalar_mut(-cutoff - f32::EPSILON);
    u.add_scalar_mut(cutoff + f32::EPSILON);
    (l, u)
}

pub fn distance_search_double<T, C>(
    cutoff: f32,
    data1: &impl PosProvider,
    data2: &impl PosProvider,
    ids1: impl Iterator<Item = usize>,
    ids2: impl Iterator<Item = usize>,
) -> C
where
    T: SearchOutputType + Send + Sync,
    C: FromIterator<T> + FromParallelIterator<T>,
{
    // Compute the extents
    let (lower, upper) = compute_bounding_box_double(cutoff, data1, data2);

    let mut grid1 = Grid::from_cutoff_and_min_max(cutoff, &lower, &upper);
    let mut grid2 = Grid::new_with_dims(grid1.get_dims());

    grid1.populate(data1.iter_pos(), ids1, &lower, &upper);
    grid2.populate(data2.iter_pos(), ids2, &lower, &upper);

    let plan = search_plan(&grid1, Some(&grid2), PBC_NONE);

    // Cycle over search plan and perform search for each cell pair
    plan.into_par_iter()
        .with_min_len(3)
        .map(|pair| {
            let mut found = Vec::new();
            search_cell_pair_double(
                cutoff * cutoff,
                &grid1,
                &grid2,
                pair,
                &mut found
            );
            search_cell_pair_double(
                cutoff * cutoff,
                &grid1,
                &grid2,
                (pair.1, pair.0, pair.2),
                &mut found,
            );
            found
        })
        .flatten()
        .collect()
}

pub fn distance_search_double_pbc<T, C>(
    cutoff: f32,
    data1: &(impl PosProvider + BoxProvider),
    data2: &impl PosProvider,
    ids1: impl Iterator<Item = usize>,
    ids2: impl Iterator<Item = usize>,
    pbox: &PeriodicBox,
    pbc_dims: PbcDims,
) -> C
where
    T: SearchOutputType + Send + Sync,
    C: FromIterator<T> + FromParallelIterator<T>,
{
    let mut grid1 = Grid::from_cutoff_and_box(cutoff, pbox);
    let mut grid2 = Grid::new_with_dims(grid1.get_dims());

    let mut buf1 = vec![];
    let mut buf2 = vec![];

    grid1.populate_pbc(data1.iter_pos(), ids1, pbox, pbc_dims, &mut buf1);
    grid2.populate_pbc(data2.iter_pos(), ids2, pbox, pbc_dims, &mut buf2);

    let plan = search_plan(&grid1, Some(&grid2), pbc_dims);

    // Cycle over search plan and perform search for each cell pair
    plan.into_par_iter()
        .with_min_len(3)
        .map(|pair| {
            let mut found = Vec::new();
            search_cell_pair_double_pbc(
                cutoff * cutoff,
                &grid1,
                &grid2,
                pair,
                pbox,
                &mut found
            );
            search_cell_pair_double_pbc(
                cutoff * cutoff,
                &grid1,
                &grid2,
                (pair.1, pair.0, pair.2),
                pbox,
                &mut found,
            );
            found
        })
        .flatten()
        .collect()
}

pub fn distance_search_double_vdw<T, C>(
    data1: &impl PosProvider,
    data2: &impl PosProvider,
    ids1: impl Iterator<Item = usize>,
    ids2: impl Iterator<Item = usize>,
    vdw1: &Vec<f32>,
    vdw2: &Vec<f32>,
) -> C
where
    T: SearchOutputType + Send + Sync,
    C: FromIterator<T> + FromParallelIterator<T>,
{
    // We need to find the largest VdW distance pair to get the grid extents
    let cutoff = vdw1.iter().cloned().reduce(f32::max).unwrap()
        + vdw2.iter().cloned().reduce(&f32::max).unwrap()
        + f32::EPSILON;

    // Compute the extents
    let (lower, upper) = compute_bounding_box_double(cutoff, data1, data2);

    let mut grid1 = Grid::from_cutoff_and_min_max(cutoff, &lower, &upper);
    let mut grid2 = Grid::new_with_dims(grid1.get_dims());

    grid1.populate(data1.iter_pos(), ids1, &lower, &upper);
    grid2.populate(data2.iter_pos(), ids2, &lower, &upper);

    let plan = search_plan(&grid1, Some(&grid2), PBC_NONE);

    // Cycle over search plan and perform search for each cell pair
    plan.into_par_iter()
        .with_min_len(3)
        .map(|pair| {
            let mut found = Vec::new();
            search_cell_pair_double_vdw(
                &grid1,
                &grid2,
                pair,
                vdw1,
                vdw2,
                &mut found
            );
            search_cell_pair_double_vdw(
                &grid1,
                &grid2,
                (pair.1, pair.0, pair.2),
                vdw1,
                vdw2,
                &mut found,
            );
            found
        })
        .flatten()
        .collect()
}

pub fn distance_search_double_vdw_pbc<T, C>(
    data1: &(impl PosProvider + BoxProvider),
    data2: &impl PosProvider,
    ids1: impl Iterator<Item = usize>,
    ids2: impl Iterator<Item = usize>,
    vdw1: &Vec<f32>,
    vdw2: &Vec<f32>,
    pbox: &PeriodicBox,
    pbc_dims: PbcDims,
) -> C
where
    T: SearchOutputType + Send + Sync,
    C: FromIterator<T> + FromParallelIterator<T>,
{
    // We need to find the largest VdW distance pair to get the grid extents
    let cutoff = vdw1.iter().cloned().reduce(f32::max).unwrap()
        + vdw2.iter().cloned().reduce(&f32::max).unwrap()
        + f32::EPSILON;

    let mut grid1 = Grid::from_cutoff_and_box(cutoff, pbox);
    let mut grid2 = Grid::new_with_dims(grid1.get_dims());

    let mut buf1 = vec![];
    let mut buf2 = vec![];

    grid1.populate_pbc(data1.iter_pos(), ids1, pbox, pbc_dims, &mut buf1);
    grid2.populate_pbc(data2.iter_pos(), ids2, pbox, pbc_dims, &mut buf2);

    let plan = search_plan(&grid1, Some(&grid2), pbc_dims);

    // Cycle over search plan and perform search for each cell pair
    plan.into_par_iter()
        .with_min_len(3)
        .map(|pair| {
            let mut found = Vec::new();
            search_cell_pair_double_vdw_pbc(
                &grid1,
                &grid2,
                pair,
                vdw1,
                vdw2,
                pbox,
                &mut found
            );
            search_cell_pair_double_vdw_pbc(
                &grid1,
                &grid2,
                (pair.1, pair.0, pair.2),
                vdw1,
                vdw2,
                pbox,
                &mut found,
            );
            found
        })
        .flatten()
        .collect()
}

//-------------------------------------------------------------------------

pub fn distance_search_single<T, C>(
    cutoff: f32,
    data: &impl PosProvider,
    ids: impl Iterator<Item = usize>,
) -> C
where
    T: SearchOutputType + Send + Sync,
    C: FromIterator<T> + FromParallelIterator<T>,
{
    // Compute the extents
    let (lower, upper) = compute_bounding_box_single(cutoff, data);

    let mut grid = Grid::from_cutoff_and_min_max(cutoff, &lower, &upper);
    grid.populate(data.iter_pos(), ids, &lower, &upper);

    let plan = search_plan(&grid, None, PBC_NONE);

    // Cycle over search plan and perform search for each cell pair
    plan.into_par_iter()
        .with_min_len(3)
        .map(|pair| search_cell_pair_single(cutoff * cutoff, &grid, pair))
        .flatten()
        .collect()
}

pub fn distance_search_single_pbc<T, C>(
    cutoff: f32,
    data: &(impl PosProvider + BoxProvider + ?Sized),
    ids: impl Iterator<Item = usize>,
    pbox: &PeriodicBox,
    pbc_dims: PbcDims,
) -> C
where
    T: SearchOutputType + Send + Sync,
    C: FromIterator<T> + FromParallelIterator<T>,
{
    let mut grid = Grid::from_cutoff_and_box(cutoff, pbox);

    let mut buf = vec![];

    grid.populate_pbc(data.iter_pos(), ids, pbox, pbc_dims, &mut buf);

    let plan = search_plan(&grid, None, pbc_dims);

    // Cycle over search plan and perform search for each cell pair
    plan.into_par_iter()
        .with_min_len(3)
        .map(|pair| search_cell_pair_single_pbc(cutoff * cutoff, &grid, pair, pbox))
        .flatten()
        .collect()
}

//-------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use crate::core::Source;

    use super::WritableToFile;

    #[test]
    fn within_plan_test() -> anyhow::Result<()> {
        let src = Source::serial_from_file("tests/albumin.pdb")?;
        let t = std::time::Instant::now();
        src.select_str("within 4.0 of resid 10:300")?;
        println!("elapsed {}", t.elapsed().as_secs_f32());
        Ok(())
    }

    #[test]
    fn within_pbc() -> anyhow::Result<()> {
        let src = Source::serial_from_file("tests/albumin.pdb")?;
        let sel = src.select_str("within 2.0 pbc yyy of (resindex 16894 and name OW)")?;
        sel.save("target/pbc_sel.pdb")?;
        Ok(())
    }
}
