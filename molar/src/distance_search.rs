use crate::par::*;
use crate::prelude::*;
use num_traits::clamp_min;

/// Trait for the results of distance seacrh
pub trait DistanceSearchOutput {
    fn from_ijd(i: usize, j: usize, d: Float) -> Self;
}

impl DistanceSearchOutput for usize {
    fn from_ijd(i: usize, _j: usize, _d: Float) -> Self {
        i
    }
}

impl DistanceSearchOutput for (usize, usize) {
    fn from_ijd(i: usize, j: usize, _d: Float) -> Self {
        (i, j)
    }
}

impl DistanceSearchOutput for (usize, usize, Float) {
    fn from_ijd(i: usize, j: usize, d: Float) -> Self {
        (i, j, d)
    }
}

//--------------------------------------------------------------------------------

// For periodic selections Grid is self-referencial and thus
// must never be moved! Gris is only used
// inside search_* functions, so this should never be a problem.
struct Grid<'a> {
    cells: Vec<Vec<(usize, &'a Pos)>>,
    dims: [usize; 3],
    wrapped_pos: Vec<Pos>,
    fractional_lower: Vector3f,
    fractional_span: Vector3f,
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

    pub(crate) fn new_with_dims(dims: [usize; 3]) -> Self {
        Self {
            cells: vec![vec![]; dims[0] * dims[1] * dims[2]],
            dims,
            wrapped_pos: vec![],
            fractional_lower: Vector3f::zeros(),
            fractional_span: Vector3f::repeat(1.0),
        }
    }

    #[inline(always)]
    fn loc_to_ind(&self, loc: &[usize; 3]) -> usize {
        loc[0] + loc[1] * self.dims[0] + loc[2] * self.dims[0] * self.dims[1]
    }

    // pub fn get_loc_mut(&mut self, loc: [usize; 3]) -> &mut Vec<(usize,&'a Pos)> {
    //     let i = self.loc_to_ind(loc);
    //     &mut self.cells[i]
    // }

    fn push_loc(&mut self, loc: &[usize; 3], data: (usize, &'a Pos)) {
        let i = self.loc_to_ind(loc);
        self.cells[i].push(data);
    }

    unsafe fn push_ptr(&mut self, ind: usize, data: (usize, *const Pos)) {
        unsafe {
            self.cells[ind].push((data.0, &*data.1));
        }
    }

    fn from_cutoff_and_extents(cutoff: Float, extents: &Vector3f) -> Self {
        let mut sz = [0, 0, 0];
        // Cell size should be >= cutoff for all dimentions
        for d in 0..3 {
            sz[d] = clamp_min((extents[d] / cutoff).floor() as usize, 1);
        }
        Self::new_with_dims(sz)
    }

    pub(crate) fn from_cutoff_and_min_max(cutoff: Float, min: &Vector3f, max: &Vector3f) -> Self {
        Self::from_cutoff_and_extents(cutoff, &(max - min))
    }

    pub(crate) fn from_cutoff_and_box(cutoff: Float, box_: &PeriodicBox, pbc: PbcDims) -> Self {
        let mut spacings = box_.face_spacings();
        // Start nonperiodic directions with one bin. fit_nonperiodic expands
        // them to the data bounds before positions are inserted.
        for d in 0..3 {
            if !pbc.get_dim(d) {
                spacings[d] = 0.0;
            }
        }
        Self::from_cutoff_and_extents(cutoff, &spacings)
    }

    /// Set finite grid bounds from the data along nonperiodic directions.
    /// The unit-cell faces are not limits on nonperiodic coordinates.
    fn fit_nonperiodic<'b>(
        &mut self,
        cutoff: Float,
        box_: &PeriodicBox,
        pbc: PbcDims,
        points: impl Iterator<Item = &'b Pos>,
    ) {
        if pbc == PBC_FULL {
            return;
        }
        let mut lower = Vector3f::repeat(Float::INFINITY);
        let mut upper = Vector3f::repeat(Float::NEG_INFINITY);
        for p in points {
            let f = box_.to_box_coords(&p.coords);
            lower = lower.inf(&f);
            upper = upper.sup(&f);
        }
        let heights = box_.face_spacings();
        for d in 0..3 {
            if !pbc.get_dim(d) && lower[d].is_finite() {
                self.fractional_lower[d] = lower[d];
                // Include the upper endpoint in the final cell. Use one bin
                // for a flat coordinate distribution.
                let span = (upper[d] - lower[d]).max(cutoff / heights[d]);
                self.fractional_span[d] = span;
                self.dims[d] = ((span * heights[d] / cutoff).floor() as usize).max(1);
            }
        }
        self.cells = vec![vec![]; self.dims.iter().product()];
    }

    fn empty_like(&self) -> Self {
        let mut grid = Self::new_with_dims(self.dims);
        grid.fractional_lower = self.fractional_lower;
        grid.fractional_span = self.fractional_span;
        grid
    }

    pub(crate) fn populate(
        &mut self,
        data: impl Iterator<Item = &'a Pos>,
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
                let n = (self.dims[d] as Float * (pos[d] - lower[d]) / dim_sz[d]).floor() as isize;
                if n < 0 || n >= self.dims[d] as isize {
                    continue 'outer;
                } else {
                    loc[d] = n as usize;
                }
            }
            self.push_loc(&loc, (id, pos));
        }
    }

    pub(crate) fn populate_pbc(
        &mut self,
        data: impl Iterator<Item = &'a Pos>,
        ids: impl Iterator<Item = usize>,
        box_: &PeriodicBox,
        pbc_dims: PbcDims,
    ) {
        let mut wrapped_ind = vec![];
        self.wrapped_pos.clear();

        for (id, pos) in ids.zip(data) {
            let mut rel = box_.to_box_coords(&pos.coords);
            let mut loc = [0usize; 3];
            let mut shifted = false;
            for d in 0..3 {
                if pbc_dims.get_dim(d) {
                    if rel[d] < 0.0 || rel[d] >= 1.0 {
                        rel[d] -= rel[d].floor();
                        shifted = true;
                    }
                }
                let f = (rel[d] - self.fractional_lower[d]) / self.fractional_span[d];
                loc[d] = ((f * self.dims[d] as Float).floor() as usize).min(self.dims[d] - 1);
            }

            if shifted {
                self.wrapped_pos.push(Pos::from(box_.to_lab_coords(&rel)));
                wrapped_ind.push((self.loc_to_ind(&loc), id));
            } else {
                self.push_loc(&loc, (id, pos));
            }
        }

        // Add wrapped points to the grid if any
        for i in 0..wrapped_ind.len() {
            // We need to unsafely get a self-reference to wrapped_pos
            unsafe {
                let ptr = self.wrapped_pos.as_ptr().add(i);
                self.push_ptr(wrapped_ind[i].0, (wrapped_ind[i].1, ptr));
            }
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
                'mask: for (v1, v2) in &MASK {
                    let mut c = [
                        [x + v1[0], y + v1[1], z + v1[2]],
                        [x + v2[0], y + v2[1], z + v2[2]],
                    ];
                    // we only go to the right, so need to check the right edge
                    // Use the requested lattice for distance calculations,
                    // including within-cell pairs and grids with few bins.
                    for i in 0..=1 {
                        for d in 0..3 {
                            if c[i][d] == grid1.dims[d] {
                                if pbc_dims.get_dim(d) {
                                    c[i][d] = 0;
                                } else {
                                    // Drop point for non-periodic dimension
                                    continue 'mask;
                                }
                            }
                        }
                    }
                    // If we are here we need to add the cell pair to the plan
                    let i1 = grid1.loc_to_ind(&c[0]);
                    let i2 = grid1.loc_to_ind(&c[1]);

                    // Check if there are points in both cells
                    // This is different in single and double grid cases
                    if let Some(grid2) = grid2 {
                        if (grid1.cells[i1].len() > 0 && grid2.cells[i2].len() > 0)
                            || (grid2.cells[i1].len() > 0 && grid1.cells[i2].len() > 0)
                        {
                            plan.push((i1, i2, pbc_dims));
                        }
                    } else if grid1.cells[i1].len() > 0 && grid1.cells[i2].len() > 0 {
                        plan.push((i1, i2, pbc_dims));
                    }
                }
            }
        }
    }
    // With one or two bins, different stencil entries name the same pair.
    // Compare each unordered cell pair once, including the within-cell case.
    if pbc_dims.any() && grid1.dims.iter().any(|&n| n <= 2) {
        for (i, j, _) in &mut plan {
            if *i > *j {
                std::mem::swap(i, j);
            }
        }
        plan.sort_unstable_by_key(|&(i, j, _)| (i, j));
        plan.dedup_by_key(|pair| (pair.0, pair.1));
    }
    plan
}

fn search_cell_pair_within(
    cutoff2: Float,
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
    cutoff2: Float,
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

fn search_cell_pair_double<T: DistanceSearchOutput>(
    cutoff2: Float,
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

fn search_cell_pair_double_pbc<T: DistanceSearchOutput>(
    cutoff2: Float,
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

fn search_cell_pair_double_vdw<T: DistanceSearchOutput>(
    grid1: &Grid,
    grid2: &Grid,
    pair: (usize, usize, PbcDims),
    vdw1: &Vec<Float>,
    vdw2: &Vec<Float>,
    found: &mut Vec<T>,
) {
    let n1 = grid1.cells[pair.0].len();
    let n2 = grid2.cells[pair.1].len();

    for i in 0..n1 {
        let (ind1, pos1) = grid1.cells[pair.0][i];
        for j in 0..n2 {
            let (ind2, pos2) = grid2.cells[pair.1][j];

            let d2 = (pos2 - pos1).norm_squared();
            let cutoff = vdw1[ind1] + vdw2[ind2] + Float::EPSILON;
            if d2 <= cutoff * cutoff {
                found.push(T::from_ijd(ind1, ind2, d2.sqrt()));
            }
        }
    }
}

fn search_cell_pair_double_vdw_pbc<T: DistanceSearchOutput>(
    grid1: &Grid,
    grid2: &Grid,
    pair: (usize, usize, PbcDims),
    vdw1: &Vec<Float>,
    vdw2: &Vec<Float>,
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

            let cutoff = vdw1[ind1] + vdw2[ind2] + Float::EPSILON;

            if d2 <= cutoff * cutoff {
                found.push(T::from_ijd(ind1, ind2, d2.sqrt()));
            }
        }
    }
}

fn search_cell_pair_single<T: DistanceSearchOutput>(
    cutoff2: Float,
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

fn search_cell_pair_single_pbc<T: DistanceSearchOutput>(
    cutoff2: Float,
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
    cutoff: Float,
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
            search_cell_pair_within(cutoff * cutoff, &grid1, &grid2, pair, &mut found);
            if pair.0 != pair.1 {
                search_cell_pair_within(
                    cutoff * cutoff,
                    &grid1,
                    &grid2,
                    (pair.1, pair.0, pair.2),
                    &mut found,
                );
            }
            found
        })
        .flatten()
        .collect()
}

pub(crate) fn distance_search_within_pbc<C>(
    cutoff: Float,
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
    let mut grid1 = Grid::from_cutoff_and_box(cutoff, pbox, pbc_dims);
    grid1.fit_nonperiodic(
        cutoff,
        pbox,
        pbc_dims,
        data1.iter_pos().chain(data2.iter_pos()),
    );
    let mut grid2 = grid1.empty_like();

    grid1.populate_pbc(data1.iter_pos(), ids1, pbox, pbc_dims);
    grid2.populate_pbc(data2.iter_pos(), ids2, pbox, pbc_dims);

    let plan = search_plan(&grid1, Some(&grid2), pbc_dims);

    // Cycle over search plan and perform search for each cell pair
    plan.into_par_iter()
        .with_min_len(3)
        .map(|pair| {
            let mut found = Vec::new();
            search_cell_pair_within_pbc(cutoff * cutoff, &grid1, &grid2, pair, pbox, &mut found);
            if pair.0 != pair.1 {
                search_cell_pair_within_pbc(
                    cutoff * cutoff,
                    &grid1,
                    &grid2,
                    (pair.1, pair.0, pair.2),
                    pbox,
                    &mut found,
                );
            }
            found
        })
        .flatten()
        .collect()
}

//-------------------------------------------------------------------------

fn compute_min_max<'a>(data: impl Iterator<Item = &'a Pos>) -> (Vector3f, Vector3f) {
    let mut lower = Vector3f::zeros();
    let mut upper = Vector3f::zeros();
    for p in data {
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

fn compute_bounding_box_double<'a>(
    cutoff: Float,
    data1: impl Iterator<Item = &'a Pos>,
    data2: impl Iterator<Item = &'a Pos>,
) -> (Vector3f, Vector3f) {
    let (l1, u1) = compute_min_max(data1);
    let (l2, u2) = compute_min_max(data2);

    let mut l = Vector3f::zeros();
    let mut u = Vector3f::zeros();
    for d in 0..3 {
        l[d] = l1[d].min(l2[d]);
        u[d] = u1[d].max(u2[d]);
    }

    l.add_scalar_mut(-cutoff - Float::EPSILON);
    u.add_scalar_mut(cutoff + Float::EPSILON);
    (l, u)
}

fn compute_bounding_box_single<'a>(
    cutoff: Float,
    data: impl Iterator<Item = &'a Pos>,
) -> (Vector3f, Vector3f) {
    let (mut l, mut u) = compute_min_max(data);
    l.add_scalar_mut(-cutoff - Float::EPSILON);
    u.add_scalar_mut(cutoff + Float::EPSILON);
    (l, u)
}

/// Performs distance search between two sets of points within a given cutoff distance
///
/// # Arguments
/// * `cutoff` - Maximum distance between points to be considered neighbors
/// * `data1` - Iterator providing positions for first set of points
/// * `data2` - Iterator providing positions for second set of points
/// * `ids1` - Iterator providing indices for first set of points
/// * `ids2` - Iterator providing indices for second set of points
///
/// # Returns
/// Collection of matched elements as specified by type parameters T and C
pub fn distance_search_double<T, C>(
    cutoff: Float,
    data1: &impl PosProvider,
    data2: &impl PosProvider,
    ids1: impl Iterator<Item = usize>,
    ids2: impl Iterator<Item = usize>,
) -> C
where
    T: DistanceSearchOutput + Send + Sync,
    C: FromIterator<T> + FromParallelIterator<T>,
{
    // Compute the extents
    let (lower, upper) = compute_bounding_box_double(cutoff, data1.iter_pos(), data2.iter_pos());

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
            search_cell_pair_double(cutoff * cutoff, &grid1, &grid2, pair, &mut found);
            if pair.0 != pair.1 {
                search_cell_pair_double(
                    cutoff * cutoff,
                    &grid1,
                    &grid2,
                    (pair.1, pair.0, pair.2),
                    &mut found,
                );
            }
            found
        })
        .flatten()
        .collect()
}

/// Performs distance search between two sets of points within periodic boundary conditions
///
/// # Arguments
/// * `cutoff` - Maximum distance between points to be considered neighbors
/// * `data1` - Iterator providing positions for first set of points
/// * `data2` - Iterator providing positions for second set of points
/// * `ids1` - Iterator providing indices for first set of points
/// * `ids2` - Iterator providing indices for second set of points
/// * `pbox` - Periodic box definition
/// * `pbc_dims` - Enabled box vectors; other directions use the bounds of the data
///
/// # Returns
/// Collection of matched elements as specified by type parameters T and C
pub fn distance_search_double_pbc<'a, T, C>(
    cutoff: Float,
    data1: impl Iterator<Item = &'a Pos>,
    data2: impl Iterator<Item = &'a Pos>,
    ids1: impl Iterator<Item = usize>,
    ids2: impl Iterator<Item = usize>,
    pbox: &PeriodicBox,
    pbc_dims: PbcDims,
) -> C
where
    T: DistanceSearchOutput + Send + Sync,
    C: FromIterator<T> + FromParallelIterator<T>,
{
    let mut data1 = data1;
    let mut data2 = data2;
    let buffered1: Vec<_> = if pbc_dims != PBC_FULL {
        data1.by_ref().collect()
    } else {
        Vec::new()
    };
    let buffered2: Vec<_> = if pbc_dims != PBC_FULL {
        data2.by_ref().collect()
    } else {
        Vec::new()
    };
    let mut grid1 = Grid::from_cutoff_and_box(cutoff, pbox, pbc_dims);
    grid1.fit_nonperiodic(
        cutoff,
        pbox,
        pbc_dims,
        buffered1.iter().chain(&buffered2).copied(),
    );
    let mut grid2 = grid1.empty_like();
    let data1 = buffered1.into_iter().chain(data1);
    let data2 = buffered2.into_iter().chain(data2);

    grid1.populate_pbc(data1, ids1, pbox, pbc_dims);
    grid2.populate_pbc(data2, ids2, pbox, pbc_dims);
    // At this point grids are self-referencial. We should not
    // move or mutate it until it is dorpped.

    let plan = search_plan(&grid1, Some(&grid2), pbc_dims);

    // Cycle over search plan and perform search for each cell pair
    plan.into_par_iter()
        .with_min_len(3)
        .map(|pair| {
            let mut found = Vec::new();
            search_cell_pair_double_pbc(cutoff * cutoff, &grid1, &grid2, pair, pbox, &mut found);
            if pair.0 != pair.1 {
                search_cell_pair_double_pbc(
                    cutoff * cutoff,
                    &grid1,
                    &grid2,
                    (pair.1, pair.0, pair.2),
                    pbox,
                    &mut found,
                );
            }
            found
        })
        .flatten()
        .collect()
}

// This always returns local idexes
/// Performs distance search between two sets of points using van der Waals radii
///
/// # Arguments
/// * `data1` - Iterator providing positions for first set of points
/// * `data2` - Iterator providing positions for second set of points
/// * `vdw1` - Van der Waals radii for first set of points
/// * `vdw2` - Van der Waals radii for second set of points
///
/// # Returns
/// Collection of matched elements as specified by type parameters T and C
pub fn distance_search_double_vdw<'a, T, C>(
    data1: &impl PosProvider,
    data2: &impl PosProvider,
    vdw1: &Vec<Float>,
    vdw2: &Vec<Float>,
) -> C
where
    T: DistanceSearchOutput + Send + Sync,
    C: FromIterator<T> + FromParallelIterator<T>,
{
    // We need to find the largest VdW distance pair to get the grid extents
    // let cutoff = vdw1.par_iter().cloned().reduce(|| Float::MIN, Float::max)
    //     + vdw2.par_iter().cloned().reduce(|| Float::MIN, Float::max)
    //     + Float::EPSILON;
    let cutoff = vdw1.iter().cloned().reduce(Float::max).unwrap()
        + vdw2.iter().cloned().reduce(Float::max).unwrap()
        + Float::EPSILON;

    // Compute the extents
    let (lower, upper) = compute_bounding_box_double(cutoff, data1.iter_pos(), data2.iter_pos());

    let mut grid1 = Grid::from_cutoff_and_min_max(cutoff, &lower, &upper);
    let mut grid2 = Grid::new_with_dims(grid1.get_dims());

    grid1.populate(data1.iter_pos(), 0..vdw1.len(), &lower, &upper);
    grid2.populate(data2.iter_pos(), 0..vdw2.len(), &lower, &upper);

    let plan = search_plan(&grid1, Some(&grid2), PBC_NONE);

    // Cycle over search plan and perform search for each cell pair
    plan.into_par_iter()
        .with_min_len(3)
        .map(|pair| {
            let mut found = Vec::new();
            search_cell_pair_double_vdw(&grid1, &grid2, pair, vdw1, vdw2, &mut found);
            if pair.0 != pair.1 {
                search_cell_pair_double_vdw(
                    &grid1,
                    &grid2,
                    (pair.1, pair.0, pair.2),
                    vdw1,
                    vdw2,
                    &mut found,
                );
            }
            found
        })
        .flatten()
        .collect()
}

// This always returns local idexes
/// Performs distance search between two sets of points using van der Waals radii with periodic boundaries
///
/// # Arguments
/// * `data1` - Iterator providing positions for first set of points
/// * `data2` - Iterator providing positions for second set of points
/// * `vdw1` - Van der Waals radii for first set of points
/// * `vdw2` - Van der Waals radii for second set of points
/// * `pbox` - Periodic box definition
/// * `pbc_dims` - Enabled box vectors; other directions use the bounds of the data
///
/// # Returns
/// Collection of matched elements as specified by type parameters T and C
pub fn distance_search_double_vdw_pbc<'a, T, C>(
    data1: impl Iterator<Item = &'a Pos>,
    data2: impl Iterator<Item = &'a Pos>,
    vdw1: &Vec<Float>,
    vdw2: &Vec<Float>,
    pbox: &PeriodicBox,
    pbc_dims: PbcDims,
) -> C
where
    T: DistanceSearchOutput + Send + Sync,
    C: FromIterator<T> + FromParallelIterator<T>,
{
    // We need to find the largest VdW distance pair to get the grid extents
    // let cutoff = vdw1.par_iter().cloned().reduce(|| Float::MIN, Float::max)
    //     + vdw2.par_iter().cloned().reduce(|| Float::MIN, Float::max)
    //     + Float::EPSILON;
    let cutoff = vdw1.iter().cloned().reduce(Float::max).unwrap()
        + vdw2.iter().cloned().reduce(Float::max).unwrap()
        + Float::EPSILON;

    let mut data1 = data1;
    let mut data2 = data2;
    let buffered1: Vec<_> = if pbc_dims != PBC_FULL {
        data1.by_ref().collect()
    } else {
        Vec::new()
    };
    let buffered2: Vec<_> = if pbc_dims != PBC_FULL {
        data2.by_ref().collect()
    } else {
        Vec::new()
    };
    let mut grid1 = Grid::from_cutoff_and_box(cutoff, pbox, pbc_dims);
    grid1.fit_nonperiodic(
        cutoff,
        pbox,
        pbc_dims,
        buffered1.iter().chain(&buffered2).copied(),
    );
    let mut grid2 = grid1.empty_like();
    let data1 = buffered1.into_iter().chain(data1);
    let data2 = buffered2.into_iter().chain(data2);

    grid1.populate_pbc(data1, 0..vdw1.len(), pbox, pbc_dims);
    grid2.populate_pbc(data2, 0..vdw2.len(), pbox, pbc_dims);

    // At this point grids are self-referencial.
    // Now on we should not move or mutate it until it is dorpped!

    let plan = search_plan(&grid1, Some(&grid2), pbc_dims);

    // Cycle over search plan and perform search for each cell pair
    plan.into_par_iter()
        .with_min_len(3)
        .map(|pair| {
            let mut found = Vec::new();
            search_cell_pair_double_vdw_pbc(&grid1, &grid2, pair, vdw1, vdw2, pbox, &mut found);
            if pair.0 != pair.1 {
                search_cell_pair_double_vdw_pbc(
                    &grid1,
                    &grid2,
                    (pair.1, pair.0, pair.2),
                    vdw1,
                    vdw2,
                    pbox,
                    &mut found,
                );
            }
            found
        })
        .flatten()
        .collect()
}

//-------------------------------------------------------------------------

/// Performs distance search within a single set of points
///
/// # Arguments
/// * `cutoff` - Maximum distance between points to be considered neighbors
/// * `data` - Iterator providing positions for points
/// * `ids` - Iterator providing indices for points
///
/// # Returns
/// Collection of matched elements as specified by type parameters T and C
pub fn distance_search_single<T, C>(
    cutoff: Float,
    data: &impl PosProvider,
    ids: impl Iterator<Item = usize>,
) -> C
where
    T: DistanceSearchOutput + Send + Sync,
    C: FromIterator<T> + FromParallelIterator<T>,
{
    // Compute the extents
    let (lower, upper) = compute_bounding_box_single(cutoff, data.iter_pos());

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

/// Performs distance search within a single set of points with periodic boundaries
///
/// # Arguments
/// * `cutoff` - Maximum distance between points to be considered neighbors
/// * `data` - Iterator providing positions for points
/// * `ids` - Iterator providing indices for points
/// * `pbox` - Periodic box definition
/// * `pbc_dims` - Enabled box vectors; other directions use the bounds of the data
///
/// # Returns
/// Collection of matched elements as specified by type parameters T and C
pub fn distance_search_single_pbc<'a, T, C>(
    cutoff: Float,
    //data: &(impl PosProvider + ?Sized),
    data: impl Iterator<Item = &'a Pos>,
    ids: impl Iterator<Item = usize>,
    pbox: &PeriodicBox,
    pbc_dims: PbcDims,
) -> C
where
    T: DistanceSearchOutput + Send + Sync,
    C: FromIterator<T> + FromParallelIterator<T>,
{
    let mut data = data;
    let buffered: Vec<_> = if pbc_dims != PBC_FULL {
        data.by_ref().collect()
    } else {
        Vec::new()
    };
    let mut grid = Grid::from_cutoff_and_box(cutoff, pbox, pbc_dims);
    grid.fit_nonperiodic(cutoff, pbox, pbc_dims, buffered.iter().copied());
    let data = buffered.into_iter().chain(data);
    grid.populate_pbc(data, ids, pbox, pbc_dims);
    // At this point grid is self-referencial. We pin it on the stack
    // to ensure that we don't move or mutate it until it is dorpped.
    // Normally this should never happen, but better to be safe

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
    use super::SaveTopologyState;
    use crate::prelude::*;

    #[test]
    fn within_plan_test() -> anyhow::Result<()> {
        let src = System::from_file("tests/albumin.pdb")?;
        let t = web_time::Instant::now();
        src.select("within 4.0 of resid 10:300")?;
        println!("elapsed {}", t.elapsed().as_secs_f32());
        Ok(())
    }

    #[test]
    fn within_pbc() -> anyhow::Result<()> {
        let src = System::from_file("tests/albumin.pdb")?;
        let sel = src.select_bound("within 2.0 pbc yyy of (resindex 16894 and name OW)")?;
        sel.save("../target/pbc_sel.pdb")?;
        Ok(())
    }

    #[test]
    fn triclinic_grid_finds_pair_across_two_old_bins() {
        let b =
            PeriodicBox::from_matrix(Matrix3f::new(10., 4., 0., 0., 10., 0., 0., 0., 10.)).unwrap();
        let points = [Pos::new(2.7, 5., 5.), Pos::new(3.5, 5., 5.)];
        let pairs: Vec<(usize, usize)> =
            super::distance_search_single_pbc(1., points.iter(), 0..2, &b, PBC_FULL);
        assert_eq!(pairs.len(), 1);
        assert!(pairs[0] == (0, 1) || pairs[0] == (1, 0));
    }

    #[test]
    fn periodic_grids_match_direct_search_all_masks_and_cutoffs() {
        use std::collections::BTreeSet;
        let matrices = [
            Matrix3f::new(10., 4., -4., 0., 10., 2., 0., 0., 10.),
            Matrix3f::new(10., 0., -2., 0., 10., 0., 0., 0., 1.),
            Matrix3f::from_diagonal(&Vector3f::new(3., 4., 5.)),
        ];
        let mut seed = 317_u64;
        for m in matrices {
            let b = PeriodicBox::from_matrix(m).unwrap();
            let mut points: Vec<_> = (0..36)
                .map(|_| {
                    Pos::from(
                        m * Vector3f::from_fn(|_, _| {
                            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                            ((seed >> 32) as u32 as Float / u32::MAX as Float - 0.5) * 3.0
                        }),
                    )
                })
                .collect();
            // Boundary points, plus points far outside nonperiodic faces.
            points.extend([
                Pos::origin(),
                Pos::from(m.column(0).into_owned()),
                Pos::from(m * Vector3f::new(-1.01, 0.49, 4.0)),
                Pos::from(m * Vector3f::new(0.01, 0.51, 4.0)),
            ]);
            for mask in 0..8 {
                let pbc = PbcDims::new(mask & 1 != 0, mask & 2 != 0, mask & 4 != 0);
                for cutoff in [0.7, 2.6, 12.0] {
                    let mut expected = BTreeSet::new();
                    for i in 0..points.len() {
                        for j in i + 1..points.len() {
                            if b.distance_squared(&points[i], &points[j], pbc) <= cutoff * cutoff {
                                expected.insert((i, j));
                            }
                        }
                    }
                    let got: Vec<(usize, usize)> = super::distance_search_single_pbc(
                        cutoff,
                        points.iter(),
                        0..points.len(),
                        &b,
                        pbc,
                    );
                    let unique: BTreeSet<_> =
                        got.iter().map(|&(i, j)| (i.min(j), i.max(j))).collect();
                    assert_eq!(
                        got.len(),
                        unique.len(),
                        "duplicate single pairs, mask {mask}, cutoff {cutoff}"
                    );
                    assert_eq!(
                        unique, expected,
                        "single mask {mask}, cutoff {cutoff}, box {m:?}"
                    );

                    let n = points.len() / 2;
                    let expected: BTreeSet<_> = expected
                        .into_iter()
                        .filter(|&(i, j)| i < n && j >= n)
                        .collect();
                    let got: Vec<(usize, usize)> = super::distance_search_double_pbc(
                        cutoff,
                        points[..n].iter(),
                        points[n..].iter(),
                        0..n,
                        n..points.len(),
                        &b,
                        pbc,
                    );
                    assert_eq!(
                        got.len(),
                        expected.len(),
                        "double count mask {mask}, cutoff {cutoff}"
                    );
                    assert_eq!(got.into_iter().collect::<BTreeSet<_>>(), expected);

                    let radii = vec![cutoff * 0.5; n];
                    let got: Vec<(usize, usize)> = super::distance_search_double_vdw_pbc(
                        points[..n].iter(),
                        points[n..].iter(),
                        &radii,
                        &radii,
                        &b,
                        pbc,
                    );
                    assert_eq!(
                        got.len(),
                        expected.len(),
                        "vdw count mask {mask}, cutoff {cutoff}"
                    );
                    assert_eq!(
                        got.into_iter()
                            .map(|(i, j)| (i, j + n))
                            .collect::<BTreeSet<_>>(),
                        expected
                    );
                }
            }
        }
    }

    #[test]
    fn within_partial_pbc_keeps_points_outside_nonperiodic_faces() {
        let b =
            PeriodicBox::from_matrix(Matrix3f::new(10., 4., 0., 0., 10., 0., 0., 0., 10.)).unwrap();
        let pbc = PbcDims::new(true, true, false);
        let points = [Pos::new(-0.1, 0., 100.), Pos::new(10.1, 0., 100.)];
        let mut grid1 = super::Grid::from_cutoff_and_box(0.5, &b, pbc);
        let mut grid2 = super::Grid::from_cutoff_and_box(0.5, &b, pbc);
        grid1.populate_pbc(points[..1].iter(), 0..1, &b, pbc);
        grid2.populate_pbc(points[1..].iter(), 1..2, &b, pbc);
        let mut found = Vec::new();
        for pair in super::search_plan(&grid1, Some(&grid2), pbc) {
            super::search_cell_pair_within_pbc(0.25, &grid1, &grid2, pair, &b, &mut found);
            if pair.0 != pair.1 {
                super::search_cell_pair_within_pbc(
                    0.25,
                    &grid1,
                    &grid2,
                    (pair.1, pair.0, pair.2),
                    &b,
                    &mut found,
                );
            }
        }
        assert_eq!(found, vec![0]);
    }
}
