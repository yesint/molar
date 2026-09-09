use crate::par::*;
use crate::prelude::*;
use thiserror::Error;

/// Invalid input to a distance search.
#[derive(Debug, Error, PartialEq)]
pub enum DistanceSearchError {
    #[error("cutoff and its square must be finite and positive: {0}")]
    InvalidCutoff(Float),
    #[error("distance search coordinates and grid bounds must be finite")]
    InvalidCoordinates,
    #[error("distance search requires finite, positive box face spacings")]
    InvalidBox,
    #[error("position and index/radius counts differ")]
    LengthMismatch,
    #[error("van der Waals radius must be finite and nonnegative: {0}")]
    InvalidRadius(Float),
}

fn validate_cutoff(cutoff: Float) -> Result<(), DistanceSearchError> {
    let squared = cutoff * cutoff;
    if cutoff <= 0.0 || !squared.is_finite() || squared == 0.0 {
        return Err(DistanceSearchError::InvalidCutoff(cutoff));
    }
    Ok(())
}

fn validate_pos(pos: &Pos) -> Result<(), DistanceSearchError> {
    if pos.coords.iter().all(|v| v.is_finite()) {
        Ok(())
    } else {
        Err(DistanceSearchError::InvalidCoordinates)
    }
}

fn vdw_cutoff(vdw1: &[Float], vdw2: &[Float]) -> Result<Float, DistanceSearchError> {
    let mut maxima = [0.0 as Float; 2];
    for (radii, maximum) in [vdw1, vdw2].into_iter().zip(&mut maxima) {
        for &radius in radii {
            if !radius.is_finite() || radius < 0.0 {
                return Err(DistanceSearchError::InvalidRadius(radius));
            }
            *maximum = maximum.max(radius);
        }
    }
    let cutoff = maxima[0] + maxima[1] + Float::EPSILON;
    validate_cutoff(cutoff)?;
    Ok(cutoff)
}

// Bound dense-grid metadata even for tiny cutoffs or very sparse inputs.
// Reducing dimensions only enlarges cells, so the neighbor stencil stays valid.
const MAX_GRID_CELLS: usize = 1 << 20;

fn bounded_dims(mut dims: [usize; 3], limit: usize) -> [usize; 3] {
    let limit = limit.clamp(1, MAX_GRID_CELLS);
    for n in &mut dims {
        *n = (*n).clamp(1, limit);
    }
    while dims
        .iter()
        .try_fold(1usize, |n, &d| n.checked_mul(d))
        .is_none_or(|n| n > limit)
    {
        let d = (0..3).max_by_key(|&d| dims[d]).unwrap();
        dims[d] = dims[d].div_ceil(2);
    }
    dims
}

// An upper size hint also works for filtered iterators. Unknown sizes use
// the absolute limit; no extra position pass or coordinate copy is needed.
fn cell_limit(points: impl Iterator) -> usize {
    points
        .size_hint()
        .1
        .unwrap_or(MAX_GRID_CELLS)
        .saturating_mul(8)
        .clamp(1, MAX_GRID_CELLS)
}

/// Trait for the results of distance search
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

type CellEntry<'a> = (usize, &'a Pos);

// Entries in each cell form one contiguous slice. Offsets include an end sentinel.
struct Cells<'a> {
    entries: Vec<CellEntry<'a>>,
    offsets: Vec<usize>,
}

impl<'a> Cells<'a> {
    fn new(n: usize) -> Self {
        Self {
            entries: Vec::new(),
            offsets: vec![0; n + 1],
        }
    }

    fn len(&self) -> usize {
        self.offsets.len() - 1
    }
}

impl<'a> std::ops::Index<usize> for Cells<'a> {
    type Output = [CellEntry<'a>];
    #[inline]
    fn index(&self, cell: usize) -> &Self::Output {
        &self.entries[self.offsets[cell]..self.offsets[cell + 1]]
    }
}

// Only owned values enter this cache. References to input coordinates never do.
#[derive(Default)]
struct GridScratch {
    wrapped_pos: Vec<Pos>,
    wrapped_entries: Vec<usize>,
    destinations: Vec<usize>,
    cursors: Vec<usize>,
}

const MAX_CACHED_SCRATCH_BYTES: usize = 8 * 1024 * 1024;
std::thread_local! {
    // At most two grids and 16 MiB of retained capacity per calling thread.
    static GRID_SCRATCH: std::cell::RefCell<[Option<GridScratch>; 2]> =
        const { std::cell::RefCell::new([None, None]) };
}

impl GridScratch {
    fn take_cached() -> Self {
        GRID_SCRATCH
            .try_with(|cache| {
                let mut cache = cache.try_borrow_mut().ok()?;
                cache.iter_mut().find_map(Option::take)
            })
            .ok()
            .flatten()
            .unwrap_or_default()
    }

    fn recycle(mut self) {
        let bytes = self
            .wrapped_pos
            .capacity()
            .saturating_mul(std::mem::size_of::<Pos>())
            .saturating_add(
                self.wrapped_entries
                    .capacity()
                    .saturating_mul(std::mem::size_of::<usize>()),
            )
            .saturating_add(
                self.destinations
                    .capacity()
                    .saturating_mul(std::mem::size_of::<usize>()),
            )
            .saturating_add(
                self.cursors
                    .capacity()
                    .saturating_mul(std::mem::size_of::<usize>()),
            );
        if bytes == 0 || bytes > MAX_CACHED_SCRATCH_BYTES {
            return;
        }
        self.wrapped_pos.clear();
        self.wrapped_entries.clear();
        self.destinations.clear();
        self.cursors.clear();
        // A nested search or thread teardown must not cause a destructor panic.
        let _ = GRID_SCRATCH.try_with(|cache| {
            if let Ok(mut cache) = cache.try_borrow_mut() {
                if let Some(slot) = cache.iter_mut().find(|slot| slot.is_none()) {
                    *slot = Some(self);
                }
            }
        });
    }
}

// Cells may refer to scratch.wrapped_pos. This buffer is frozen after population.
// Drop clears those references before returning the owned scratch buffers.
struct Grid<'a> {
    cells: Cells<'a>,
    dims: [usize; 3],
    cell_limit: usize,
    scratch: GridScratch,
    fractional_lower: Vector3f,
    fractional_span: Vector3f,
}

impl Drop for Grid<'_> {
    fn drop(&mut self) {
        self.cells.entries.clear();
        std::mem::take(&mut self.scratch).recycle();
    }
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
    pub(crate) fn new_with_dims(dims: [usize; 3]) -> Self {
        let dims = bounded_dims(dims, MAX_GRID_CELLS);
        Self {
            cells: Cells::new(dims[0] * dims[1] * dims[2]),
            dims,
            cell_limit: MAX_GRID_CELLS,
            scratch: GridScratch::take_cached(),
            fractional_lower: Vector3f::zeros(),
            fractional_span: Vector3f::repeat(1.0),
        }
    }

    #[inline(always)]
    fn loc_to_ind(&self, loc: &[usize; 3]) -> usize {
        loc[0] + loc[1] * self.dims[0] + loc[2] * self.dims[0] * self.dims[1]
    }

    fn from_cutoff_and_extents(
        cutoff: Float,
        extents: &Vector3f,
        cell_limit: usize,
    ) -> Result<Self, DistanceSearchError> {
        validate_cutoff(cutoff)?;
        if extents.iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(DistanceSearchError::InvalidCoordinates);
        }
        let mut sz = [0, 0, 0];
        // Cell size should be >= cutoff for all dimentions
        for d in 0..3 {
            sz[d] = ((extents[d] / cutoff).floor() as usize).max(1);
        }
        let mut grid = Self::new_with_dims(bounded_dims(sz, cell_limit));
        grid.cell_limit = cell_limit.clamp(1, MAX_GRID_CELLS);
        Ok(grid)
    }

    pub(crate) fn from_cutoff_and_min_max(
        cutoff: Float,
        min: &Vector3f,
        max: &Vector3f,
        cell_limit: usize,
    ) -> Result<Self, DistanceSearchError> {
        Self::from_cutoff_and_extents(cutoff, &(max - min), cell_limit)
    }

    pub(crate) fn from_cutoff_and_box(
        cutoff: Float,
        box_: &PeriodicBox,
        pbc: PbcDims,
        cell_limit: usize,
    ) -> Result<Self, DistanceSearchError> {
        let mut spacings = box_.face_spacings();
        if spacings.iter().any(|v| !v.is_finite() || *v <= 0.0) {
            return Err(DistanceSearchError::InvalidBox);
        }
        // Start nonperiodic directions with one bin. fit_nonperiodic expands
        // them to the data bounds before positions are inserted.
        for d in 0..3 {
            if !pbc.get_dim(d) {
                spacings[d] = 0.0;
            }
        }
        Self::from_cutoff_and_extents(cutoff, &spacings, cell_limit)
    }

    /// Set finite grid bounds from the data along nonperiodic directions.
    /// The unit-cell faces are not limits on nonperiodic coordinates.
    fn fit_nonperiodic(
        &mut self,
        cutoff: Float,
        box_: &PeriodicBox,
        lower: Vector3f,
        upper: Vector3f,
        pbc: PbcDims,
    ) -> Result<(), DistanceSearchError> {
        let heights = box_.face_spacings();
        for d in 0..3 {
            if !pbc.get_dim(d) && lower[d].is_finite() {
                self.fractional_lower[d] = lower[d];
                // Include the upper endpoint in the final cell. Use one bin
                // for a flat coordinate distribution.
                let span = (upper[d] - lower[d]).max(cutoff / heights[d]);
                if !span.is_finite() {
                    return Err(DistanceSearchError::InvalidCoordinates);
                }
                self.fractional_span[d] = span;
                self.dims[d] = ((span * heights[d] / cutoff).floor() as usize).max(1);
            }
        }
        self.dims = bounded_dims(self.dims, self.cell_limit);
        self.cells
            .offsets
            .resize(self.dims.iter().product::<usize>() + 1, 0);
        Ok(())
    }

    fn empty_like(&self) -> Self {
        let mut grid = Self::new_with_dims(self.dims);
        grid.cell_limit = self.cell_limit;
        grid.fractional_lower = self.fractional_lower;
        grid.fractional_span = self.fractional_span;
        grid
    }

    // Load once into the final entry array. Partial PBC uses these references to
    // find its bounds before binning, so it needs no separate input buffer.
    fn load(
        &mut self,
        data: impl Iterator<Item = &'a Pos>,
        mut ids: impl Iterator<Item = usize>,
    ) -> Result<(), DistanceSearchError> {
        debug_assert!(self.cells.entries.is_empty());
        self.cells
            .entries
            .reserve(data.size_hint().0.min(MAX_GRID_CELLS));
        for pos in data {
            let id = ids.next().ok_or(DistanceSearchError::LengthMismatch)?;
            validate_pos(pos)?;
            self.cells.entries.push((id, pos));
        }
        if ids.next().is_some() {
            return Err(DistanceSearchError::LengthMismatch);
        }
        Ok(())
    }

    // Stable counting permutation, performed in place. Only references and IDs
    // move. The destinations and cell cursors are reused by subsequent searches.
    fn arrange_cells(&mut self) {
        let offsets = &mut self.cells.offsets;
        offsets.fill(0);
        for &cell in &self.scratch.destinations {
            offsets[cell + 1] += 1;
        }
        for i in 1..offsets.len() {
            offsets[i] += offsets[i - 1];
        }
        self.scratch.cursors.clear();
        self.scratch
            .cursors
            .extend_from_slice(&offsets[..offsets.len() - 1]);
        for destination in &mut self.scratch.destinations {
            let cursor = &mut self.scratch.cursors[*destination];
            *destination = *cursor;
            *cursor += 1;
        }
        for i in 0..self.cells.entries.len() {
            while self.scratch.destinations[i] != i {
                let j = self.scratch.destinations[i];
                self.cells.entries.swap(i, j);
                self.scratch.destinations.swap(i, j);
            }
        }
    }

    pub(crate) fn populate(
        &mut self,
        data: impl Iterator<Item = &'a Pos>,
        ids: impl Iterator<Item = usize>,
        lower: &Vector3f,
        upper: &Vector3f,
    ) -> Result<(), DistanceSearchError> {
        self.load(data, ids)?;
        let dim_sz = upper - lower;
        let dims = self.dims;
        let destinations = &mut self.scratch.destinations;
        destinations.clear();
        destinations.reserve(self.cells.entries.len());
        self.cells.entries.retain(|&(_, pos)| {
            let mut loc = [0; 3];
            for d in 0..3 {
                if pos[d] < lower[d] || pos[d] > upper[d] {
                    return false;
                }
                let f = if dim_sz[d] == 0.0 {
                    0.0
                } else {
                    (pos[d] - lower[d]) / dim_sz[d]
                };
                loc[d] = ((f * dims[d] as Float).floor() as usize).min(dims[d] - 1);
            }
            destinations.push(loc[0] + dims[0] * (loc[1] + dims[1] * loc[2]));
            true
        });
        self.arrange_cells();
        Ok(())
    }

    fn bin_pbc(&mut self, box_: &PeriodicBox, pbc: PbcDims) -> Result<(), DistanceSearchError> {
        self.scratch.destinations.clear();
        self.scratch.destinations.reserve(self.cells.entries.len());
        for i in 0..self.cells.entries.len() {
            let pos = self.cells.entries[i].1;
            let mut rel = box_.to_box_coords(&pos.coords);
            validate_pos(&Pos::from(rel))?;
            let mut loc = [0; 3];
            let mut shifted = false;
            for d in 0..3 {
                if pbc.get_dim(d) && (rel[d] < 0.0 || rel[d] >= 1.0) {
                    rel[d] -= rel[d].floor();
                    shifted = true;
                }
                let f = (rel[d] - self.fractional_lower[d]) / self.fractional_span[d];
                loc[d] = ((f * self.dims[d] as Float).floor() as usize).min(self.dims[d] - 1);
            }
            self.scratch.destinations.push(self.loc_to_ind(&loc));
            if shifted {
                let wrapped = Pos::from(box_.to_lab_coords(&rel));
                validate_pos(&wrapped)?;
                self.scratch.wrapped_pos.push(wrapped);
                self.scratch.wrapped_entries.push(i);
            }
        }
        // All wrapped coordinates now have stable addresses. These references
        // remain private to this grid. Drop clears them before caching the buffer.
        for (wrapped, &entry) in self.scratch.wrapped_entries.iter().enumerate() {
            self.cells.entries[entry].1 =
                unsafe { &*self.scratch.wrapped_pos.as_ptr().add(wrapped) };
        }
        self.arrange_cells();
        Ok(())
    }

    pub fn get_dims(&self) -> [usize; 3] {
        self.dims
    }
}

fn fractional_bounds<'a>(
    data: impl Iterator<Item = &'a Pos>,
    box_: &PeriodicBox,
) -> Result<(Vector3f, Vector3f), DistanceSearchError> {
    let mut lower = Vector3f::repeat(Float::INFINITY);
    let mut upper = Vector3f::repeat(Float::NEG_INFINITY);
    for pos in data {
        let f = box_.to_box_coords(&pos.coords);
        validate_pos(&Pos::from(f))?;
        lower = lower.inf(&f);
        upper = upper.sup(&f);
    }
    Ok((lower, upper))
}

fn periodic_grid<'a>(
    cutoff: Float,
    mut data: impl Iterator<Item = &'a Pos>,
    ids: impl Iterator<Item = usize>,
    pbox: &PeriodicBox,
    pbc: PbcDims,
) -> Result<Grid<'a>, DistanceSearchError> {
    let mut grid = Grid::from_cutoff_and_box(cutoff, pbox, pbc, cell_limit(data.by_ref()))?;
    // The iterator's upper hint bounds storage even when it is filtered.
    grid.load(data, ids)?;
    if pbc != PBC_FULL {
        let (lower, upper) = fractional_bounds(grid.cells.entries.iter().map(|&(_, p)| p), pbox)?;
        grid.fit_nonperiodic(cutoff, pbox, lower, upper, pbc)?;
    }
    grid.bin_pbc(pbox, pbc)?;
    Ok(grid)
}

fn periodic_grids<'a>(
    cutoff: Float,
    mut data1: impl Iterator<Item = &'a Pos>,
    mut data2: impl Iterator<Item = &'a Pos>,
    ids1: impl Iterator<Item = usize>,
    ids2: impl Iterator<Item = usize>,
    pbox: &PeriodicBox,
    pbc: PbcDims,
) -> Result<(Grid<'a>, Grid<'a>), DistanceSearchError> {
    let limit = cell_limit(data1.by_ref().chain(data2.by_ref()));
    let mut grid1 = Grid::from_cutoff_and_box(cutoff, pbox, pbc, limit)?;
    let mut grid2 = grid1.empty_like();
    grid1.load(data1, ids1)?;
    grid2.load(data2, ids2)?;
    if pbc != PBC_FULL {
        let (lower, upper) = fractional_bounds(
            grid1
                .cells
                .entries
                .iter()
                .chain(&grid2.cells.entries)
                .map(|&(_, p)| p),
            pbox,
        )?;
        grid1.fit_nonperiodic(cutoff, pbox, lower, upper, pbc)?;
        grid2.dims = grid1.dims;
        grid2.fractional_lower = grid1.fractional_lower;
        grid2.fractional_span = grid1.fractional_span;
        grid2.cells.offsets.resize(grid1.cells.offsets.len(), 0);
    }
    grid1.bin_pbc(pbox, pbc)?;
    grid2.bin_pbc(pbox, pbc)?;
    Ok((grid1, grid2))
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

// Reuse a result buffer across several cell pairs. Small plans and wasm
// avoid thread scheduling; native large plans use one buffer per Rayon fold.
fn collect_search<T, C>(count: usize, search: impl Fn(usize, &mut Vec<T>) + Send + Sync) -> C
where
    T: Send,
    C: FromIterator<T> + FromParallelIterator<T>,
{
    #[cfg(not(target_arch = "wasm32"))]
    if count > 32 {
        return (0..count)
            .into_par_iter()
            .fold(Vec::new, |mut found, pair| {
                search(pair, &mut found);
                found
            })
            .flatten()
            .collect();
    }
    let mut found = Vec::new();
    for pair in 0..count {
        search(pair, &mut found);
    }
    found.into_iter().collect()
}

fn collect_cell_pairs<T, C>(
    plan: Vec<(usize, usize, PbcDims)>,
    search: impl Fn((usize, usize, PbcDims), &mut Vec<T>) + Send + Sync,
) -> C
where
    T: Send,
    C: FromIterator<T> + FromParallelIterator<T>,
{
    collect_search(plan.len(), |i, found| search(plan[i], found))
}

// Center first, then distinct adjacent cells. One- and two-bin periodic axes
// must not visit the same cell twice.
fn adjacent_axis(center: usize, n: usize, periodic: bool) -> ([usize; 3], usize) {
    let mut cells = [center; 3];
    let mut count = 1;
    if center > 0 {
        cells[count] = center - 1;
        count += 1;
    } else if periodic && n > 1 {
        cells[count] = n - 1;
        count += 1;
    }
    let next = if center + 1 < n {
        Some(center + 1)
    } else if periodic {
        Some(0)
    } else {
        None
    };
    if let Some(next) = next {
        if !cells[..count].contains(&next) {
            cells[count] = next;
            count += 1;
        }
    }
    (cells, count)
}

fn search_within<C>(
    grid1: &Grid,
    grid2: &Grid,
    pbc: PbcDims,
    within: impl Fn(&Pos, &Pos) -> bool + Send + Sync,
) -> C
where
    C: FromIterator<usize> + FromParallelIterator<usize>,
{
    collect_search(grid1.cells.len(), |cell, found| {
        let query = &grid1.cells[cell];
        if query.is_empty() {
            return;
        }
        let x = cell % grid1.dims[0];
        let yz = cell / grid1.dims[0];
        let (xs, nx) = adjacent_axis(x, grid1.dims[0], pbc.get_dim(0));
        let (ys, ny) = adjacent_axis(yz % grid1.dims[1], grid1.dims[1], pbc.get_dim(1));
        let (zs, nz) = adjacent_axis(yz / grid1.dims[1], grid1.dims[2], pbc.get_dim(2));
        // Stack storage replaces the within search plan. Build this list once
        // per query cell, then stop each point at its first match in any cell.
        let mut neighbors = [0; 27];
        let mut count = 0;
        for &z in &zs[..nz] {
            for &y in &ys[..ny] {
                for &x in &xs[..nx] {
                    let neighbor = grid2.loc_to_ind(&[x, y, z]);
                    if !grid2.cells[neighbor].is_empty() {
                        neighbors[count] = neighbor;
                        count += 1;
                    }
                }
            }
        }
        'points: for &(id, pos) in query {
            for &neighbor in &neighbors[..count] {
                for &(_, target) in &grid2.cells[neighbor] {
                    if within(pos, target) {
                        found.push(id);
                        continue 'points;
                    }
                }
            }
        }
    })
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
    vdw1: &[Float],
    vdw2: &[Float],
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
    vdw1: &[Float],
    vdw2: &[Float],
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
    found: &mut Vec<T>,
) {
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
}

fn search_cell_pair_single_pbc<T: DistanceSearchOutput>(
    cutoff2: Float,
    grid: &Grid,
    pair: (usize, usize, PbcDims),
    pbox: &PeriodicBox,
    found: &mut Vec<T>,
) {
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
}

pub(crate) fn distance_search_within<'a, C>(
    cutoff: Float,
    data1: &impl PosProvider,
    data2: &impl PosProvider,
    ids1: impl Iterator<Item = usize>,
    ids2: impl Iterator<Item = usize>,
    lower: &Vector3f,
    upper: &Vector3f,
) -> Result<C, DistanceSearchError>
where
    C: FromIterator<usize> + FromParallelIterator<usize>,
{
    let limit = cell_limit(data1.iter_pos().chain(data2.iter_pos()));
    let mut grid1 = Grid::from_cutoff_and_min_max(cutoff, lower, upper, limit)?;
    let mut grid2 = Grid::new_with_dims(grid1.get_dims());

    grid1.populate(data1.iter_pos(), ids1, lower, upper)?;
    grid2.populate(data2.iter_pos(), ids2, lower, upper)?;

    //grid1.debug();

    Ok(search_within(&grid1, &grid2, PBC_NONE, |pos, target| {
        (pos - target).norm_squared() <= cutoff * cutoff
    }))
}

pub(crate) fn distance_search_within_pbc<C>(
    cutoff: Float,
    data1: &impl PosProvider,
    data2: &impl PosProvider,
    ids1: impl Iterator<Item = usize>,
    ids2: impl Iterator<Item = usize>,
    pbox: &PeriodicBox,
    pbc_dims: PbcDims,
) -> Result<C, DistanceSearchError>
where
    C: FromIterator<usize> + FromParallelIterator<usize>,
{
    let (grid1, grid2) = periodic_grids(
        cutoff,
        data1.iter_pos(),
        data2.iter_pos(),
        ids1,
        ids2,
        pbox,
        pbc_dims,
    )?;

    Ok(search_within(&grid1, &grid2, pbc_dims, |pos, target| {
        pbox.distance_squared(pos, target, pbc_dims) <= cutoff * cutoff
    }))
}

//-------------------------------------------------------------------------

fn compute_min_max<'a>(
    mut data: impl Iterator<Item = &'a Pos>,
) -> Result<(Vector3f, Vector3f), DistanceSearchError> {
    let Some(first) = data.next() else {
        return Ok((Vector3f::zeros(), Vector3f::zeros()));
    };
    validate_pos(first)?;
    let mut lower = first.coords;
    let mut upper = first.coords;
    for p in data {
        validate_pos(p)?;
        for d in 0..3 {
            if p[d] < lower[d] {
                lower[d] = p[d];
            }
            if p[d] > upper[d] {
                upper[d] = p[d];
            }
        }
    }
    Ok((lower, upper))
}

fn compute_bounding_box_double<'a>(
    cutoff: Float,
    data1: impl Iterator<Item = &'a Pos>,
    data2: impl Iterator<Item = &'a Pos>,
) -> Result<(Vector3f, Vector3f), DistanceSearchError> {
    validate_cutoff(cutoff)?;
    let (mut l, mut u) = compute_min_max(data1.chain(data2))?;

    l.add_scalar_mut(-cutoff - Float::EPSILON);
    u.add_scalar_mut(cutoff + Float::EPSILON);
    Ok((l, u))
}

fn compute_bounding_box_single<'a>(
    cutoff: Float,
    data: impl Iterator<Item = &'a Pos>,
) -> Result<(Vector3f, Vector3f), DistanceSearchError> {
    validate_cutoff(cutoff)?;
    let (mut l, mut u) = compute_min_max(data)?;
    l.add_scalar_mut(-cutoff - Float::EPSILON);
    u.add_scalar_mut(cutoff + Float::EPSILON);
    Ok((l, u))
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
/// Collection of matched elements as specified by type parameters T and C.
/// Empty inputs return an empty collection. Each ID/radius iterator must have
/// exactly one entry per position. Radii are nonnegative and use nanometers.
///
/// # Errors
/// Returns [`DistanceSearchError`] for invalid cutoffs, nonfinite coordinates,
/// invalid periodic boxes, invalid radii, or unequal input lengths. A numeric
/// cutoff and its square must both be finite and positive.
pub fn distance_search_double<T, C>(
    cutoff: Float,
    data1: &impl PosProvider,
    data2: &impl PosProvider,
    ids1: impl Iterator<Item = usize>,
    ids2: impl Iterator<Item = usize>,
) -> Result<C, DistanceSearchError>
where
    T: DistanceSearchOutput + Send + Sync,
    C: FromIterator<T> + FromParallelIterator<T>,
{
    let limit = cell_limit(data1.iter_pos().chain(data2.iter_pos()));
    // Compute the extents
    let (lower, upper) = compute_bounding_box_double(cutoff, data1.iter_pos(), data2.iter_pos())?;

    let mut grid1 = Grid::from_cutoff_and_min_max(cutoff, &lower, &upper, limit)?;
    let mut grid2 = Grid::new_with_dims(grid1.get_dims());

    grid1.populate(data1.iter_pos(), ids1, &lower, &upper)?;
    grid2.populate(data2.iter_pos(), ids2, &lower, &upper)?;

    let plan = search_plan(&grid1, Some(&grid2), PBC_NONE);

    // Cycle over search plan and perform search for each cell pair
    Ok(collect_cell_pairs(plan, |pair, found| {
        search_cell_pair_double(cutoff * cutoff, &grid1, &grid2, pair, found);
        if pair.0 != pair.1 {
            search_cell_pair_double(
                cutoff * cutoff,
                &grid1,
                &grid2,
                (pair.1, pair.0, pair.2),
                found,
            );
        }
    }))
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
/// Collection of matched elements as specified by type parameters T and C.
/// Empty inputs return an empty collection. Each ID/radius iterator must have
/// exactly one entry per position. Radii are nonnegative and use nanometers.
///
/// # Errors
/// Returns [`DistanceSearchError`] for invalid cutoffs, nonfinite coordinates,
/// invalid periodic boxes, invalid radii, or unequal input lengths. A numeric
/// cutoff and its square must both be finite and positive.
pub fn distance_search_double_pbc<'a, T, C>(
    cutoff: Float,
    data1: impl Iterator<Item = &'a Pos>,
    data2: impl Iterator<Item = &'a Pos>,
    ids1: impl Iterator<Item = usize>,
    ids2: impl Iterator<Item = usize>,
    pbox: &PeriodicBox,
    pbc_dims: PbcDims,
) -> Result<C, DistanceSearchError>
where
    T: DistanceSearchOutput + Send + Sync,
    C: FromIterator<T> + FromParallelIterator<T>,
{
    let (grid1, grid2) = periodic_grids(cutoff, data1, data2, ids1, ids2, pbox, pbc_dims)?;

    let plan = search_plan(&grid1, Some(&grid2), pbc_dims);

    // Cycle over search plan and perform search for each cell pair
    Ok(collect_cell_pairs(plan, |pair, found| {
        search_cell_pair_double_pbc(cutoff * cutoff, &grid1, &grid2, pair, pbox, found);
        if pair.0 != pair.1 {
            search_cell_pair_double_pbc(
                cutoff * cutoff,
                &grid1,
                &grid2,
                (pair.1, pair.0, pair.2),
                pbox,
                found,
            );
        }
    }))
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
/// Collection of matched elements as specified by type parameters T and C.
/// Empty inputs return an empty collection. Each ID/radius iterator must have
/// exactly one entry per position. Radii are nonnegative and use nanometers.
///
/// # Errors
/// Returns [`DistanceSearchError`] for invalid cutoffs, nonfinite coordinates,
/// invalid periodic boxes, invalid radii, or unequal input lengths. A numeric
/// cutoff and its square must both be finite and positive.
pub fn distance_search_double_vdw<'a, T, C>(
    data1: &impl PosProvider,
    data2: &impl PosProvider,
    vdw1: &[Float],
    vdw2: &[Float],
) -> Result<C, DistanceSearchError>
where
    T: DistanceSearchOutput + Send + Sync,
    C: FromIterator<T> + FromParallelIterator<T>,
{
    let limit = cell_limit(data1.iter_pos().chain(data2.iter_pos()));
    let cutoff = vdw_cutoff(vdw1, vdw2)?;

    // Compute the extents
    let (lower, upper) = compute_bounding_box_double(cutoff, data1.iter_pos(), data2.iter_pos())?;

    let mut grid1 = Grid::from_cutoff_and_min_max(cutoff, &lower, &upper, limit)?;
    let mut grid2 = Grid::new_with_dims(grid1.get_dims());

    grid1.populate(data1.iter_pos(), 0..vdw1.len(), &lower, &upper)?;
    grid2.populate(data2.iter_pos(), 0..vdw2.len(), &lower, &upper)?;

    let plan = search_plan(&grid1, Some(&grid2), PBC_NONE);

    // Cycle over search plan and perform search for each cell pair
    Ok(collect_cell_pairs(plan, |pair, found| {
        search_cell_pair_double_vdw(&grid1, &grid2, pair, vdw1, vdw2, found);
        if pair.0 != pair.1 {
            search_cell_pair_double_vdw(
                &grid1,
                &grid2,
                (pair.1, pair.0, pair.2),
                vdw1,
                vdw2,
                found,
            );
        }
    }))
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
/// Collection of matched elements as specified by type parameters T and C.
/// Empty inputs return an empty collection. Each ID/radius iterator must have
/// exactly one entry per position. Radii are nonnegative and use nanometers.
///
/// # Errors
/// Returns [`DistanceSearchError`] for invalid cutoffs, nonfinite coordinates,
/// invalid periodic boxes, invalid radii, or unequal input lengths. A numeric
/// cutoff and its square must both be finite and positive.
pub fn distance_search_double_vdw_pbc<'a, T, C>(
    data1: impl Iterator<Item = &'a Pos>,
    data2: impl Iterator<Item = &'a Pos>,
    vdw1: &[Float],
    vdw2: &[Float],
    pbox: &PeriodicBox,
    pbc_dims: PbcDims,
) -> Result<C, DistanceSearchError>
where
    T: DistanceSearchOutput + Send + Sync,
    C: FromIterator<T> + FromParallelIterator<T>,
{
    let cutoff = vdw_cutoff(vdw1, vdw2)?;
    let (grid1, grid2) = periodic_grids(
        cutoff,
        data1,
        data2,
        0..vdw1.len(),
        0..vdw2.len(),
        pbox,
        pbc_dims,
    )?;

    let plan = search_plan(&grid1, Some(&grid2), pbc_dims);

    // Cycle over search plan and perform search for each cell pair
    Ok(collect_cell_pairs(plan, |pair, found| {
        search_cell_pair_double_vdw_pbc(&grid1, &grid2, pair, vdw1, vdw2, pbox, found);
        if pair.0 != pair.1 {
            search_cell_pair_double_vdw_pbc(
                &grid1,
                &grid2,
                (pair.1, pair.0, pair.2),
                vdw1,
                vdw2,
                pbox,
                found,
            );
        }
    }))
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
/// Collection of matched elements as specified by type parameters T and C.
/// Empty inputs return an empty collection. Each ID/radius iterator must have
/// exactly one entry per position. Radii are nonnegative and use nanometers.
///
/// # Errors
/// Returns [`DistanceSearchError`] for invalid cutoffs, nonfinite coordinates,
/// invalid periodic boxes, invalid radii, or unequal input lengths. A numeric
/// cutoff and its square must both be finite and positive.
pub fn distance_search_single<T, C>(
    cutoff: Float,
    data: &impl PosProvider,
    ids: impl Iterator<Item = usize>,
) -> Result<C, DistanceSearchError>
where
    T: DistanceSearchOutput + Send + Sync,
    C: FromIterator<T> + FromParallelIterator<T>,
{
    let limit = cell_limit(data.iter_pos());
    // Compute the extents
    let (lower, upper) = compute_bounding_box_single(cutoff, data.iter_pos())?;

    let mut grid = Grid::from_cutoff_and_min_max(cutoff, &lower, &upper, limit)?;
    grid.populate(data.iter_pos(), ids, &lower, &upper)?;

    let plan = search_plan(&grid, None, PBC_NONE);

    // Cycle over search plan and perform search for each cell pair
    Ok(collect_cell_pairs(plan, |pair, found| {
        search_cell_pair_single(cutoff * cutoff, &grid, pair, found);
    }))
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
/// Collection of matched elements as specified by type parameters T and C.
/// Empty inputs return an empty collection. Each ID/radius iterator must have
/// exactly one entry per position. Radii are nonnegative and use nanometers.
///
/// # Errors
/// Returns [`DistanceSearchError`] for invalid cutoffs, nonfinite coordinates,
/// invalid periodic boxes, invalid radii, or unequal input lengths. A numeric
/// cutoff and its square must both be finite and positive.
pub fn distance_search_single_pbc<'a, T, C>(
    cutoff: Float,
    //data: &(impl PosProvider + ?Sized),
    data: impl Iterator<Item = &'a Pos>,
    ids: impl Iterator<Item = usize>,
    pbox: &PeriodicBox,
    pbc_dims: PbcDims,
) -> Result<C, DistanceSearchError>
where
    T: DistanceSearchOutput + Send + Sync,
    C: FromIterator<T> + FromParallelIterator<T>,
{
    let grid = periodic_grid(cutoff, data, ids, pbox, pbc_dims)?;

    let plan = search_plan(&grid, None, pbc_dims);

    // Cycle over search plan and perform search for each cell pair
    Ok(collect_cell_pairs(plan, |pair, found| {
        search_cell_pair_single_pbc(cutoff * cutoff, &grid, pair, pbox, found);
    }))
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
            super::distance_search_single_pbc(1., points.iter(), 0..2, &b, PBC_FULL).unwrap();
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
                    )
                    .unwrap();
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
                    )
                    .unwrap();
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
                    )
                    .unwrap();
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
        let left = state(&points[..1]);
        let right = state(&points[1..]);
        let found: Vec<usize> =
            super::distance_search_within_pbc(0.5, &left, &right, 0..1, 1..2, &b, pbc).unwrap();
        assert_eq!(found, vec![0]);
    }
    fn state(points: &[Pos]) -> State {
        State {
            coords: points.to_vec(),
            ..State::default()
        }
    }

    #[test]
    fn invalid_search_inputs_return_errors() {
        use super::*;
        let data = state(&[Pos::origin(), Pos::new(0.25, 0., 0.)]);
        let b = PeriodicBox::from_matrix(Matrix3f::identity()).unwrap();
        for cutoff in [
            0.,
            -0.5,
            Float::NAN,
            Float::INFINITY,
            Float::MAX,
            Float::MIN_POSITIVE,
        ] {
            assert!(matches!(
                distance_search_single::<(usize, usize), Vec<_>>(cutoff, &data, 0..2),
                Err(DistanceSearchError::InvalidCutoff(_))
            ));
            assert!(matches!(
                distance_search_single_pbc::<(usize, usize), Vec<_>>(
                    cutoff,
                    data.iter_pos(),
                    0..2,
                    &b,
                    PBC_FULL
                ),
                Err(DistanceSearchError::InvalidCutoff(_))
            ));
            assert!(
                distance_search_double::<usize, Vec<_>>(cutoff, &data, &data, 0..2, 0..2).is_err()
            );
            assert!(
                distance_search_double_pbc::<usize, Vec<_>>(
                    cutoff,
                    data.iter_pos(),
                    data.iter_pos(),
                    0..2,
                    0..2,
                    &b,
                    PBC_FULL
                )
                .is_err()
            );
            assert!(
                distance_search_within::<Vec<_>>(
                    cutoff,
                    &data,
                    &data,
                    0..2,
                    0..2,
                    &Vector3f::zeros(),
                    &Vector3f::repeat(1.)
                )
                .is_err()
            );
            assert!(
                distance_search_within_pbc::<Vec<_>>(
                    cutoff,
                    &data,
                    &data,
                    0..2,
                    0..2,
                    &b,
                    PBC_FULL
                )
                .is_err()
            );
        }
        for n in [1, 3] {
            assert_eq!(
                distance_search_single::<usize, Vec<_>>(0.5, &data, 0..n),
                Err(DistanceSearchError::LengthMismatch)
            );
            assert_eq!(
                distance_search_single_pbc::<usize, Vec<_>>(
                    0.5,
                    data.iter_pos(),
                    0..n,
                    &b,
                    PBC_FULL
                ),
                Err(DistanceSearchError::LengthMismatch)
            );
            assert_eq!(
                distance_search_double_vdw::<usize, Vec<_>>(
                    &data,
                    &data,
                    &vec![0.25; n],
                    &[0.25; 2]
                ),
                Err(DistanceSearchError::LengthMismatch)
            );
            assert_eq!(
                distance_search_double_vdw_pbc::<usize, Vec<_>>(
                    data.iter_pos(),
                    data.iter_pos(),
                    &[0.25; 2],
                    &vec![0.25; n],
                    &b,
                    PBC_FULL
                ),
                Err(DistanceSearchError::LengthMismatch)
            );
        }
        for value in [-1., Float::NAN, Float::INFINITY] {
            assert!(matches!(
                distance_search_double_vdw::<usize, Vec<_>>(&data, &data, &[value; 2], &[0.25; 2]),
                Err(DistanceSearchError::InvalidRadius(_))
            ));
        }
        for value in [Float::NAN, Float::INFINITY, Float::NEG_INFINITY] {
            let bad = state(&[Pos::new(value, 0., 0.)]);
            assert_eq!(
                distance_search_single::<usize, Vec<_>>(0.5, &bad, 0..1),
                Err(DistanceSearchError::InvalidCoordinates)
            );
            assert_eq!(
                distance_search_single_pbc::<usize, Vec<_>>(
                    0.5,
                    bad.iter_pos(),
                    0..1,
                    &b,
                    PBC_FULL
                ),
                Err(DistanceSearchError::InvalidCoordinates)
            );
        }
        assert_eq!(
            distance_search_single_pbc::<usize, Vec<_>>(
                0.5,
                data.iter_pos(),
                0..2,
                &PeriodicBox::default(),
                PBC_FULL
            ),
            Err(DistanceSearchError::InvalidBox)
        );
    }

    #[test]
    fn empty_searches_are_empty_and_still_check_lengths() {
        use super::*;
        let empty = State::default();
        let data = state(&[Pos::new(1000., 1000., 1000.)]);
        let b = PeriodicBox::from_matrix(Matrix3f::identity()).unwrap();
        assert!(
            distance_search_single::<usize, Vec<_>>(0.5, &empty, 0..0)
                .unwrap()
                .is_empty()
        );
        assert!(
            distance_search_double::<usize, Vec<_>>(0.5, &empty, &data, 0..0, 0..1)
                .unwrap()
                .is_empty()
        );
        assert!(
            distance_search_double_vdw::<usize, Vec<_>>(&empty, &data, &[], &[0.25])
                .unwrap()
                .is_empty()
        );
        assert!(
            distance_search_double_vdw_pbc::<usize, Vec<_>>(
                empty.iter_pos(),
                empty.iter_pos(),
                &[],
                &[],
                &b,
                PBC_FULL
            )
            .unwrap()
            .is_empty()
        );
        assert_eq!(
            distance_search_double_vdw::<usize, Vec<_>>(&data, &data, &[], &[]),
            Err(DistanceSearchError::LengthMismatch)
        );
        assert_eq!(
            distance_search_single::<usize, Vec<_>>(0.5, &empty, 0..1),
            Err(DistanceSearchError::LengthMismatch)
        );
    }

    #[test]
    fn translated_and_sparse_inputs_preserve_pairs() {
        use super::*;
        let points = [Pos::origin(), Pos::new(0.5, 0., 0.), Pos::new(0., 2., 0.)];
        for shift in [0., 1000., -1000.] {
            let data = state(&points.map(|p| p + Vector3f::repeat(shift)));
            let (lo, hi) = compute_bounding_box_single(1., data.iter_pos()).unwrap();
            assert!((hi - lo).max() < 5.);
            let pairs: Vec<(usize, usize)> = distance_search_single(1., &data, 0..3).unwrap();
            assert_eq!(pairs.len(), 1);
            assert!(pairs[0] == (0, 1) || pairs[0] == (1, 0));
        }
        let data = state(&[
            Pos::origin(),
            Pos::new(0.5, 0., 0.),
            Pos::from(Vector3f::repeat(1.0e10)),
        ]);
        let (lo, hi) = compute_bounding_box_single(1., data.iter_pos()).unwrap();
        let grid =
            Grid::from_cutoff_and_min_max(1., &lo, &hi, cell_limit(data.iter_pos())).unwrap();
        assert!(grid.cells.len() <= 24);
        let pairs: Vec<(usize, usize)> = distance_search_single(1., &data, 0..3).unwrap();
        assert_eq!(pairs, vec![(0, 1)]);
        assert!(
            bounded_dims([usize::MAX; 3], MAX_GRID_CELLS)
                .iter()
                .product::<usize>()
                <= MAX_GRID_CELLS
        );
    }

    #[test]
    fn cutoff_and_upper_grid_faces_are_inclusive() {
        use super::*;
        // Absolute EPSILON cannot pad the upper bound at this translation.
        let data = state(&[Pos::new(1024., 0., 0.)]);
        let target = Pos::new(1023., 0., 0.);
        let lower = target.coords.add_scalar(-1. - Float::EPSILON);
        let upper = target.coords.add_scalar(1. + Float::EPSILON);
        let found: Vec<usize> =
            distance_search_within(1., &data, &target, 0..1, 0..1, &lower, &upper).unwrap();
        assert_eq!(found, vec![0]);
        for distance in [1.0 as Float - Float::EPSILON, 1., 1. + Float::EPSILON] {
            let data = state(&[Pos::origin(), Pos::new(distance, 0., 0.)]);
            let found: Vec<(usize, usize)> = distance_search_single(1., &data, 0..2).unwrap();
            assert_eq!(found.len(), usize::from(distance <= 1.));
        }
    }

    // Independent reference: enumerate images directly, without PeriodicBox's
    // minimum-image routine. These bounded coordinates need only nearby images.
    fn image_distance2(a: &Pos, b: &Pos, matrix: &Matrix3f, pbc: PbcDims) -> Float {
        let mut best = Float::INFINITY;
        for x in -3..=3 {
            for y in -3..=3 {
                for z in -3..=3 {
                    let shift = [x, y, z];
                    if (0..3).any(|d| !pbc.get_dim(d) && shift[d] != 0) {
                        continue;
                    }
                    let delta = b - a + matrix * Vector3f::new(x as Float, y as Float, z as Float);
                    best = best.min(delta.norm_squared());
                }
            }
        }
        best
    }

    #[test]
    fn searches_match_independent_images_and_unequal_radii() {
        use super::*;
        use std::collections::BTreeSet;
        let matrix = Matrix3f::new(3., 0.75, 0.25, 0., 4., 0.5, 0., 0., 5.);
        let b = PeriodicBox::from_matrix(matrix).unwrap();
        let mut seed = 991_u64;
        let points: Vec<_> = (0..30)
            .map(|_| {
                Pos::from(
                    matrix
                        * Vector3f::from_fn(|_, _| {
                            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                            ((seed >> 32) as u32 as Float / u32::MAX as Float) * 2. - 0.5
                        }),
                )
            })
            .collect();
        let left = state(&points[..15]);
        let right = state(&points[15..]);
        let radii1: Vec<_> = (0..15).map(|i| 0.1 + (i % 5) as Float * 0.2).collect();
        let radii2: Vec<_> = (0..15).map(|i| 0.05 + (i % 3) as Float * 0.15).collect();
        for mask in 0..8 {
            let pbc = PbcDims::new(mask & 1 != 0, mask & 2 != 0, mask & 4 != 0);
            let mut expected_vdw = BTreeSet::new();
            for i in 0..15 {
                for j in 0..15 {
                    let cutoff = radii1[i] + radii2[j] + Float::EPSILON;
                    if image_distance2(&points[i], &points[j + 15], &matrix, pbc) <= cutoff * cutoff
                    {
                        expected_vdw.insert((i, j));
                    }
                }
            }
            let got: Vec<(usize, usize)> = distance_search_double_vdw_pbc(
                left.iter_pos(),
                right.iter_pos(),
                &radii1,
                &radii2,
                &b,
                pbc,
            )
            .unwrap();
            assert_eq!(got.len(), expected_vdw.len());
            assert_eq!(got.into_iter().collect::<BTreeSet<_>>(), expected_vdw);
            if pbc == PBC_NONE {
                let got: Vec<(usize, usize)> =
                    distance_search_double_vdw(&left, &right, &radii1, &radii2).unwrap();
                assert_eq!(got.into_iter().collect::<BTreeSet<_>>(), expected_vdw);
            }
            for cutoff in [0.5, 1.5, 6.] {
                let mut expected = BTreeSet::new();
                for i in 0..15 {
                    for j in 0..15 {
                        if image_distance2(&points[i], &points[j + 15], &matrix, pbc)
                            <= cutoff * cutoff
                        {
                            expected.insert((i, j));
                        }
                    }
                }
                let got: Vec<(usize, usize)> = distance_search_double_pbc(
                    cutoff,
                    left.iter_pos(),
                    right.iter_pos(),
                    0..15,
                    0..15,
                    &b,
                    pbc,
                )
                .unwrap();
                assert_eq!(got.len(), expected.len());
                assert_eq!(got.into_iter().collect::<BTreeSet<_>>(), expected);
                let expected_within: BTreeSet<_> = expected.iter().map(|&(i, _)| i).collect();
                let got: Vec<usize> =
                    distance_search_within_pbc(cutoff, &left, &right, 0..15, 0..15, &b, pbc)
                        .unwrap();
                assert_eq!(got.len(), expected_within.len());
                assert_eq!(got.into_iter().collect::<BTreeSet<_>>(), expected_within);
                if pbc == PBC_NONE {
                    let got: Vec<(usize, usize)> =
                        distance_search_double(cutoff, &left, &right, 0..15, 0..15).unwrap();
                    assert_eq!(got.into_iter().collect::<BTreeSet<_>>(), expected);
                }
            }
        }
    }
}
