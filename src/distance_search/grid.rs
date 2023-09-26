use crate::core::{PbcDims, PeriodicBox, Pos, Vector3f};
use nalgebra::Vector3;

pub type CellLoc = nalgebra::Vector3<usize>;

// Grid cell with points and indexes
#[derive(Debug, Clone, Default)]
pub struct GridCellData {
    pub ids: Vec<usize>,
    pub coords: Vec<Pos>,
}

impl GridCellData {
    pub fn add(&mut self, id: usize, coord: &Pos) {
        self.ids.push(id);
        self.coords.push(coord.clone());
    }

    pub fn len(&self) -> usize {
        self.ids.len()
    }
}

//================================================================================================

#[derive(Debug, Clone)]
struct CellPair {
    c1: CellLoc,
    c2: CellLoc,
    wrapped: PbcDims,
}

// Iterator over cell pairs in the grid of given size and periodicity
struct CellPairIter {
    grid_size: [usize; 3],
    periodic_dims: PbcDims,
    // Current grid location
    grid_loc: CellLoc,
    // Current buffer of pairs
    buf_iter: <Vec<CellPair> as IntoIterator>::IntoIter,
}

pub static STENCIL: [[Vector3<usize>; 2]; 14] = [
    // Center
    [Vector3::new(0, 0, 0), Vector3::new(0, 0, 0)],
    // Edges
    [Vector3::new(0, 0, 0), Vector3::new(1, 0, 0)], //X
    [Vector3::new(0, 0, 0), Vector3::new(0, 1, 0)], //Y
    [Vector3::new(0, 0, 0), Vector3::new(0, 0, 1)], //Z
    // Face angles
    [Vector3::new(0, 0, 0), Vector3::new(1, 1, 0)], //XY
    [Vector3::new(0, 0, 0), Vector3::new(1, 0, 1)], //XZ
    [Vector3::new(0, 0, 0), Vector3::new(0, 1, 1)], //YZ
    // Far angle
    [Vector3::new(0, 0, 0), Vector3::new(1, 1, 1)], //XYZ
    // Face-diagonals
    [Vector3::new(1, 0, 0), Vector3::new(0, 1, 0)], // XY
    [Vector3::new(1, 0, 0), Vector3::new(0, 0, 1)], // XZ
    [Vector3::new(0, 1, 0), Vector3::new(0, 0, 1)], // YZ
    // Cross-diagonals
    [Vector3::new(1, 1, 0), Vector3::new(0, 0, 1)], // XY-Z
    [Vector3::new(1, 0, 1), Vector3::new(0, 1, 0)], // XZ-Y
    [Vector3::new(0, 1, 1), Vector3::new(1, 0, 0)], // YZ-X
];

impl CellPairIter {
    fn new(grid_size: &[usize; 3], periodic_dims: &PbcDims) -> Self {
        Self {
            grid_size: *grid_size,
            periodic_dims: *periodic_dims,
            buf_iter: Vec::new().into_iter(),
            grid_loc: CellLoc::zeros(),
        }
    }

    fn validate_cell_pair(&self, mut pair: CellPair) -> Option<CellPair> {
        for d in 0..3 {
            let sz = self.grid_size[d];
            let pbc = self.periodic_dims[d];
            match sz {
                // Corner cases:
                1 => {
                    // Always ignore any pairs beyond limits (only 0 allowed)
                    if pair.c1[d] != 0 || pair.c2[d] != 0 {
                        return None;
                    }
                    // For only one posisble valid pair 0:0 set periodicity
                    if pbc {
                        pair.wrapped[d] = true;
                    }
                }
                2 => {
                    // Always ignore any pairs beyond limits
                    if pair.c1[d] >= sz || pair.c2[d] >= sz {
                        return None;
                    }
                    // Set periodicity for 0:1 and 1:0 valid pairs
                    // For 0:0 and 1:1 pairs no periodicity is needed
                    if pbc && pair.c1[d] != pair.c2[d] {
                        pair.wrapped[d] = true;
                    }
                }
                // Usual case
                _ => {
                    if pair.c1[d] >= sz {
                        // point beyond the right edge
                        if pbc {
                            pair.c1[d] = pair.c1[d] % sz; // Wrap this dimension
                            pair.wrapped[d] = true;
                        } else {
                            return None; // don't use this pair
                        }
                    }

                    if pair.c2[d] >= sz {
                        // point beyond the right edge
                        if pbc {
                            pair.c2[d] = pair.c2[d] % sz; // Wrap this dimension
                            pair.wrapped[d] = true;
                        } else {
                            return None; // don't use this pair
                        }
                    }
                }
            };
        }
        Some(pair)
    }

    fn generate_pairs_around_cell(&mut self) {
        let mut buf = Vec::<CellPair>::with_capacity(STENCIL.len());
        // Loop over all stencil offsets
        for [o1, o2] in &STENCIL {
            let pp = CellPair {
                c1: self.grid_loc + o1,
                c2: self.grid_loc + o2,
                wrapped: [false, false, false],
            };
            if let Some(pair) = self.validate_cell_pair(pp) {
                // Pair is valid.
                buf.push(pair);
            }
        }
        // Convert buf to iterator
        self.buf_iter = buf.into_iter();
    }

    fn next_cell(&self) -> Option<CellLoc> {
        if self.grid_loc[0] == self.grid_size[0] - 1
            && self.grid_loc[1] == self.grid_size[1] - 1
            && self.grid_loc[2] == self.grid_size[2] - 1
        {
            return None;
        }
        Some(CellLoc::from_fn(|d, _| {
            (self.grid_loc[d] + 1) % self.grid_size[d]
        }))
    }
}

impl Iterator for CellPairIter {
    type Item = CellPair;

    fn next(&mut self) -> Option<Self::Item> {
        let ret = self.buf_iter.next();
        if let None = ret {
            // End of buffer reached, go to the next cell
            match self.next_cell() {
                Some(c) => self.grid_loc = c,
                None => return None, // The end of grid reached
            }
            // Generate new buffer
            self.generate_pairs_around_cell();
        };
        ret
    }
}

#[derive(Debug, Clone)]
pub struct Grid<T> {
    data: ndarray::Array3<T>,
    pbc: Option<(PeriodicBox, PbcDims)>,
}

pub trait IdPosIterator<'a>: ExactSizeIterator<Item = (usize, &'a Pos)> {}
impl<'a, T> IdPosIterator<'a> for T where T: ExactSizeIterator<Item = (usize, &'a Pos)> {}

impl<T> Grid<T>
where
    T: Default,
{
    pub fn new(sz: [usize; 3]) -> Self {
        Self {
            data: ndarray::Array3::<T>::from_shape_simple_fn(sz, || T::default()),
            pbc: None,
        }
    }

    pub fn dim(&self) -> [usize; 3] {
        let d = self.data.dim();
        [d.0, d.1, d.2]
    }

    pub fn cell_pair_iter(&self) -> CellPairIter {
        CellPairIter::new(
            &self.dim(),
            match self.pbc {
                Some((_, pbc_dims)) => &pbc_dims,
                None => &[false, false, false],
            },
        )
    }
}

impl Grid<GridCellData> {
    pub fn populate<'a>(
        &mut self,
        id_pos: impl IdPosIterator<'a>,
        lower: &Vector3f,
        upper: &Vector3f,
    ) {
        let dim = self.dim();
        let dim_sz = upper - lower;
        'outer: for (id, pos) in id_pos {
            let mut ind = [0usize, 0, 0];
            for d in 0..3 {
                let n = (dim[d] as f32 * (pos[d] - lower[d]) / dim_sz[d]).floor() as isize;
                if n < 0 || n >= dim[d] as isize {
                    continue 'outer;
                } else {
                    ind[d] = n as usize;
                }
            }
            self.data[ind].add(id, pos);
        }
    }

    pub fn populate_periodic<'a>(
        &mut self,
        id_pos: impl IdPosIterator<'a>,
        box_: &PeriodicBox,
        pbc_dims: &PbcDims,
    ) {
        let dim = self.dim();

        'outer: for (id, pos) in id_pos {
            let wrapped = box_.wrap_vector_dims(&pos.coords, &pbc_dims);
            let rel = box_.to_box_coords(&wrapped);
            let mut ind = [0usize, 0, 0];
            for d in 0..3 {
                let n = (dim[d] as f32 * rel[d]).floor() as isize;
                // If dimension in not periodic and
                // out of bounds - skip the point
                if pbc_dims[d] && (n > dim[d] as isize || n < 0) {
                    continue 'outer;
                }
                // Correct for possible minor numeric errors
                ind[d] = n.clamp(0, dim[d] as isize - 1) as usize;
            }
            self.data[ind].add(id, pos);
        }
        self.pbc = Some((box_.clone(),pbc_dims.clone()));
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

#[test]
fn test_grid() {
    use crate::io::*;
    use std::iter::zip;
    let mut r = FileHandler::new_reader("tests/no_ATP.pdb").unwrap();
    let st = r.read_next_state().unwrap().unwrap();

    let mut gr = Grid::new([10, 10, 10]);
    gr.populate_periodic(
        zip(0..st.coords.len(), st.coords.iter()),
        &st.box_,
        &[true, true, true],
    );
}
