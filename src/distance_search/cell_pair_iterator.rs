use crate::core::{PbcDims, PBC_NONE};
use super::grid::{CellLoc, CellPair};
use nalgebra::Vector3;

// Search stencil
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

// Iterator over cell pairs in the grid of given size and periodicity
pub struct CellPairIter {
    grid_size: [usize; 3],
    periodic_dims: PbcDims,
    // Current grid location
    grid_loc: CellLoc,
    // Iterator over buffer of pairs for current central cell
    //buf_iter: <Vec<CellPair> as IntoIterator>::IntoIter,
    buf_iter: std::vec::IntoIter<CellPair>,
}

impl CellPairIter {
    pub fn new(grid_size: &[usize; 3], periodic_dims: &PbcDims) -> Self {
        let mut ret = Self {
            grid_size: *grid_size,
            periodic_dims: *periodic_dims,
            buf_iter: Vec::new().into_iter(),
            grid_loc: CellLoc::zeros(),
        };
        ret.gen_buf_iter();
        ret
    }

    fn cell_pair_from_stencil(&self, mut c1: CellLoc, mut c2: CellLoc) -> Option<CellPair> {
        let mut wrapped = PBC_NONE;
        for d in 0..3 {
            let sz = self.grid_size[d];
            let pbc = self.periodic_dims[d];
            match sz {
                // Corner cases:
                1 => {
                    // Always ignore any pairs beyond limits (only 0 allowed)
                    if c1[d] != 0 || c2[d] != 0 {
                        return None;
                    }
                    // For only one posisble valid pair 0:0 set periodicity
                    if pbc {
                        wrapped[d] = true;
                    }
                }
                2 => {
                    // Always ignore any pairs beyond limits
                    if c1[d] >= sz || c2[d] >= sz {
                        return None;
                    }
                    // Set periodicity for 0:1 and 1:0 valid pairs
                    // For 0:0 and 1:1 pairs no periodicity is needed
                    if pbc && c1[d] != c2[d] {
                        wrapped[d] = true;
                    }
                }
                // Usual case
                _ => {
                    if c1[d] >= sz {
                        // point beyond the right edge
                        if pbc {
                            c1[d] = c1[d] % sz; // Wrap this dimension
                            wrapped[d] = true;
                        } else {
                            return None; // don't use this pair
                        }
                    }

                    if c2[d] >= sz {
                        // point beyond the right edge
                        if pbc {
                            c2[d] = c2[d] % sz; // Wrap this dimension
                            wrapped[d] = true;
                        } else {
                            return None; // don't use this pair
                        }
                    }
                }
            };
        }
        Some(CellPair { c1, c2, wrapped })
    }

    fn gen_buf_iter(&mut self) {
        let mut buf = Vec::<CellPair>::with_capacity(STENCIL.len());
        // Loop over all stencil offsets
        for [o1, o2] in &STENCIL {
            if let Some(pair) = self.cell_pair_from_stencil(self.grid_loc + o1, self.grid_loc + o2) {
                // Pair is valid.
                buf.push(pair);
            }
        }
        // Convert buf to iterator and return it
        self.buf_iter = buf.into_iter();
    }

    fn next_cell(&self) -> Option<CellLoc> {
        let mut loc = self.grid_loc;
        loc[2] += 1;
        if loc[2] == self.grid_size[2] {
            loc[2] = 0;
            loc[1] += 1;
            if loc[1] == self.grid_size[1] {
                loc[1] = 0;
                loc[0] += 1;
                if loc[0] == self.grid_size[0] {
                    return None;
                }
            }
        }
        Some(loc)
    }
}

impl Iterator for CellPairIter {
    type Item = CellPair;

    fn next(&mut self) -> Option<Self::Item> {
        let mut ret = self.buf_iter.next();
        while let None = ret {
            // End of buffer reached, go to the next cell
            match self.next_cell() {
                Some(c) => self.grid_loc = c,
                None => return None, // The end of grid reached
            }
            // Generate new buffer and return its iterator
            self.gen_buf_iter();
            ret = self.buf_iter.next();
        }
        ret
    }
}
