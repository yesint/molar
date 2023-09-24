use super::grid::*;
use crate::core::{IndexIterator, PbcDims, PeriodicBox, Vector3f};
use nalgebra::Vector3;
use num_traits::clamp_min;

struct CellPair {
    c1: CellLoc,
    c2: CellLoc,
    wrapped: PbcDims,
}

trait GridSearchable {
    fn process_cell_pair(&self, pair: CellPair);
}

struct SearcherSingleGrid<'a> {
    grid: Grid,
    driver: SearchDriver<'a, Self>,
}

struct SearcherDoubleGrid<'a> {
    grid1: Grid,
    grid2: Grid,
    driver: SearchDriver<'a, Self>,
}

impl<'a> GridSearchable for SearcherSingleGrid<'a> {
    fn process_cell_pair(&self, pair: CellPair) {
        // Search within a single grid
        if pair.c1 == pair.c2 {
            // Inside cell
            Self::search_inside_cell(&self.grid[&pair.c1]);
        } else {
            // Between cells
            Self::search_between_cells(&self.grid[&pair.c1], &self.grid[&pair.c2]);
        }
    }
}

impl<'a> SearcherSingleGrid<'a> {
    fn search_inside_cell(cell: &GridCell){
        todo!()
    }

    fn search_between_cells(cell1: &GridCell, cell2: &GridCell){
        todo!()
    }
}

impl<'a> GridSearchable for SearcherDoubleGrid<'a> {
    fn process_cell_pair(&self, pair: CellPair) {
        // Search in two grids
        // First search c1->grid1, c2->grid2
        Self::search_between_cells(&self.grid1[&pair.c1], &self.grid2[&pair.c2]);
        if pair.c1 != pair.c2 {
            // Then search c2->grid1, c1->grid2
            Self::search_between_cells(&self.grid1[&pair.c2], &self.grid2[&pair.c1]);
        }
    }
}

impl<'a> SearcherDoubleGrid<'a> {
    fn search_between_cells(cell1: &GridCell, cell2: &GridCell){
        todo!()
    }
}

struct SearchDriver<'a, T>
where
    T: GridSearchable,
{
    parent_searcher: &'a T,
    cutoff: f32,
    grid_size: [usize;3],
    periodic_dims: PbcDims,
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

impl<'a, T> SearchDriver<'a, T>
where
    T: GridSearchable + 'a,
{
    fn grid_size_from_cutoff_and_extents(cutoff: f32, extents: &Vector3f) -> [usize;3] {
        let mut res = [0,0,0];
        // Cell size should be >= cutoff for all dimentions
        for d in 0..3 {
            res[d] = clamp_min((extents[d] / cutoff).floor() as usize, 1);
        }
        res
    }

    fn grid_size(cutoff: f32, min: &Vector3f, max: &Vector3f) -> [usize;3] {
        Self::grid_size_from_cutoff_and_extents(cutoff,&(max - min))
    }

    // Periodic variant
    fn grid_size_periodic(cutoff: f32, box_: &PeriodicBox) -> [usize;3] {
        Self::grid_size_from_cutoff_and_extents(cutoff,&box_.get_extents())
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

    fn index_to_pos(&self, ind: usize) -> Vector3<usize> {
        // i = z+Nz*y+Nz*Ny*x
        Vector3::new(
            (ind / self.grid_size[2] / self.grid_size[1]) % self.grid_size[0],
            (ind / self.grid_size[2]) % self.grid_size[1],
            ind % self.grid_size[2],
        )
    }

    // Search around given cell
    // Each cell in the grid should be searched just once
    fn search_around_cell(&self, cell: &CellLoc) {
        // Loop over all stencil offsets
        for [o1, o2] in STENCIL {
            let pp = CellPair {
                c1: cell + o1,
                c2: cell + o2,
                wrapped: [false, false, false],
            };
            if let Some(pair) = self.validate_cell_pair(pp) {
                // Pair is valid. Call the method of parent searcher to process it
                self.parent_searcher.process_cell_pair(pair);
            }
        }
    }
}
