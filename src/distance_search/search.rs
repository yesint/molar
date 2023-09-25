use super::grid::*;
use crate::core::{IndexIterator, PbcDims, PeriodicBox, State, Vector3f, Pos};
use nalgebra::Vector3;
use num_traits::clamp_min;

struct CellPair {
    c1: CellLoc,
    c2: CellLoc,
    wrapped: PbcDims,
}

trait GridSearchable {
    fn process_cell_pair(&mut self, pair: CellPair);
}

struct SearcherSingleGrid {
    grid: Grid<GridCellData>,
    found: Vec<[usize;2]>,
    cutoff: f32,
    periodic_dims: PbcDims,
    box_: PeriodicBox,    
}

struct SearcherDoubleGrid {
    grid1: Grid<GridCellData>,
    grid2: Grid<GridCellData>,
    //found: Grid<Vec<[usize;2]>>,
}

impl SearcherSingleGrid {
    fn from_state_subset_periodic(cutoff: f32, state: &State, subset: impl IndexIterator) -> Self {
        let natoms = subset.len();
        // Get grid dimensions
        let grid_sz = grid_size_periodic(cutoff, &state.box_);
        // Create an instance
        let mut instance = Self {
            grid: Grid::<GridCellData>::new(grid_sz),
            found: Vec::<[usize;2]>::new(),
            cutoff,
            periodic_dims: [true, true, true],
            box_: state.box_.clone(),
        };
        
        // Fill the data grid
        instance.grid.populate_periodic(
            subset.map(|i| (i, &state.coords[i])),
            &state.box_,
            &instance.periodic_dims,
        );

        instance
    }

    fn search(&mut self) {
        // Create a driver
        let mut driver = SearchDriver {
            cutoff: self.cutoff,
            grid_size: self.grid.dim(),
            periodic_dims: self.periodic_dims,
            parent_searcher: self
        };
        // Perform a search
        driver.search();
    }

    fn search_inside_cell(&mut self,pair: CellPair) {
        let cell = &self.grid[&pair.c1];
        let n = cell.len();
        let cutoff2 = self.cutoff * self.cutoff;

        if pair.wrapped[0] || pair.wrapped[1] || pair.wrapped[2] {
            // Periodic variant
            for i in 0..n-1 {
                for j in i+1..n {
                    let dist = self.box_.distance_squared(&cell.coords[i], &cell.coords[j], &self.periodic_dims);
                    if dist <= cutoff2 {
                        self.found.push([cell.ids[i],cell.ids[j]]);
                    }
                }
            }
        } else {
            // Non-periodic variant
            for i in 0..n-1 {
                for j in i+1..n {
                    let dist = (cell.coords[i]-cell.coords[j]).norm_squared();
                    if dist <= cutoff2 {
                        self.found.push([cell.ids[i],cell.ids[j]]);
                    }
                }
            }
        }
    }

    fn search_between_cells(&mut self,pair: CellPair) {
        let cell1 = &self.grid[&pair.c1];
        let cell2 = &self.grid[&pair.c2];
        let n1 = cell1.len();
        let n2 = cell2.len();
        let cutoff2 = self.cutoff * self.cutoff;

        if pair.wrapped[0] || pair.wrapped[1] || pair.wrapped[2] {
            // Periodic variant
            for i in 0..n1 {
                for j in 0..n2 {
                    let dist = self.box_.distance_squared(&cell1.coords[i], &cell2.coords[j], &self.periodic_dims);
                    if dist <= cutoff2 {
                        self.found.push([cell1.ids[i],cell2.ids[j]]);
                    }
                }
            }
        } else {
            // Non-periodic variant
            for i in 0..n1 {
                for j in 0..n2 {
                    let dist = (cell1.coords[i]-cell2.coords[j]).norm_squared();
                    if dist <= cutoff2 {
                        self.found.push([cell1.ids[i],cell2.ids[j]]);
                    }
                }
            }
        }
    }
}

impl GridSearchable for SearcherSingleGrid {
    fn process_cell_pair(&mut self, pair: CellPair) {
        if self.grid[&pair.c1].len() == 0 || self.grid[&pair.c1].len() == 0 {
            // Nothing to do
            return
        }

        // Search within a single grid
        if pair.c1 == pair.c2 {
            // Inside cell
            self.search_inside_cell(pair);
        } else {
            // Between cells
            self.search_between_cells(pair);
        }
    }
}


impl GridSearchable for SearcherDoubleGrid {
    fn process_cell_pair(&mut self, pair: CellPair) {
        // Search in two grids
        // First search c1->grid1, c2->grid2
        Self::search_between_cells(&self.grid1[&pair.c1], &self.grid2[&pair.c2]);
        if pair.c1 != pair.c2 {
            // Then search c2->grid1, c1->grid2
            Self::search_between_cells(&self.grid1[&pair.c2], &self.grid2[&pair.c1]);
        }
    }
}

impl SearcherDoubleGrid {
    fn search_between_cells(cell1: &GridCellData, cell2: &GridCellData) {
        todo!()
    }
}

struct SearchDriver<'a, T>
where
    T: GridSearchable,
{
    parent_searcher: &'a mut T,
    cutoff: f32,
    grid_size: [usize; 3],
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

fn grid_size_from_cutoff_and_extents(cutoff: f32, extents: &Vector3f) -> [usize; 3] {
    let mut res = [0, 0, 0];
    // Cell size should be >= cutoff for all dimentions
    for d in 0..3 {
        res[d] = clamp_min((extents[d] / cutoff).floor() as usize, 1);
    }
    res
}

pub fn grid_size(cutoff: f32, min: &Vector3f, max: &Vector3f) -> [usize; 3] {
    grid_size_from_cutoff_and_extents(cutoff, &(max - min))
}

// Periodic variant
pub fn grid_size_periodic(cutoff: f32, box_: &PeriodicBox) -> [usize; 3] {
    grid_size_from_cutoff_and_extents(cutoff, &box_.get_extents())
}

impl<'a, T> SearchDriver<'a, T>
where
    T: GridSearchable + 'a,
{
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

    // Main search function
    fn search(&mut self) {
        // Cycle over all grid cells and search each one
        for i in 0..self.grid_size[0] {
            for j in 0..self.grid_size[1] {
                for k in 0..self.grid_size[2] {
                    println!("{i} {j} {k}");
                    self.search_around_cell(&Vector3::new(i,j,k))
                }
            }    
        }
    }

    // Search around given cell
    // Each cell in the grid should be searched just once
    fn search_around_cell(&mut self, cell: &CellLoc) {
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


#[test]
fn test_single() {
    use crate::io::*;
    use std::iter::zip;
    let mut r = FileHandler::new_reader("tests/no_ATP.pdb").unwrap();
    
    let st = r.read_next_state().unwrap().unwrap();
    //let _ = r.read_structure().unwrap();
    

    let mut searcher = SearcherSingleGrid::from_state_subset_periodic(0.3, &st, 0..100);
    searcher.search();
}
