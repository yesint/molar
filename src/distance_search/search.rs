use super::grid::*;
use crate::core::{IndexIterator, PbcDims, PeriodicBox, Pos, State, Vector3f};
use nalgebra::Vector3;
use num_traits::clamp_min;
use rayon::prelude::*;

struct SearcherSingleGrid {
    grid: Grid<GridCellData>,
    found: Vec<[usize; 2]>,
    cutoff: f32,
}

struct SearcherDoubleGrid {
    grid1: Grid<GridCellData>,
    grid2: Grid<GridCellData>,
    found: Vec<[usize; 2]>,
    cutoff: f32,
}

impl SearcherSingleGrid {
    fn from_state_subset_periodic(
        cutoff: f32,
        state: &State,
        subset: impl IndexIterator,
        periodic_dims: &PbcDims,
    ) -> Self {
        let natoms = subset.len();
        // Get grid dimensions
        let grid_sz = grid_size_periodic(cutoff, &state.box_);
        let mut grid = Grid::<GridCellData>::new(grid_sz);
        grid.populate_periodic(
            subset.map(|i| (i, &state.coords[i])),
            &state.box_,
            &periodic_dims,
        );
        // Create an instance
        Self {
            grid,
            found: Vec::<[usize; 2]>::new(),
            cutoff,
        }
    }

    fn search(&mut self) {
        let grid_pbc = self.grid.pbc.clone();
        println!("{:?}",self.grid.dim());
        
        //self.grid.cell_pair_iter().par_bridge().map(|el| println!("{:?}",el)).collect::<()>();

        let mut it = self.grid.cell_pair_iter().par_bridge().map(|pair|{
                //println!("{:?}",pair);
                match grid_pbc.as_ref() {
                    Some(pbc) => {
                        self.search_pair(pair,|p1:&Pos,p2:&Pos| pbc.box_.distance_squared(p1,p2,&pbc.dims));
                    },
                    None => self.search_pair(pair,|p1:&Pos,p2:&Pos| (p1-p2).norm_squared()),
                };
            }
        );
        it.for_each(drop);
    }

    fn search_pair<F>(&self, pair: CellPair, dist_func: F)
    where F: Fn(&Pos,&Pos)->f32
    {
        let n1 = self.grid[&pair.c1].len();
        let n2 = self.grid[&pair.c2].len();

        if n1*n2 == 0 {return}
        println!("{:?}",pair);
        
        let cutoff2 = self.cutoff.powi(2);

        if pair.c1==pair.c2 { // Same cell
            for i in 0..n1 - 1 {
                let p1 = &self.grid[&pair.c1].coords[i];
                for j in i + 1..n1 {
                    let p2 = &self.grid[&pair.c1].coords[j];
                    let dist = dist_func(p1,p2);
                    
                    if dist <= cutoff2 {
                        //self.found.push([self.grid[&pair.c1].ids[i], self.grid[&pair.c1].ids[j]]);
                    }
                }
            }
        } else { // Different cells
            for i in 0..n1 {
                let p1 = &self.grid[&pair.c1].coords[i];
                for j in 0..n2 {
                    let p2 = &self.grid[&pair.c2].coords[j];
                    let dist = dist_func(p1,p2);
                    if dist <= cutoff2 {
                        //self.found.push([self.grid[&pair.c1].ids[i], self.grid[&pair.c2].ids[j]]);
                    }
                }
            }
        }
    }
}
    

impl SearcherDoubleGrid {
    fn search_between_cells(cell1: &GridCellData, cell2: &GridCellData) {
        todo!()
    }
}

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


#[test]
fn test_single() {
    use crate::io::*;
    use std::iter::zip;
    let mut r = FileHandler::new_reader("tests/no_ATP.pdb").unwrap();

    let st = r.read_next_state().unwrap().unwrap();
    //let _ = r.read_structure().unwrap();

    let mut searcher = SearcherSingleGrid::from_state_subset_periodic(0.3, &st, 0..st.coords.len(),&[true,true,true]);
    searcher.search();
    println!("{:?}",searcher.found.len())
}
