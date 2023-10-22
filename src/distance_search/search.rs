use super::grid::*;
use crate::{
    core::{IndexIterator, ParticleIterator, PbcDims, PeriodicBox, Pos, State, Vector3f},
    distance_search::cell_pair_iterator::CellPairIter,
};
use nalgebra::Vector3;
use num_traits::clamp_min;
use rayon::prelude::*;
use std::{collections::HashMap, borrow::Cow};

pub struct ValidPair {
    pub i: usize,
    pub j: usize,
    pub d: f32,
}

#[derive(Debug,Default)]
pub struct SearchConnectivity (HashMap<usize,Vec<usize>>);

impl SearchConnectivity {
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl FromIterator<(usize,usize)> for SearchConnectivity {
    fn from_iter<T: IntoIterator<Item = (usize,usize)>>(iter: T) -> Self {
        let mut res = Self::default();
        for (i,j) in iter {
            if res.0.contains_key(&i) {
                res.0.get_mut(&i).unwrap().push(j);
            } else {
                res.0.insert(i,vec![j]);
            }
            
            if res.0.contains_key(&j) {
                res.0.get_mut(&j).unwrap().push(i);
            } else {
                res.0.insert(j,vec![i]);
            }
        }
        res
    }
}


pub struct SearcherSingleGrid {
    grid: Grid<GridCellData>,
    cutoff: f32,
}

impl SearcherSingleGrid {
    pub fn from_state_subset(
        cutoff: f32,
        state: &State,
        subset: impl IndexIterator,
        lower: &Vector3f,
        upper: &Vector3f,
    ) -> Self {
        // Get grid dimensions
        let grid_sz = grid_size(cutoff, lower, upper);
        let mut grid = Grid::<GridCellData>::new(grid_sz);

        grid.populate(subset.map(|i| (i, &state.coords[i])), lower, upper);
        // Create an instance
        Self { grid, cutoff }
    }

    pub fn from_state_subset_periodic(
        cutoff: f32,
        state: &State,
        subset: impl IndexIterator,
        periodic_dims: &PbcDims,
    ) -> Self {
        // Get grid dimensions
        let grid_sz = grid_size_periodic(cutoff, &state.box_.as_ref().unwrap());
        let mut grid = Grid::<GridCellData>::new(grid_sz);

        grid.populate_periodic(
            subset.map(|i| (i, &state.coords[i])),
            &state.box_.as_ref().unwrap(),
            &periodic_dims,
        );
        // Create an instance
        Self { grid, cutoff }
    }

    pub fn from_particles_periodic<'a>(
        cutoff: f32,
        particles: impl ParticleIterator<'a>,
        box_: &PeriodicBox,
        periodic_dims: &PbcDims,
    ) -> Self {
        // Get grid dimensions
        let grid_sz = grid_size_periodic(cutoff, box_);
        let mut grid = Grid::<GridCellData>::new(grid_sz);

        grid.populate_periodic(
            particles.map(|p| (p.id, p.pos)),
            box_,
            &periodic_dims,
        );
        // Create an instance
        Self { grid, cutoff }
    }

    #[inline]
    fn dist_periodic(&self, p1: &Pos, p2: &Pos) -> f32 {
        let pbc = self.grid.pbc.as_ref().unwrap();
        pbc.box_.distance_squared(p1, p2, &pbc.dims)
    }

    #[inline]
    fn dist_non_periodic(&self, p1: &Pos, p2: &Pos) -> f32 {
        (p1 - p2).norm_squared()
    }

    // Main search function
    pub fn search<T: SearchOutputType, C: FromIterator<T>>(&self) -> C {
        //pub fn search(&self) -> Vec<ValidPair> {
        // Get periodic dimensions for iterator
        let periodic_dims = match self.grid.pbc.as_ref() {
            Some(pbc) => pbc.dims,
            None => [false, false, false],
        };

        // Get periodic or non-periodic distance function
        let dist_func = match self.grid.pbc.as_ref() {
            Some(pbc) => Self::dist_periodic,
            None => Self::dist_non_periodic,
        };

        // Get iterator over cell pairs
        let iter = CellPairIter::new(&self.grid.dim(), &periodic_dims);

        iter.map(|pair| self.search_cell_pair(pair, dist_func))
            .flatten()
            .collect::<C>()
    }

    fn search_cell_pair<T: SearchOutputType>(
        &self,
        pair: CellPair,
        dist_func: fn(&Self, &Pos, &Pos) -> f32,
    ) -> Vec<T> {
        let n1 = self.grid[&pair.c1].len();
        let n2 = self.grid[&pair.c2].len();

        let mut found = Vec::<T>::new();

        // Nothing to do if cell is empty
        if n1 * n2 == 0 {
            return found;
        }

        let cutoff2 = self.cutoff.powi(2);

        if pair.c1 == pair.c2 {
            // Same cell
            for i in 0..n1 - 1 {
                let p1 = &self.grid[&pair.c1].coords[i];
                for j in i + 1..n1 {
                    let p2 = &self.grid[&pair.c1].coords[j];
                    let dist = dist_func(self, p1, p2);

                    if dist <= cutoff2 {
                        found.push(T::from_search_results(
                            self.grid[&pair.c1].ids[i],
                            self.grid[&pair.c1].ids[j],
                            dist.sqrt(),
                        ));
                    }
                }
            }
        } else {
            // Different cells
            for i in 0..n1 {
                let p1 = &self.grid[&pair.c1].coords[i];
                for j in 0..n2 {
                    let p2 = &self.grid[&pair.c2].coords[j];
                    let dist = dist_func(self, p1, p2);
                    if dist <= cutoff2 {
                        found.push(T::from_search_results(
                            self.grid[&pair.c1].ids[i],
                            self.grid[&pair.c2].ids[j],
                            dist.sqrt(),
                        ));
                    }
                }
            }
        }
        found
    }
}

//=============================================================================
pub struct SearcherDoubleGrid {
    grid1: Grid<GridCellData>,
    grid2: Grid<GridCellData>,
    cutoff: f32,
}

impl SearcherDoubleGrid {
    pub fn from_state_subset(
        cutoff: f32,
        state1: &State,
        subset1: impl IndexIterator,
        state2: &State,
        subset2: impl IndexIterator,
        lower: &Vector3f,
        upper: &Vector3f,
    ) -> Self {
        // Get grid dimensions
        let grid_sz = grid_size(cutoff, lower, upper);
        let mut grid1 = Grid::<GridCellData>::new(grid_sz);
        let mut grid2 = Grid::<GridCellData>::new(grid_sz);

        grid1.populate(subset1.map(|i| (i, &state1.coords[i])), lower, upper);
        grid2.populate(subset2.map(|i| (i, &state2.coords[i])), lower, upper);
        // Create an instance
        Self {
            grid1,
            grid2,
            cutoff,
        }
    }

    pub fn from_state_subset_periodic(
        cutoff: f32,
        state1: &State,
        subset1: impl IndexIterator,
        state2: &State,
        subset2: impl IndexIterator,
        periodic_dims: &PbcDims,
    ) -> Self {
        // Get grid dimensions
        let grid_sz = grid_size_periodic(cutoff, &state1.box_.as_ref().unwrap());
        let mut grid1 = Grid::<GridCellData>::new(grid_sz);
        let mut grid2 = Grid::<GridCellData>::new(grid_sz);

        grid1.populate_periodic(
            subset1.map(|i| (i, &state1.coords[i])),
            &state1.box_.as_ref().unwrap(),
            &periodic_dims,
        );

        grid2.populate_periodic(
            subset2.map(|i| (i, &state2.coords[i])),
            &state1.box_.as_ref().unwrap(), // The same box as the first grid!
            &periodic_dims,
        );
        // Create an instance
        Self {
            grid1,
            grid2,
            cutoff,
        }
    }

    #[inline]
    fn dist_periodic(&self, p1: &Pos, p2: &Pos) -> f32 {
        let pbc = self.grid1.pbc.as_ref().unwrap(); // First grid as reference
        pbc.box_.distance_squared(p1, p2, &pbc.dims)
    }

    #[inline]
    fn dist_non_periodic(&self, p1: &Pos, p2: &Pos) -> f32 {
        (p1 - p2).norm_squared()
    }

    // Main search function
    pub fn search<T: SearchOutputType, C: FromIterator<T>>(&self) -> C {
        // Get periodic dimensions for iterator
        let periodic_dims = match self.grid1.pbc.as_ref() {
            Some(pbc) => pbc.dims,
            None => [false, false, false],
        };

        // Get periodic or non-periodic distance function
        let dist_func = match self.grid1.pbc.as_ref() {
            Some(_) => Self::dist_periodic,
            None => Self::dist_non_periodic,
        };

        // Get iterator over cell pairs
        let iter = CellPairIter::new(&self.grid1.dim(), &periodic_dims);

        // We first search for pair c1->grid1; c2->grid2
        // then for the swapped pair c2->grid1; c1->grid2
        // grid1 always goes first
        iter.map(|pair| {
            let swapped_pair = pair.swaped();
            let mut v1 = self.search_cell_pair(pair, dist_func);
            let v2 = self.search_cell_pair(swapped_pair, dist_func);
            v1.extend(v2);
            v1
        })
        .flatten()
        .collect::<C>()
    }

    fn search_cell_pair<T: SearchOutputType>(
        &self,
        pair: CellPair,
        dist_func: fn(&Self, &Pos, &Pos) -> f32,
    ) -> Vec<T> {
        let n1 = self.grid1[&pair.c1].len();
        let n2 = self.grid2[&pair.c2].len();

        let mut found = Vec::<T>::new();

        // Nothing to do if cell is empty
        if n1 * n2 == 0 {
            return found;
        }

        let cutoff2 = self.cutoff.powi(2);

        for i in 0..n1 {
            let p1 = &self.grid1[&pair.c1].coords[i];
            for j in 0..n2 {
                let p2 = &self.grid2[&pair.c2].coords[j];
                let dist = dist_func(self, p1, p2);
                if dist <= cutoff2 {
                    found.push(T::from_search_results(
                        self.grid1[&pair.c1].ids[i],
                        self.grid2[&pair.c2].ids[j],
                        dist,
                    ));
                }
            }
        }

        found
    }
}

pub trait SearchOutputType {
    fn from_search_results(i: usize, j: usize, d: f32) -> Self;
}

impl SearchOutputType for usize {
    fn from_search_results(i: usize, j: usize, d: f32) -> Self {
        i
    }
}

impl SearchOutputType for ValidPair {
    fn from_search_results(i: usize, j: usize, d: f32) -> Self {
        Self { i, j, d }
    }
}

impl SearchOutputType for (usize, usize) {
    fn from_search_results(i: usize, j: usize, d: f32) -> Self {
        (i, j)
    }
}

//==================================================================

fn grid_size_from_cutoff_and_extents(cutoff: f32, extents: &Vector3f) -> [usize; 3] {
    let mut res = [0, 0, 0];
    // Cell size should be >= cutoff for all dimentions
    for d in 0..3 {
        res[d] = clamp_min((extents[d] / cutoff).floor() as usize, 1);
    }
    res
}

fn grid_size(cutoff: f32, min: &Vector3f, max: &Vector3f) -> [usize; 3] {
    grid_size_from_cutoff_and_extents(cutoff, &(max - min))
}

// Periodic variant
fn grid_size_periodic(cutoff: f32, box_: &PeriodicBox) -> [usize; 3] {
    grid_size_from_cutoff_and_extents(cutoff, &box_.get_extents())
}

//============================================================================
// Tests

#[test]
fn test_single_periodic() {
    use crate::io::*;
    use std::iter::zip;
    let mut r = FileHandler::new_reader("tests/no_ATP.pdb").unwrap();
    let st = r.read_next_state().unwrap().unwrap();
    //let _ = r.read_structure().unwrap();

    let mut searcher = SearcherSingleGrid::from_state_subset_periodic(
        0.3,
        &st,
        0..st.coords.len(),
        &[true, true, true],
    );
    let found: SearchConnectivity = searcher.search();
    println!("{:?}", found.len())
}

#[test]
fn test_single_non_periodic() {
    use crate::io::*;
    use std::iter::zip;
    let mut r = FileHandler::new_reader("tests/no_ATP.pdb").unwrap();
    let st = r.read_next_state().unwrap().unwrap();
    //let _ = r.read_structure().unwrap();

    let mut searcher = SearcherSingleGrid::from_state_subset(
        0.3,
        &st,
        0..st.coords.len(),
        &Vector3f::new(0.0, 0.0, 0.0),
        &Vector3f::new(1.0, 1.0, 1.0),
    );
    let found: SearchConnectivity = searcher.search();
    println!("{:?}", found.len())
}

#[test]
fn test_double_periodic() {
    use crate::io::*;
    use std::iter::zip;
    let mut r = FileHandler::new_reader("tests/no_ATP.pdb").unwrap();
    let st = r.read_next_state().unwrap().unwrap();

    let mut searcher = SearcherDoubleGrid::from_state_subset_periodic(
        0.3,
        &st,
        0..st.coords.len(),
        &st,
        0..st.coords.len(),
        &[true, true, true],
    );
    let found: Vec<usize> = searcher.search();
    println!("{:?}", found.len())
}
