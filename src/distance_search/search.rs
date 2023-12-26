use super::grid::*;
use crate::{
    core::{PbcDims, PeriodicBox, Pos, Vector3f, IdPosIterator, MeasurePos},
    distance_search::cell_pair_iterator::CellPairIter,
};
//use rayon::prelude::*;
use std::collections::HashMap;

pub struct FoundPair {
    pub i: usize,
    pub j: usize,
    pub d: f32,
}

#[derive(Debug,Default)]
pub struct SearchConnectivity (HashMap<usize,Vec<usize>>);

impl SearchConnectivity {
    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn iter(&self) -> SearchConnectivityIter {
        SearchConnectivityIter(self.0.iter())
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

pub struct SearchConnectivityIter<'a>(
    std::collections::hash_map::Iter<'a,usize,Vec<usize>>
);

impl<'a> Iterator for SearchConnectivityIter<'a> {
    type Item = (&'a usize,&'a Vec<usize>);
    fn next(&mut self) -> Option<Self::Item> {
        self.0.next()
    }
}

impl IntoIterator for SearchConnectivity {
    type Item = (usize,Vec<usize>);
    type IntoIter = std::collections::hash_map::IntoIter<usize,Vec<usize>>;
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

//========================================================
pub struct DistanceSearcherSingle {
    grid: Grid<GridCellData>,
    cutoff: f32,
}

impl DistanceSearcherSingle {
    pub fn new<'a>(
        cutoff: f32,
        data: impl IdPosIterator<'a>,
        lower: &Vector3f,
        upper: &Vector3f,
    ) -> Self {
        // Create grid
        let mut grid = Grid::<GridCellData>::from_cuoff_and_min_max(cutoff,lower,upper);
        grid.populate(data, lower, upper);
        // Create an instance
        Self { grid, cutoff }
    }

    pub fn new_periodic<'a>(
        cutoff: f32,
        data: impl IdPosIterator<'a>,
        box_: &PeriodicBox,
        periodic_dims: &PbcDims,
    ) -> Self {
        // Create grid
        let mut grid = Grid::<GridCellData>::from_cuoff_and_box(cutoff,box_);
        grid.populate_periodic(
            data,
            box_,
            &periodic_dims,
        );
        // Create an instance
        Self { grid, cutoff }
    }

    fn dist_periodic(&self, p1: &Pos, p2: &Pos) -> f32 {
        let pbc = self.grid.pbc.as_ref().unwrap();
        pbc.box_.distance_squared(p1, p2, &pbc.dims)
    }

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
            Some(_) => Self::dist_periodic,
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
pub struct DistanceSearcherDouble {
    grid1: Grid<GridCellData>,
    grid2: Grid<GridCellData>,
    cutoff: f32,
}

impl DistanceSearcherDouble {
    pub fn new<'a>(
        cutoff: f32,
        data1: impl IdPosIterator<'a>,
        data2: impl IdPosIterator<'a>,
        lower: &Vector3f,
        upper: &Vector3f,
    ) -> Self {
        // Create grids
        let mut grid1 = Grid::<GridCellData>::from_cuoff_and_min_max(cutoff, lower, upper);
        let mut grid2 = Grid::<GridCellData>::new(grid1.dim());

        grid1.populate(data1, lower, upper);
        grid2.populate(data2, lower, upper);
        // Create an instance
        Self {
            grid1,
            grid2,
            cutoff,
        }
    }

    pub fn new_periodic<'a>(
        cutoff: f32,
        data1: impl IdPosIterator<'a>,
        data2: impl IdPosIterator<'a>,
        box_: &PeriodicBox,
        periodic_dims: &PbcDims,
    ) -> Self {
        // Create grids
        let mut grid1 = Grid::<GridCellData>::from_cuoff_and_box(cutoff, box_);
        let mut grid2 = Grid::<GridCellData>::new(grid1.dim());

        grid1.populate_periodic(data1,box_,&periodic_dims);
        grid2.populate_periodic(data2,box_,&periodic_dims);
        // Create an instance
        Self {
            grid1,
            grid2,
            cutoff,
        }
    }

    fn dist_periodic(&self, p1: &Pos, p2: &Pos) -> f32 {
        let pbc = self.grid1.pbc.as_ref().unwrap(); // First grid as reference
        pbc.box_.distance_squared(p1, p2, &pbc.dims)
    }

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
    fn from_search_results(i: usize, _j: usize, _d: f32) -> Self {
        i
    }
}

impl SearchOutputType for FoundPair {
    fn from_search_results(i: usize, j: usize, d: f32) -> Self {
        Self { i, j, d }
    }
}

impl SearchOutputType for (usize, usize) {
    fn from_search_results(i: usize, j: usize, _d: f32) -> Self {
        (i, j)
    }
}

impl SearchOutputType for (usize, usize, f32) {
    fn from_search_results(i: usize, j: usize, d: f32) -> Self {
        (i, j, d)
    }
}

//==================================================================
// Tests

#[test]
fn test_single_periodic() {
    use crate::io::*;
    let mut r = FileHandler::new_reader("tests/no_ATP.pdb").unwrap();
    let st = r.read_next_state().unwrap().unwrap();

    let searcher = DistanceSearcherSingle::new_periodic(
        0.3,
        (0..st.coords.len()).map(|i| (i,&st.coords[i])),
        st.box_.as_ref().unwrap(),
        &[true, true, true],
    );
    let found: SearchConnectivity = searcher.search();
    println!("{:?}", found.len())
}

#[test]
fn test_single_non_periodic() {
    use crate::io::*;
    let mut r = FileHandler::new_reader("tests/no_ATP.pdb").unwrap();
    let st = r.read_next_state().unwrap().unwrap();

    let (lower,upper) = st.min_max();
    let searcher = DistanceSearcherSingle::new(
        0.3,
        (0..st.coords.len()).map(|i| (i,&st.coords[i])),
        &lower.coords,
        &upper.coords,
    );
    let found: Vec<(usize,usize)> = searcher.search();
    println!("{:?}", found.len())
}

#[test]
fn test_double_periodic() {
    use crate::io::*;
    let mut r = FileHandler::new_reader("tests/no_ATP.pdb").unwrap();
    let st = r.read_next_state().unwrap().unwrap();

    let searcher = DistanceSearcherDouble::new_periodic(
        0.3,
        (0..st.coords.len()).map(|i| (i,&st.coords[i])),
        (0..st.coords.len()).map(|i| (i,&st.coords[i])),
        &st.box_.as_ref().unwrap(),
        &[true, true, true],
    );
    let found: Vec<usize> = searcher.search();
    println!("{:?}", found.len())
}
