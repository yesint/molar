use super::grid::*;
use crate::{
    core::{PbcDims, PeriodicBox, Pos, Vector3f, PBC_NONE},
    distance_search::cell_pair_iterator::CellPairIter,
};
use rayon::prelude::*;
//use std::collections::HashMap;

pub static SERIAL_LIMIT: usize = 2;
const INIT_BUF_SIZE: usize = 1000;

pub struct FoundPair {
    pub i: usize,
    pub j: usize,
    pub d: f32,
}

#[derive(Debug, Default)]
pub struct SearchConnectivity(rustc_hash::FxHashMap<usize, Vec<usize>>);

impl SearchConnectivity {
    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn iter(&self) -> SearchConnectivityIter {
        SearchConnectivityIter(self.0.iter())
    }
}

impl FromIterator<(usize, usize)> for SearchConnectivity {
    fn from_iter<T: IntoIterator<Item = (usize, usize)>>(iter: T) -> Self {
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

impl FromParallelIterator<(usize, usize)> for SearchConnectivity {
    fn from_par_iter<I>(par_iter: I) -> Self
    where
        I: IntoParallelIterator<Item = (usize, usize)>,
    {
        let v = par_iter.into_par_iter().collect::<Vec<_>>();
        let mut res = Self::default();
        for (i, j) in v.iter().cloned() {
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

//========================================================
pub struct DistanceSearcherSingle<I> {
    grid: Grid<Vec<I>>,
    cutoff: f32,
    serial_limit: usize,
}

/*
impl DistanceSearcherSingle<(usize,&Pos,f32)> {
    pub fn search_vdw<T, C>(&self) -> C
    where
        T: SearchOutputType + Send + Sync + Clone,
        C: FromParallelIterator<T> + FromIterator<T> + Send,
    {
        let cutoff2 = self.cutoff.powi(2);
        // Get periodic or non-periodic distance function
        match self.grid.pbc.as_ref() {
            Some(pbc) => self.search_internal(&pbc.dims, |p1, p2| {
                let d2 = pbc.pbox.distance_squared(p1, p2, &pbc.dims);
                if d2 <= cutoff2 {
                    (true, d2.sqrt())
                } else {
                    (false, 0.0)
                }
            }),
            None => self.search_internal(&PBC_NONE, |p1, p2| {
                let d2 = (p1 - p2).norm_squared();
                if d2 <= cutoff2 {
                    (true, d2.sqrt())
                } else {
                    (false, 0.0)
                }
            }),
        }
    }
}
*/

impl<I: GridItem + Sync + Send> DistanceSearcherSingle<I> {
    pub fn new<'a>(
        cutoff: f32,
        data: impl Iterator<Item = I>,
        lower: &Vector3f,
        upper: &Vector3f,
    ) -> Self {
        // Create grid
        let mut grid = Grid::from_cutoff_and_min_max(cutoff, lower, upper);
        grid.populate(data, lower, upper);
        // Create an instance
        Self {
            grid,
            cutoff,
            serial_limit: SERIAL_LIMIT,
        }
    }

    pub fn new_periodic<'a>(
        cutoff: f32,
        data: impl Iterator<Item = I>,
        box_: &PeriodicBox,
        periodic_dims: &PbcDims,
    ) -> Self {
        // Create grid
        let mut grid = Grid::from_cutoff_and_box(cutoff, box_);
        grid.populate_periodic(data, box_, &periodic_dims);
        // Create an instance
        Self {
            grid,
            cutoff,
            serial_limit: SERIAL_LIMIT,
        }
    }

    pub fn set_serial_limit(&mut self, limit: usize) {
        self.serial_limit = limit;
    }

    // Main search function
    pub fn search<T, C>(&self) -> C
    where
        T: SearchOutputType + Send + Sync + Clone,
        C: FromParallelIterator<T> + FromIterator<T> + Send,
    {
        let cutoff2 = self.cutoff.powi(2);
        // Get periodic or non-periodic distance function
        match self.grid.pbc.as_ref() {
            Some(pbc) => self.search_internal(&pbc.dims, |p1, p2| {
                let d2 = pbc.pbox.distance_squared(p1, p2, &pbc.dims);
                if d2 <= cutoff2 {
                    (true, d2.sqrt())
                } else {
                    (false, 0.0)
                }
            }),
            None => self.search_internal(&PBC_NONE, |p1, p2| {
                let d2 = (p1 - p2).norm_squared();
                if d2 <= cutoff2 {
                    (true, d2.sqrt())
                } else {
                    (false, 0.0)
                }
            }),
        }
    }

    fn search_internal<T, C, F>(&self, periodic_dims: &PbcDims, accept_func: F) -> C
    where
        T: SearchOutputType + Send + Sync + Clone,
        C: FromParallelIterator<T> + FromIterator<T> + Send,
        F: Fn(&Pos, &Pos) -> (bool, f32) + Send + Sync,
    {
        let nt = rayon::current_num_threads();
        if self.grid.dim()[0] >= nt && self.grid.dim()[0] >= self.serial_limit {
            // Parallel
            let chunk_size = self.grid.dim()[0] / nt;

            // List of result buffers
            let mut buffers = vec![Vec::<T>::with_capacity(INIT_BUF_SIZE); nt];

            buffers.par_iter_mut().enumerate().for_each(|(ch, buf)| {
                // Create a cell pair iterator
                for pair in CellPairIter::new_chunk(
                    &self.grid.dim(),
                    &periodic_dims,
                    ch * chunk_size,
                    if ch < nt - 1 {
                        (ch + 1) * chunk_size
                    } else {
                        self.grid.dim()[0]
                    },
                ) {
                    self.search_cell_pair(pair, &accept_func, buf);
                }
            });

            buffers.into_iter().flatten().collect::<C>()
        } else {
            // Serial
            let mut buf = Vec::<T>::with_capacity(INIT_BUF_SIZE);
            for pair in CellPairIter::new(&self.grid.dim(), &periodic_dims) {
                self.search_cell_pair(pair, &accept_func, &mut buf);
            }
            buf.into_iter().collect::<C>()
        }
    }

    fn search_cell_pair<T: SearchOutputType>(
        &self,
        pair: CellPair,
        accept_func: impl Fn(&Pos, &Pos) -> (bool, f32),
        found: &mut Vec<T>,
    ) {
        let n1 = self.grid[&pair.c1].len();
        let n2 = self.grid[&pair.c2].len();

        // Nothing to do if cell is empty
        if n1 * n2 == 0 {
            return;
        }

        if pair.c1 == pair.c2 {
            // Same cell
            for i in 0..n1 - 1 {
                let at1 = &self.grid[&pair.c1][i];
                for j in i + 1..n1 {
                    let at2 = &self.grid[&pair.c1][j];
                    if let (true, d) = accept_func(&at1.get_pos(), &at2.get_pos()) {
                        found.push(T::from_search_results(at1.get_id(), at2.get_id(), d));
                    }
                }
            }
        } else {
            // Different cells
            for i in 0..n1 {
                let at1 = &self.grid[&pair.c1][i];
                for j in 0..n2 {
                    let at2 = &self.grid[&pair.c2][j];
                    if let (true, d) = accept_func(&at1.get_pos(), &at2.get_pos()) {
                        found.push(T::from_search_results(at1.get_id(), at2.get_id(), d));
                    }
                }
            }
        }
    }
}

//=============================================================================
pub struct DistanceSearcherDouble<I> {
    grid1: Grid<Vec<I>>,
    grid2: Grid<Vec<I>>,
    cutoff: f32,
    serial_limit: usize,
}

impl<I: GridItem + Send + Sync> DistanceSearcherDouble<I> {
    pub fn new<'a>(
        cutoff: f32,
        data1: impl Iterator<Item = I>,
        data2: impl Iterator<Item = I>,
        lower: &Vector3f,
        upper: &Vector3f,
    ) -> Self {
        // Create grids
        let mut grid1 = Grid::from_cutoff_and_min_max(cutoff, lower, upper);
        let mut grid2 = Grid::new(grid1.dim());

        grid1.populate(data1, lower, upper);
        grid2.populate(data2, lower, upper);
        // Create an instance
        Self {
            grid1,
            grid2,
            cutoff,
            serial_limit: SERIAL_LIMIT,
        }
    }

    pub fn new_periodic<'a>(
        cutoff: f32,
        data1: impl Iterator<Item = I>,
        data2: impl Iterator<Item = I>,
        box_: &PeriodicBox,
        periodic_dims: &PbcDims,
    ) -> Self {
        // Create grids
        let mut grid1 = Grid::from_cutoff_and_box(cutoff, box_);
        let mut grid2 = Grid::new(grid1.dim());

        grid1.populate_periodic(data1, box_, &periodic_dims);
        grid2.populate_periodic(data2, box_, &periodic_dims);
        // Create an instance
        Self {
            grid1,
            grid2,
            cutoff,
            serial_limit: SERIAL_LIMIT,
        }
    }

    // Main search function
    pub fn search<T, C>(&self) -> C
    where
        T: SearchOutputType + Send + Clone,
        C: FromIterator<T>,
    {        
        let cutoff2 = self.cutoff.powi(2);
        // Get periodic or non-periodic distance function
        match self.grid1.pbc.as_ref() {
            Some(pbc) => self.search_internal(&pbc.dims, |p1, p2| {
                let d2 = pbc.pbox.distance_squared(p1, p2, &pbc.dims);
                if d2 <= cutoff2 {
                    (true, d2.sqrt())
                } else {
                    (false, 0.0)
                }
            }),
            None => self.search_internal(&PBC_NONE, |p1, p2| {
                let d2 = (p1 - p2).norm_squared();
                if d2 <= cutoff2 {
                    (true, d2.sqrt())
                } else {
                    (false, 0.0)
                }
            }),
        }
    }

    fn search_internal<T, C, F>(&self, periodic_dims: &PbcDims, accept_func: F) -> C
    where
        T: SearchOutputType + Send + Clone,
        C: FromIterator<T>,
        F: Fn(&Pos, &Pos) -> (bool, f32) + Send + Sync,
    {
        let nt = rayon::current_num_threads();
        if self.grid1.dim()[0] >= nt && self.grid1.dim()[0] >= self.serial_limit {
            // Parallel
            let chunk_size = self.grid1.dim()[0] / nt;

            // List of result buffers
            let mut buffers = vec![Vec::<T>::with_capacity(INIT_BUF_SIZE); nt];

            buffers.par_iter_mut().enumerate().for_each(|(ch, buf)| {
                // Create a cell pair iterator
                for pair in CellPairIter::new_chunk(
                    &self.grid1.dim(),
                    &periodic_dims,
                    ch * chunk_size,
                    if ch < nt - 1 {
                        (ch + 1) * chunk_size
                    } else {
                        self.grid1.dim()[0]
                    },
                ) {
                    // We first search for pair c1->grid1; c2->grid2
                    // then for the swapped pair c2->grid1; c1->grid2
                    // grid1 always goes first
                    let swapped_pair = pair.swaped();
                    self.search_cell_pair(pair, &accept_func, buf);
                    self.search_cell_pair(swapped_pair, &accept_func, buf);
                }
            });

            buffers.into_iter().flatten().collect::<C>()
        } else {
            // Serial
            let mut buf = Vec::<T>::with_capacity(INIT_BUF_SIZE);
            for pair in CellPairIter::new(&self.grid1.dim(), &periodic_dims) {
                let swapped_pair = pair.swaped();
                self.search_cell_pair(pair, &accept_func, &mut buf);
                self.search_cell_pair(swapped_pair, &accept_func, &mut buf);
            }
            buf.into_iter().collect::<C>()
        }
    }

    fn search_cell_pair<T: SearchOutputType>(
        &self,
        pair: CellPair,
        accept_func: impl Fn(&Pos, &Pos) -> (bool,f32),
        found: &mut Vec<T>,
    ) {
        let n1 = self.grid1[&pair.c1].len();
        let n2 = self.grid2[&pair.c2].len();

        // Nothing to do if cell is empty
        if n1 * n2 == 0 {
            return;
        }        

        for i in 0..n1 {
            let at1 = &self.grid1[&pair.c1][i];
            for j in 0..n2 {
                let at2 = &self.grid2[&pair.c2][j];
                if let (true,d) = accept_func(&at1.get_pos(), &at2.get_pos()) {                
                    found.push(T::from_search_results(at1.get_id(), at2.get_id(), d));
                }
            }
        }
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
/*

#[test]
fn test_single_periodic() {
    use crate::io::*;
    use crate::core::PBC_FULL;

    let mut r = FileHandler::open("tests/protein.pdb").unwrap();
    let st = r.read_state().unwrap().unwrap();

    let searcher = DistanceSearcherSingle::new_periodic(
        0.3,
        (0..st.coords.len()).map(|i| (i,&st.coords[i])),
        st.box_.as_ref().unwrap(),
        &PBC_FULL,
    );
    let found: SearchConnectivity = searcher.search();
    println!("{:?}", found.len())
}

#[test]
fn test_single_non_periodic() {
    use crate::io::*;

    let mut r = FileHandler::open("tests/protein.pdb").unwrap();
    let st = r.read_state().unwrap().unwrap();

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
    use crate::core::PBC_FULL;

    let mut r = FileHandler::open("tests/protein.pdb").unwrap();
    let st = r.read_state().unwrap().unwrap();

    let searcher = DistanceSearcherDouble::new_periodic(
        0.3,
        (0..st.coords.len()).map(|i| (i,&st.coords[i])),
        (0..st.coords.len()).map(|i| (i,&st.coords[i])),
        &st.box_.as_ref().unwrap(),
        &PBC_FULL,
    );
    let found: Vec<usize> = searcher.search();
    println!("{:?}", found.len())
}
*/
