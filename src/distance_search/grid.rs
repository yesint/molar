use crate::core::{PbcDims, PeriodicBox, Pos, Vector3f};

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
pub struct Grid<T>(ndarray::Array3<T>);

pub trait IdPosIterator<'a>: ExactSizeIterator<Item = (usize, &'a Pos)> {}
impl<'a, T> IdPosIterator<'a> for T where T: ExactSizeIterator<Item = (usize, &'a Pos)> {}

impl<T> Grid<T>
where
    T: Default,
{
    pub fn new(sz: [usize;3]) -> Self {
        Self(ndarray::Array3::<T>::from_shape_simple_fn(
            sz,
            || T::default(),
        ))
    }

    pub fn dim(&self) -> [usize;3] {
        let d = self.0.raw_dim();
        [d[0],d[1],d[2]]
    }
}

impl Grid<GridCellData> {
    pub fn populate<'a>(&mut self, id_pos: impl IdPosIterator<'a>, lower: &Vector3f, upper: &Vector3f) {
        let dim = self.0.raw_dim();
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
            self.0[ind].add(id, pos);
        }
    }

    pub fn populate_periodic<'a>(
        &mut self,
        id_pos: impl IdPosIterator<'a>,
        box_: &PeriodicBox,
        pbc_dims: &PbcDims,
    ) {
        let dim = self.0.raw_dim();

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
            self.0[ind].add(id, pos);
        }
    }
}

// Implement indexing by nalgebra Vector3 for generic Grid
impl<T> std::ops::Index<&CellLoc> for Grid<T> {
    type Output = T;

    fn index(&self, index: &CellLoc) -> &Self::Output {
        &self.0[(index.x, index.y, index.z)]
    }
}

impl<T> std::ops::IndexMut<&CellLoc> for Grid<T> {
    fn index_mut(&mut self, index: &CellLoc) -> &mut Self::Output {
        &mut self.0[(index.x, index.y, index.z)]
    }
}

//============================================================================

#[test]
fn test_grid() {
    use crate::io::*;
    use std::iter::zip;
    let mut r = FileHandler::new_reader("tests/topol.tpr").unwrap();
    let st = r.read_next_state().unwrap().unwrap();

    let mut gr = Grid::new([10, 10, 10]);
    gr.populate_periodic(
        zip(0..st.coords.len(), st.coords.iter()),
        &st.box_,
        &[true, true, true],
    );
}
