use ndarray::Array3;
use std::iter::zip;
use crate::core::{PeriodicBox, PbcDims, Pos, Vector3f, IndexIterator, PosIterator};


#[derive(Debug,Clone,Default)]
struct GridCell {
    ids: Vec<usize>,
    coords: Vec<Pos>,
}

impl GridCell {
    fn add(&mut self, id: usize, coord: &Pos) {
        self.ids.push(id);
        self.coords.push(coord.clone());
    }

    fn len(&self) -> usize {
        self.ids.len()
    }
}

#[derive(Debug,Clone)]
struct Grid {
    data: Array3<GridCell>,
}

pub trait IdPosIterator<'a>: ExactSizeIterator<Item = (usize,&'a Pos)> {}
impl<'a,T> IdPosIterator<'a> for T where T: ExactSizeIterator<Item = (usize,&'a Pos)> {}

impl Grid {
    fn new(x_sz: usize, y_sz: usize, z_sz: usize) -> Self {
        Self {
            data: Array3::<GridCell>::from_shape_simple_fn(
                (x_sz,y_sz,z_sz),
                || GridCell::default()
            )
        }
    }

    fn populate<'a>(&mut self, 
        id_pos: impl IdPosIterator<'a>,
        lower: &Vector3f,
        upper: &Vector3f) 
    {
        let dim = self.data.raw_dim();
        let dim_sz = upper-lower;
        'outer: for (id, pos) in id_pos {
            let mut ind = [0usize,0,0];
            for d in 0..3 {
                let n = (dim[d] as f32 * (pos[d]-lower[d])/dim_sz[d]).floor() as isize;
                if n<0 || n>=dim[d] as isize {
                    continue 'outer;
                } else {
                    ind[d] = n as usize;
                }
            }
            self.data[ind].add(id,pos);
        }
    }


    fn populate_periodic<'a>(&mut self, 
        id_pos: impl IdPosIterator<'a>,
        box_: &PeriodicBox, 
        pbc_dims: PbcDims)
    {
        let dim = self.data.raw_dim();
        
        'outer: for (id, pos) in id_pos {
            let wrapped = box_.wrap_vector_dims(&pos.coords,&pbc_dims);
            let rel = box_.to_box_coords(&wrapped);
            let mut ind = [0usize,0,0];
            for d in 0..3 {
                let n = (dim[d] as f32 * rel[d]).floor() as isize;
                // If dimension in not periodic and 
                // out of bounds - skip the point
                if pbc_dims[d]==0 && (n>dim[d] as isize || n<0) {
                    continue 'outer;
                }
                // Correct for possible minor numeric errors
                ind[d] = n.clamp(0, dim[d] as isize - 1) as usize;
            }
            self.data[ind].add(id,pos);
        }
    }
}

#[test]
fn test_grid() {
    use crate::io::*;
    use std::iter::zip;
    let mut r = FileHandler::new_reader("tests/topol.tpr").unwrap();
    let st = r.read_next_state().unwrap().unwrap();

    let mut gr = Grid::new(10,10,10);
    gr.populate_periodic(
        zip(0..st.coords.len(),st.coords.iter()),
        &st.box_, [1u8,1,1]
    );
}