use ndarray::Array3;
use std::iter::zip;
use crate::core::{PeriodicBox, PbcDims, Point3f, Vector3f};


#[derive(Debug,Clone,Default)]
struct GridCell {
    ids: Vec<usize>,
    coords: Vec<Point3f>,
}

impl GridCell {
    fn add(&mut self, id: usize, coord: &Point3f) {
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

impl Grid {
    fn new(x_sz: usize, y_sz: usize, z_sz: usize) -> Self {
        Self {
            data: Array3::<GridCell>::from_shape_simple_fn(
                (x_sz,y_sz,z_sz),
                || GridCell::default()
            )
        }
    }

    fn populate(&mut self, 
        ids: impl ExactSizeIterator<Item = usize>,
        coords: impl ExactSizeIterator<Item = Point3f>,
        lower: &Vector3f,
        upper: &Vector3f) 
    {
        let dim = self.data.raw_dim();
        let dim_sz = upper-lower;
        'outer: for (id,ref crd) in zip(ids,coords) {
            let mut ind = [0usize,0,0];
            for d in 0..3 {
                let n = (dim[d] as f32 * (crd[d]-lower[d])/dim_sz[d]).floor() as isize;
                if n<0 || n>=dim[d] as isize {
                    continue 'outer;
                } else {
                    ind[d] = n as usize;
                }
            }
            self.data[ind].add(id,crd);
        }
    }


    fn populate_periodic(&mut self, 
        ids: impl ExactSizeIterator<Item = usize>,
        coords: impl ExactSizeIterator<Item = Point3f>,
        box_: &PeriodicBox, 
        pbc_dims: PbcDims)
    {
        let dim = self.data.raw_dim();
        
        'outer: for (id,ref crd) in zip(ids,coords) {
            let wrapped = box_.wrap_vector_dims(&crd.coords,&pbc_dims);
            let rel = box_.to_box_coords(&wrapped);
            let mut ind = [0usize,0,0];
            for d in 0..3 {
                let mut n = (dim[d] as f32 * rel[d]).floor() as isize;
                // If dimension in not periodic and 
                // out of bounds - skip the point
                if pbc_dims[d]==0 && (n>dim[d] as isize || n<0) {
                    continue 'outer;
                }
                // Correct for possible minor numeric errors
                ind[d] = n.clamp(0, dim[d] as isize - 1) as usize;
            }
            self.data[ind].add(id,crd);
        }
    }
}