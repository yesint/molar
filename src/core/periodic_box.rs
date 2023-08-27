use nalgebra::{Matrix3, Vector3, Point3};
use anyhow::{Result,bail};

pub struct PeriodicBox {
    matrix: Matrix3<f32>,
    inv: Matrix3<f32>,
    is_rectangular: bool,
}

impl PeriodicBox {
    fn new(matrix: Matrix3<f32>) -> Result<Self> {
        // Sanity check
        for col in matrix.column_iter() {
            if col.norm() == 0.0 { bail!("Three non-zero periodic box vector required! Given {}",matrix) }
        }

        let is_rectangular =
        if matrix[(0,1)] != 0.0 || matrix[(0,2)] != 0.0 
        || matrix[(1,0)] != 0.0 || matrix[(1,2)] != 0.0 
        || matrix[(2,0)] != 0.0 || matrix[(2,1)] != 0.0 { false } else { true };

        Ok(Self{
            matrix,
            inv: matrix.try_inverse().unwrap(),
            is_rectangular
        })
    }

    fn new_from_vectors_angles(a: f32, b: f32, c: f32, alpha: f32, beta: f32, gamma: f32) -> Result<Self> {
        let mut m = Matrix3::<f32>::zeros();
    
        if a==0.0 || b==0.0 || c==0.0 {
            bail!("Box vectors have to be non-zero! ({} {} {})",a,b,c);
        }

        if alpha<60.0 || b<60.0 || c<60.0 {
            bail!("Box angles should be >= 60 deg! ({} {} {})",alpha,beta,gamma);
        }

        m[(0,0)] = a;
    
        if alpha != 90.0 || beta != 90.0 || gamma != 90.0 {
            let cosa = if alpha != 90.0 {
                alpha.to_radians().cos()
            } else {
                0.0
            };
            let cosb = if beta != 90.0 {
                beta.to_radians().cos()
            } else {
                0.0
            };
            let (sing, cosg) = if gamma != 90.0 {
                gamma.to_radians().sin_cos()
            } else {
                (1.0, 0.0)
            };
            m[(0,1)] = b * cosg;
            m[(1,1)] = b * sing;
            m[(0,2)] = c * cosb;
            m[(1,2)] = c * (cosa - cosb * cosg) / sing;
            m[(2,2)] = (c * c - m[(0,2)].powf(2.0) - m[(1,2)].powf(2.0)).sqrt();
        } else {
            m[(1,1)] = b;
            m[(2,2)] = c;
        }
        
        Self::new(m)
    }

    fn wrap_vector(&self, vec: Vector3<f32>) -> Vector3<f32> {
        if self.is_rectangular {
            return vec.clone();
        } else {
            // Get vector in box fractional coordinates
            let mut box_vec = self.inv*vec;
            for i in 0..3 {
                box_vec[i] -= box_vec[i].round();
            }
            return self.matrix * box_vec;
        }
    }


    fn wrap_vector_dims(&self, vec: Vector3<f32>, pbc_dims: &[i32;3]) -> Vector3<f32> {
        if self.is_rectangular {
            return vec.clone();
        } else {
            // Get vector in box fractional coordinates
            let mut box_vec = self.inv*vec;
            for i in 0..3 {
                if pbc_dims[i] !=0 {
                    box_vec[i] -= box_vec[i].round();
                }
            }
            return self.matrix * box_vec;
        }
    }

    fn wrap_point(&self, point: Point3<f32>) -> Point3<f32> {
        if self.is_rectangular {
            return point;
        } else {
            // Get vector in box fractional coordinates
            let mut box_vec = self.inv*point;
            for i in 0..3 {
                box_vec[i] -= box_vec[i].round();
            }
            return self.matrix * box_vec;
        }
    }



    fn closest_image(&self, point: &Point3<f32>, target: &Point3<f32>) -> Point3<f32> {
        target + self.wrap_vector(point-target)
    }

    fn closest_image_dims(&self, point: &Point3<f32>, target: &Point3<f32>, pbc_dims: &[i32;3]) -> Point3<f32> {
        target + self.wrap_vector_dims(point-target, pbc_dims)
    }


}