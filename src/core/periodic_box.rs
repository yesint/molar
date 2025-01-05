use crate::core::{Matrix3f, Pos, Vector3f};
use nalgebra::Const;
use thiserror::Error;

#[derive(Debug, Default, Clone)]
pub struct PeriodicBox {
    matrix: Matrix3f,
    inv: Matrix3f,
}

#[derive(Debug,PartialEq,Clone,Copy)]
pub struct PbcDims(u8);

impl PbcDims {
    pub fn set_dim(&mut self, n: usize, val: bool) {
        if n>2 {
            panic!("pbc has only 3 dimentions")
        }
        if val {
            self.0 |= 1 << n;
        } else {
            self.0 &= !(1 << n);
        }
    }

    pub fn new(x: bool, y: bool, z: bool) -> Self {
        let mut ret = Self(0);
        ret.set_dim(0, x);
        ret.set_dim(1, y);
        ret.set_dim(2, z);
        ret
    }

    pub fn get_dim(&self, n: usize) -> bool {
        if n>2 {
            panic!("pbc has only 3 dimentions")
        }
        (self.0 & (1 << n)) != 0
    }
    
    pub fn any(&self) -> bool {
        (self.0 & (1 << 0)) != 0
        ||
        (self.0 & (1 << 1)) != 0
        ||
        (self.0 & (1 << 2)) != 0
    }
}

pub const PBC_FULL: PbcDims = PbcDims(0b0000_0111);
pub const PBC_NONE: PbcDims = PbcDims(0b0000_0000);

#[derive(Error,Debug)]
pub enum PeriodicBoxError {
    #[error("pbc operation withon periodic box")]
    NoPbc,

    #[error("zero length box vector")]
    ZeroLengthVector,
    
    #[error("inverse failed")]
    InverseFailed,
    
    #[error("box angle is <60 deg")]
    AngleTooSmall,
}

impl PeriodicBox {
    pub fn from_matrix<S>(matrix: nalgebra::Matrix<f32,Const<3>,Const<3>,S>) -> Result<Self, PeriodicBoxError> 
    where S: nalgebra::storage::Storage<f32, Const<3>, Const<3>>
    {
        // Sanity check
        for col in matrix.column_iter() {
            if col.norm() == 0.0 {
                Err(PeriodicBoxError::ZeroLengthVector)?
            }
        }

        Ok(Self {
            matrix: matrix.clone_owned(),
            inv: matrix
                .try_inverse()
                .ok_or_else(|| PeriodicBoxError::InverseFailed)?,
        })
    }

    pub fn from_vectors_angles(
        a: f32,
        b: f32,
        c: f32,
        alpha: f32,
        beta: f32,
        gamma: f32,
    ) -> Result<Self, PeriodicBoxError> {
        let mut m = Matrix3f::zeros();

        if a == 0.0 || b == 0.0 || c == 0.0 {
            Err(PeriodicBoxError::ZeroLengthVector)?;
        }

        if alpha < 60.0 || beta < 60.0 || gamma < 60.0 {
            Err(PeriodicBoxError::AngleTooSmall)?;
        }

        m[(0, 0)] = a;

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
            m[(0, 1)] = b * cosg;
            m[(1, 1)] = b * sing;
            m[(0, 2)] = c * cosb;
            m[(1, 2)] = c * (cosa - cosb * cosg) / sing;
            m[(2, 2)] = (c * c - m[(0, 2)].powf(2.0) - m[(1, 2)].powf(2.0)).sqrt();
        } else {
            m[(1, 1)] = b;
            m[(2, 2)] = c;
        }

        Self::from_matrix(m)
    }

    pub fn to_vectors_angles(&self) -> (Vector3f, Vector3f) {
        let mut vectors = Vector3f::zeros();
        let mut angles = Vector3f::zeros();

        let vx = self.matrix.column(0);
        let vy = self.matrix.column(1);
        let vz = self.matrix.column(2);

        angles[0] = if vy.norm_squared() * vz.norm_squared() != 0.0 {
            vy.angle(&vz).to_degrees()
        } else {
            90.0
        };

        angles[1] = if vx.norm_squared() * vz.norm_squared() != 0.0 {
            vx.angle(&vz).to_degrees()
        } else {
            90.0
        };

        angles[2] = if vx.norm_squared() * vy.norm_squared() != 0.0 {
            vx.angle(&vy).to_degrees()
        } else {
            90.0
        };

        vectors[0] = vx.norm();
        vectors[1] = vy.norm();
        vectors[2] = vz.norm();

        (vectors, angles)
    }

    #[inline(always)]
    pub fn shortest_vector<S>(&self, vec: &nalgebra::Vector<f32,Const<3>,S>) -> Vector3f 
    where S: nalgebra::storage::Storage<f32, Const<3>>,
    {
        // Get vector in box fractional coordinates
        let mut box_vec = self.inv * vec;
        box_vec.apply(|v| {
            *v -= v.round();
        });
        return self.matrix * box_vec;
    }

    #[inline(always)]
    pub fn shortest_vector_dims<S>(&self, vec: &nalgebra::Vector<f32,Const<3>,S>, pbc_dims: PbcDims) -> Vector3f 
    where S: nalgebra::storage::Storage<f32, Const<3>>,
    {
        // Get vector in box fractional coordinates
        let mut box_vec = self.inv * vec;
        for i in 0..3 {
            if pbc_dims.get_dim(i) {
                box_vec[i] -= box_vec[i].round();
            }
        }
        return self.matrix * box_vec;
    }

    #[inline(always)]
    pub fn closest_image(&self, point: &Pos, target: &Pos) -> Pos {
        target + self.shortest_vector(&(point - target))
    }

    #[inline(always)]
    pub fn closest_image_dims(&self, point: &Pos, target: &Pos, pbc_dims: PbcDims) -> Pos {
        target + self.shortest_vector_dims(&(point - target), pbc_dims)
    }

    #[inline(always)]
    pub fn get_matrix(&self) -> Matrix3f {
        self.matrix
    }

    #[inline(always)]
    pub fn to_box_coords<S>(&self, vec: &nalgebra::Vector<f32,Const<3>,S>) -> Vector3f
    where S: nalgebra::storage::Storage<f32, Const<3>>,
    {
        self.inv * vec
    }

    #[inline(always)]
    pub fn is_inside(&self, point: &Pos) -> bool {
        let v = self.inv * point;
        v[0]<1.0 && v[1]<1.0 && v[2]<1.0 
        && v[0]>=0.0 && v[1]>=0.0 && v[2]>=0.0
    }

    #[inline(always)]
    pub fn to_lab_coords<S>(&self, vec: &nalgebra::Vector<f32,Const<3>,S>) -> Vector3f 
    where S: nalgebra::storage::Storage<f32, Const<3>>,
    {
        self.matrix * vec
    }

    #[inline(always)]
    pub fn get_box_extents(&self) -> Vector3f {
        Vector3f::from_iterator(self.matrix.column_iter().map(|c| c.norm()))
    }

    pub fn get_lab_extents(&self) -> Vector3f {
        Vector3f::new(
            self.matrix[(0, 0)] + self.matrix[(0, 1)] + self.matrix[(0, 2)],
            self.matrix[(1, 0)] + self.matrix[(1, 1)] + self.matrix[(1, 2)],
            self.matrix[(2, 0)] + self.matrix[(2, 1)] + self.matrix[(2, 2)],
        )
    }

    #[inline(always)]
    pub fn distance_squared(&self, p1: &Pos, p2: &Pos, pbc_dims: PbcDims) -> f32 {
        self.shortest_vector_dims(&(p2 - p1), pbc_dims).norm_squared()
    }

    #[inline(always)]
    pub fn distance(&self, p1: &Pos, p2: &Pos, pbc: PbcDims) -> f32 {
        self.distance_squared(p1, p2, pbc).sqrt()
    }

    pub fn is_triclinic(&self) -> bool {
        self.matrix[(0,1)]!=0.0 || self.matrix[(0,2)]!=0.0 ||
        self.matrix[(1,0)]!=0.0 || self.matrix[(1,2)]!=0.0 ||
        self.matrix[(2,0)]!=0.0 || self.matrix[(2,1)]!=0.0
    }
    
    pub(crate) fn scale_vectors(&mut self, scale_factors: [f32; 3]) -> Result<(),PeriodicBoxError> {
        for c in 0..3 {
            self.matrix.column_mut(c).apply(|el| *el *= scale_factors[c]);
        }
        self.inv = self.matrix
                .try_inverse()
                .ok_or_else(|| PeriodicBoxError::InverseFailed)?;
        Ok(())
    }

    #[inline(always)]
    pub fn wrap_point(&self, p: &Pos) -> Pos {
        // Get vector in box fractional coordinates
        let mut bv = self.inv * p;
        for i in 0..3 {
            bv[i] = bv[i].fract();
            if bv[i] < 0.0 {
                bv[i] = 1.0 - bv[i];
            }
        }
        return self.matrix * bv;
    }
}

#[cfg(test)]
mod tests {
    use super::PeriodicBox;

    #[test]
    #[should_panic]
    fn invalid_from_vec_ang() {
        let _b = PeriodicBox::from_vectors_angles(
            10.0,0.2,15.0, 90.0, 9.0, 90.0
        ).unwrap();
    }
}