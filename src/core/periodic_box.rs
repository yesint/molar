use crate::core::{Matrix3f, Pos, Vector3f};
use anyhow::{anyhow, bail, Result};

#[derive(Debug, Default, Clone)]
pub struct PeriodicBox {
    matrix: Matrix3f,
    inv: Matrix3f,
}

pub type PbcDims = [bool; 3];
pub const PBC_FULL: PbcDims = [true,true,true];
pub const PBC_NONE: PbcDims = [false,false,false];

impl PeriodicBox {
    pub fn from_matrix(matrix: Matrix3f) -> Result<Self> {
        // Sanity check
        for col in matrix.column_iter() {
            if col.norm() == 0.0 {
                bail!(
                    "Three non-zero periodic box vector required! Given {}",
                    matrix
                )
            }
        }

        Ok(Self {
            matrix,
            inv: matrix
                .try_inverse()
                .ok_or(anyhow!("Can't invert the pbc matrix: {}",matrix))?,
        })
    }

    pub fn from_vectors_angles(
        a: f32,
        b: f32,
        c: f32,
        alpha: f32,
        beta: f32,
        gamma: f32,
    ) -> Result<Self> {
        let mut m = Matrix3f::zeros();

        if a == 0.0 || b == 0.0 || c == 0.0 {
            bail!("Box vectors have to be non-zero! ({} {} {})", a, b, c);
        }

        if alpha < 60.0 || beta < 60.0 || gamma < 60.0 {
            bail!(
                "Box angles should be >= 60 deg! ({} {} {})",
                alpha,
                beta,
                gamma
            );
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
    pub fn shortest_vector(&self, vec: &Vector3f) -> Vector3f {
        // Get vector in box fractional coordinates
        let mut box_vec = self.inv * vec;
        box_vec.apply(|v| {
            *v -= v.round();
        });
        return self.matrix * box_vec;
    }

    #[inline(always)]
    pub fn shortest_vector_dims(&self, vec: &Vector3f, pbc_dims: &PbcDims) -> Vector3f {
        // Get vector in box fractional coordinates
        let mut box_vec = self.inv * vec;
        for i in 0..3 {
            if pbc_dims[i] {
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
    pub fn closest_image_dims(&self, point: &Pos, target: &Pos, pbc_dims: &PbcDims) -> Pos {
        target + self.shortest_vector_dims(&(point - target), pbc_dims)
    }

    #[inline(always)]
    pub fn get_matrix(&self) -> Matrix3f {
        self.matrix
    }

    #[inline(always)]
    pub fn to_box_coords(&self, vec: &Vector3f) -> Vector3f {
        self.inv * vec
    }

    #[inline(always)]
    pub fn to_lab_coords(&self, vec: &Vector3f) -> Vector3f {
        self.matrix * vec
    }

    #[inline(always)]
    pub fn get_extents(&self) -> Vector3f {
        Vector3f::from_iterator(self.matrix.column_iter().map(|c| c.norm()))
    }

    #[inline(always)]
    pub fn distance_squared(&self, p1: &Pos, p2: &Pos, pbc_dims: &PbcDims) -> f32 {
        self.shortest_vector_dims(&(p2 - p1), pbc_dims).norm_squared()
    }

    #[inline(always)]
    pub fn distance(&self, p1: &Pos, p2: &Pos, pbc: &PbcDims) -> f32 {
        self.distance_squared(p1, p2, pbc).sqrt()
    }
}
