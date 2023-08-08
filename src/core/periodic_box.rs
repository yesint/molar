use nalgebra::{Matrix3, Vector3};
pub struct PeriodicBox {
    matrix: Matrix3<f32>,
    inv: Matrix3<f32>,
}

impl PeriodicBox {
    fn new(matrix: Matrix3<f32>) -> Self {
        Self{
            matrix,
            inv: matrix.try_inverse().unwrap()
        }
    }

    fn new_from_vectors_angles(a: f32, b: f32, c: f32, alpha: f32, beta: f32, gamma: f32) -> Self {
        let mut m = Matrix3::<f32>::zeros();
    
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

}