use crate::prelude::*;
use nalgebra::{Const, Matrix, storage::Storage};
use thiserror::Error;

/// Periodic box allowing working with periodicity and computing periodic distances and images.
///
/// # Matrix convention
///
/// The box is stored as a 3x3 matrix whose **columns** are the box vectors `a`, `b`, `c`
/// (not the rows). This is the opposite of the row convention used by
/// mdtraj/MDAnalysis/`unitcell_vectors`. If you are porting a matrix from one of those
/// libraries, transpose it first (or assign each vector to the corresponding column
/// explicitly).
#[derive(Debug, Default, Clone)]
pub struct PeriodicBox {
    matrix: Matrix3f,
    inv: Matrix3f,
    // Precomputed lattice-vector corrections i*a + j*b + k*c for (i,j,k) ∈
    // {-1,0,1}^3 \ {(0,0,0)} that could shorten a displacement already
    // reduced to the primary parallelepiped. Empty for orthogonal boxes
    // (zero overhead on the hot path).
    tric_corrections: Vec<Vector3f>,
}

fn build_tric_corrections(m: &Matrix3f) -> Vec<Vector3f> {
    // Orthogonal fast path: every off-diagonal is zero => no correction needed.
    if m[(0, 1)] == 0.0
        && m[(0, 2)] == 0.0
        && m[(1, 0)] == 0.0
        && m[(1, 2)] == 0.0
        && m[(2, 0)] == 0.0
        && m[(2, 1)] == 0.0
    {
        return Vec::new();
    }
    let a: Vector3f = m.column(0).into();
    let b: Vector3f = m.column(1).into();
    let c: Vector3f = m.column(2).into();
    // After fractional reduction, |dx| is bounded by half the longest
    // space-diagonal of the primary parallelepiped. A shift s can shorten
    // a displacement of norm d only if |s| < 2d. Using the upper bound
    // on d as a conservative filter keeps only candidates that could ever
    // help — for well-formed boxes this prunes 26 → a handful.
    let half_diag = 0.5
        * (a + b + c)
            .norm()
            .max((a + b - c).norm())
            .max((a - b + c).norm())
            .max((-a + b + c).norm());
    let bound2 = (2.0 * half_diag).powi(2);
    let mut out = Vec::with_capacity(26);
    for i in -1_i32..=1 {
        for j in -1_i32..=1 {
            for k in -1_i32..=1 {
                if i == 0 && j == 0 && k == 0 {
                    continue;
                }
                let s: Vector3f = (i as Float) * a + (j as Float) * b + (k as Float) * c;
                if s.norm_squared() < bound2 {
                    out.push(s);
                }
            }
        }
    }
    out
}

/// Representation of periodic dimensions in 3D space.
#[derive(Debug,PartialEq,Clone,Copy)]
pub struct PbcDims(u8);

impl PbcDims {
    /// Sets the periodic boundary condition for a specific dimension.
    /// 
    /// # Arguments
    /// * `n` - Dimension index (0=x, 1=y, 2=z)
    /// * `val` - true to enable periodicity, false to disable
    /// 
    /// # Panics
    /// Panics if n > 2 (only 3 dimensions are supported)
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

    /// Creates a new PbcDims instance with specified periodicities.
    /// 
    /// # Arguments
    /// * `x` - Periodicity in x dimension
    /// * `y` - Periodicity in y dimension
    /// * `z` - Periodicity in z dimension
    /// 
    /// # Returns
    /// A new PbcDims instance with the specified periodicities
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

/// All dimentions are periodic
pub const PBC_FULL: PbcDims = PbcDims(0b0000_0111);
/// All dimentions are non-periodic
pub const PBC_NONE: PbcDims = PbcDims(0b0000_0000);

/// Errors related to periodic boxes and periodicity
#[derive(Error,Debug)]
pub enum PeriodicBoxError {
    #[error("pbc operation withon periodic box")]
    NoPbc,

    #[error("zero length box vector")]
    ZeroLengthVector,
    
    #[error("box matrix inverse failed")]
    InverseFailed,
    
    #[error("box angle is <60 deg")]
    AngleTooSmall,
}

impl PeriodicBox {
    /// Creates a new PeriodicBox from a 3x3 matrix representing box vectors.
    /// 
    /// # Arguments
    /// * `matrix` - 3x3 matrix where **columns** are the box vectors `a`, `b`, `c`.
    ///   Note the convention: mdtraj/MDAnalysis use rows as box vectors — transpose
    ///   their matrices before passing them here.
    ///
    /// # Errors
    /// Returns error if any vector has zero length or matrix is not invertible
    pub fn from_matrix<S>(matrix: Matrix<Float,Const<3>,Const<3>,S>) -> Result<Self, PeriodicBoxError>
    where S: Storage<Float, Const<3>, Const<3>>
    {
        // Sanity check
        for col in matrix.column_iter() {
            if col.norm() == 0.0 {
                Err(PeriodicBoxError::ZeroLengthVector)?
            }
        }

        let matrix = matrix.clone_owned();
        let inv = matrix
            .try_inverse()
            .ok_or_else(|| PeriodicBoxError::InverseFailed)?;
        let tric_corrections = build_tric_corrections(&matrix);
        Ok(Self {
            matrix,
            inv,
            tric_corrections,
        })
    }

    /// Creates a new PeriodicBox from box vectors lengths and angles between them.
    /// 
    /// # Arguments
    /// * `a`, `b`, `c` - Lengths of box vectors
    /// * `alpha` - Angle between b and c vectors (degrees)
    /// * `beta` - Angle between a and c vectors (degrees)
    /// * `gamma` - Angle between a and b vectors (degrees)
    /// 
    /// # Errors
    /// Returns error if any length is zero or any angle is less than 60 degrees
    pub fn from_vectors_angles(
        a: Float,
        b: Float,
        c: Float,
        alpha: Float,
        beta: Float,
        gamma: Float,
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

    /// Returns box vectors lengths and angles between them.
    /// 
    /// # Returns
    /// Tuple containing:
    /// - Vector of lengths (a, b, c)
    /// - Vector of angles in degrees (alpha, beta, gamma)
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

    /// Computes the shortest vector between two points considering periodicity.
    #[inline(always)]
    pub fn shortest_vector<S>(&self, vec: &nalgebra::Vector<Float,Const<3>,S>) -> Vector3f
    where S: Storage<Float, Const<3>>,
    {
        self.shortest_vector_dims(vec, PBC_FULL)
    }

    /// Computes the shortest vector between two points considering periodicity only in specified dimensions.
    #[inline(always)]
    pub fn shortest_vector_dims<S>(&self, vec: &nalgebra::Vector<Float,Const<3>,S>, pbc_dims: PbcDims) -> Vector3f
    where S: Storage<Float, Const<3>>,
    {
        // Step 1: fractional-coord reduction — collapses arbitrarily
        // far-apart pairs into the primary parallelepiped cell.
        let mut box_vec = self.inv * vec;
        for i in 0..3 {
            if pbc_dims.get_dim(i) {
                box_vec[i] -= box_vec[i].round();
            }
        }
        let start = self.matrix * box_vec;

        // Step 2: GROMACS-style triclinic correction. The list is empty
        // for orthogonal boxes, so orthogonal cells pay nothing beyond
        // a `Vec::is_empty` check. Partial PBC on a triclinic box is
        // ill-defined (triclinic vectors are not axis-aligned), so we
        // only apply the correction under full PBC.
        if self.tric_corrections.is_empty() || pbc_dims != PBC_FULL {
            return start;
        }
        let mut best = start;
        let mut best2 = start.norm_squared();
        for s in &self.tric_corrections {
            let cand = start + s;
            let n2 = cand.norm_squared();
            if n2 < best2 {
                best2 = n2;
                best = cand;
            }
        }
        best
    }

    /// Returns the closest periodic image of a point relative to a target.
    #[inline(always)]
    pub fn closest_image(&self, point: &Pos, target: &Pos) -> Pos {
        target + self.shortest_vector(&(point - target))
    }

    /// Returns the closest periodic image of a point relative to a target, considering periodicity only in specified dimensions.
    #[inline(always)]
    pub fn closest_image_dims(&self, point: &Pos, target: &Pos, pbc_dims: PbcDims) -> Pos {
        target + self.shortest_vector_dims(&(point - target), pbc_dims)
    }

    /// Returns the box matrix.
    #[inline(always)]
    pub fn get_matrix(&self) -> Matrix3f {
        self.matrix
    }

    /// Converts coordinates from lab frame to box frame (fractional coordinates).
    #[inline(always)]
    pub fn to_box_coords<S>(&self, vec: &nalgebra::Vector<Float,Const<3>,S>) -> Vector3f
    where S: Storage<Float, Const<3>>,
    {
        self.inv * vec
    }

    /// Checks if a point is inside the box (coordinates between 0 and 1 in box frame).
    #[inline(always)]
    pub fn is_inside(&self, point: &Pos) -> bool {
        let v = self.inv * point;
        v[0]<1.0 && v[1]<1.0 && v[2]<1.0 
        && v[0]>=0.0 && v[1]>=0.0 && v[2]>=0.0
    }

    /// Converts coordinates from box frame to lab frame.
    #[inline(always)]
    pub fn to_lab_coords<S>(&self, vec: &nalgebra::Vector<Float,Const<3>,S>) -> Vector3f 
    where S: Storage<Float, Const<3>>,
    {
        self.matrix * vec
    }

    /// Returns the lengths of box vectors.
    #[inline(always)]
    pub fn get_box_extents(&self) -> Vector3f {
        Vector3f::from_iterator(self.matrix.column_iter().map(|c| c.norm()))
    }

    /// Returns the maximum extents of the box in lab frame coordinates.
    pub fn get_lab_extents(&self) -> Vector3f {
        Vector3f::new(
            self.matrix[(0, 0)] + self.matrix[(0, 1)] + self.matrix[(0, 2)],
            self.matrix[(1, 0)] + self.matrix[(1, 1)] + self.matrix[(1, 2)],
            self.matrix[(2, 0)] + self.matrix[(2, 1)] + self.matrix[(2, 2)],
        )
    }

    /// Computes squared distance between two points considering periodic boundary conditions.
    #[inline(always)]
    pub fn distance_squared(&self, p1: &Pos, p2: &Pos, pbc_dims: PbcDims) -> Float {
        self.shortest_vector_dims(&(p2 - p1), pbc_dims).norm_squared()
    }

    /// Computes distance between two points considering periodic boundary conditions.
    #[inline(always)]
    pub fn distance(&self, p1: &Pos, p2: &Pos, pbc: PbcDims) -> Float {
        self.distance_squared(p1, p2, pbc).sqrt()
    }

    /// Checks if the box is triclinic (has non-orthogonal vectors).
    pub fn is_triclinic(&self) -> bool {
        self.matrix[(0,1)]!=0.0 || self.matrix[(0,2)]!=0.0 ||
        self.matrix[(1,0)]!=0.0 || self.matrix[(1,2)]!=0.0 ||
        self.matrix[(2,0)]!=0.0 || self.matrix[(2,1)]!=0.0
    }
    
    /// Scales box vectors by given factors.
    pub(crate) fn scale_vectors(&mut self, scale_factors: [Float; 3]) -> Result<(),PeriodicBoxError> {
        for c in 0..3 {
            self.matrix.column_mut(c).apply(|el| *el *= scale_factors[c]);
        }
        self.inv = self.matrix
                .try_inverse()
                .ok_or_else(|| PeriodicBoxError::InverseFailed)?;
        Ok(())
    }

    /// Wraps a point into the primary unit cell.
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

    #[inline(always)]
    pub fn wrap_vec<S>(&self, vec: &nalgebra::Vector<Float,Const<3>,S>) -> Vector3f 
    where S: Storage<Float, Const<3>>,
    {
        // Get vector in box fractional coordinates
        let mut bv = self.inv * vec;
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
    use crate::prelude::*;
    use super::PeriodicBox;

    const EPSILON: Float = 1e-6;

    fn assert_vec_eq(v1: &Vector3f, v2: &Vector3f) {
        assert!((v1 - v2).norm() < EPSILON, "Vectors not equal: {:?} != {:?}", v1, v2);
    }

    #[test]
    #[should_panic]
    fn invalid_from_vec_ang() {
        let _b = PeriodicBox::from_vectors_angles(
            10.0,0.2,15.0, 90.0, 9.0, 90.0
        ).unwrap();
    }

    #[test]
    fn test_shortest_vector_dims_no_pbc() {
        let box_matrix = Matrix3f::from_diagonal(&Vector3f::new(10.0, 10.0, 10.0));
        let pbox = PeriodicBox::from_matrix(box_matrix).unwrap();
        let test_vec = Vector3f::new(8.0, 8.0, 8.0);
        
        let result = pbox.shortest_vector_dims(&test_vec, PBC_NONE);
        assert_vec_eq(&result, &test_vec);
    }

    #[test]
    fn test_shortest_vector_dims_full_pbc() {
        let box_matrix = Matrix3f::from_diagonal(&Vector3f::new(10.0, 10.0, 10.0));
        let pbox = PeriodicBox::from_matrix(box_matrix).unwrap();
        let test_vec = Vector3f::new(8.0, 8.0, 8.0);
        
        let result = pbox.shortest_vector_dims(&test_vec, PBC_FULL);
        assert_vec_eq(&result, &Vector3f::new(-2.0, -2.0, -2.0));
    }

    #[test]
    fn test_shortest_vector_dims_x_only() {
        let box_matrix = Matrix3f::from_diagonal(&Vector3f::new(10.0, 10.0, 10.0));
        let pbox = PeriodicBox::from_matrix(box_matrix).unwrap();
        let test_vec = Vector3f::new(8.0, 8.0, 8.0);
        
        let pbc_x = PbcDims::new(true, false, false);
        let result = pbox.shortest_vector_dims(&test_vec, pbc_x);
        assert_vec_eq(&result, &Vector3f::new(-2.0, 8.0, 8.0));
    }

    #[test]
    fn test_shortest_vector_dims_xy_only() {
        let box_matrix = Matrix3f::from_diagonal(&Vector3f::new(10.0, 10.0, 10.0));
        let pbox = PeriodicBox::from_matrix(box_matrix).unwrap();
        let test_vec = Vector3f::new(8.0, 8.0, 8.0);
        
        let pbc_xy = PbcDims::new(true, true, false);
        let result = pbox.shortest_vector_dims(&test_vec, pbc_xy);
        assert_vec_eq(&result, &Vector3f::new(-2.0, -2.0, 8.0));
    }

    #[test]
    fn test_closest_image_no_pbc() {
        let box_matrix = Matrix3f::from_diagonal(&Vector3f::new(10.0, 10.0, 10.0));
        let pbox = PeriodicBox::from_matrix(box_matrix).unwrap();
        let point = Pos::new(8.0, 8.0, 8.0);
        let target = Pos::origin();
        
        let result = pbox.closest_image_dims(&point, &target, PBC_NONE);
        assert_vec_eq(&result.coords, &point.coords);
    }

    #[test]
    fn test_closest_image_full_pbc() {
        let box_matrix = Matrix3f::from_diagonal(&Vector3f::new(10.0, 10.0, 10.0));
        let pbox = PeriodicBox::from_matrix(box_matrix).unwrap();
        let point = Pos::new(8.0, 8.0, 8.0);
        let target = Pos::origin();
        
        let result = pbox.closest_image_dims(&point, &target, PBC_FULL);
        assert_vec_eq(&result.coords, &Pos::new(-2.0, -2.0, -2.0).coords);
    }

    #[test]
    fn test_closest_image_x_only() {
        let box_matrix = Matrix3f::from_diagonal(&Vector3f::new(10.0, 10.0, 10.0));
        let pbox = PeriodicBox::from_matrix(box_matrix).unwrap();
        let point = Pos::new(8.0, 8.0, 8.0);
        let target = Pos::origin();
        
        let pbc_x = PbcDims::new(true, false, false);
        let result = pbox.closest_image_dims(&point, &target, pbc_x);
        assert_vec_eq(&result.coords, &Pos::new(-2.0, 8.0, 8.0).coords);
    }

    #[test]
    fn test_closest_image_xy_only() {
        let box_matrix = Matrix3f::from_diagonal(&Vector3f::new(10.0, 10.0, 10.0));
        let pbox = PeriodicBox::from_matrix(box_matrix).unwrap();
        let point = Pos::new(8.0, 8.0, 8.0);
        let target = Pos::origin();

        let pbc_xy = PbcDims::new(true, true, false);
        let result = pbox.closest_image_dims(&point, &target, pbc_xy);
        assert_vec_eq(&result.coords, &Pos::new(-2.0, -2.0, 8.0).coords);
    }

    // Orthogonal boxes must carry no correction shifts — this guards the
    // hot-path early-return in shortest_vector_dims.
    #[test]
    fn test_orthogonal_has_no_tric_corrections() {
        let box_matrix = Matrix3f::from_diagonal(&Vector3f::new(10.0, 20.0, 30.0));
        let pbox = PeriodicBox::from_matrix(box_matrix).unwrap();
        assert!(pbox.tric_corrections.is_empty());
    }

    // Regression for issue #6: a GROMACS-legal triclinic box with off-diagonal
    // components exposes the former fractional-rounding bug. Box vectors here
    // match mdtraj's row-convention interpretation of the reporter's numpy
    // input: a=(10,0,0), b=(4,10,0), c=(-4,0,10). Points are the reporter's.
    // Brute-force and mdtraj/MDAnalysis agree on 5.353627 nm; the old
    // algorithm returned 5.597546 nm.
    #[test]
    fn test_triclinic_mdtraj_box_matches_brute_force() {
        // Rows-as-matrix notation: columns are the box vectors.
        let box_matrix = Matrix3f::new(
            10.0, 4.0, -4.0,
            0.0, 10.0,  0.0,
            0.0,  0.0, 10.0,
        );
        let pbox = PeriodicBox::from_matrix(box_matrix).unwrap();
        let p1 = Pos::new( 38.9214, 40.0078, -34.0795);
        let p2 = Pos::new(-26.6187, 40.8926,  30.9709);
        let d = pbox.distance(&p1, &p2, PBC_FULL);
        assert!(
            (d - 5.353627).abs() < 1e-3,
            "expected ~5.353627, got {d} (prior buggy value ~5.597546)"
        );
    }

    // Independent sanity check: a small triclinic box with a displacement
    // deliberately placed where independent-axis fractional rounding lands
    // in a neighboring image. Brute-force minimum is computed here too.
    #[test]
    fn test_triclinic_corner_matches_brute_force() {
        // A skewed box with a large c_x shear.
        let box_matrix = Matrix3f::new(
            6.0, 0.0, 3.0,
            0.0, 6.0, 3.0,
            0.0, 0.0, 6.0,
        );
        let pbox = PeriodicBox::from_matrix(box_matrix).unwrap();
        let dx = Vector3f::new(2.9, 2.9, 2.9);

        // Reference: exhaustive 5x5x5 brute-force minimum over lattice images.
        let a: Vector3f = box_matrix.column(0).into();
        let b: Vector3f = box_matrix.column(1).into();
        let c: Vector3f = box_matrix.column(2).into();
        let mut best = Float::INFINITY;
        for i in -2..=2i32 { for j in -2..=2i32 { for k in -2..=2i32 {
            let cand = dx + (i as Float)*a + (j as Float)*b + (k as Float)*c;
            best = best.min(cand.norm());
        }}}

        let got = pbox.shortest_vector(&dx).norm();
        assert!((got - best).abs() < 1e-5, "expected {best}, got {got}");
    }

    // Far-apart points (|dx| >> box) must still be handled correctly via the
    // initial fractional-coord reduction step.
    #[test]
    fn test_triclinic_far_apart_reduction() {
        let box_matrix = Matrix3f::new(
            10.0, 4.0, -4.0,
            0.0, 10.0,  0.0,
            0.0,  0.0, 10.0,
        );
        let pbox = PeriodicBox::from_matrix(box_matrix).unwrap();
        // p1 and p2 ~6 box-lengths apart.
        let p1 = Pos::new(0.1, 0.2, 0.3);
        let p2 = Pos::new(60.1, 0.2, 0.3);
        let d = pbox.distance(&p1, &p2, PBC_FULL);
        assert!(d < 1e-4, "pure a-vector shift should collapse to ~0, got {d}");
    }
}