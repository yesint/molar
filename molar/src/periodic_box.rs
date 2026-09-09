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
    // XY, XZ, YZ and XYZ lattices. Diagonal boxes need no search data.
    lattices: Option<Box<[ImageLattice; 4]>>,
}

/// Thin QR factorization of the enabled box vectors. The orthogonal complement
/// is not part of the distance objective, so partial PBC also works for points
/// arbitrarily far from the nonperiodic unit-cell faces.
#[derive(Debug, Clone)]
struct ImageLattice {
    basis: Matrix3f,
    q: Matrix3f,
    r: Matrix3f,
    inv_diag: Vector3f,
    ndim: usize,
    cartesian: bool,
    unique_radius2: Float,
    orthogonal: bool,
    corrections: Option<Vec<ImageCorrection>>,
}

#[derive(Debug, Clone)]
struct ImageCorrection {
    projected: Vector3f,
    coefficients: Vector3f,
    norm2: Float,
}

fn lattice_qr(
    basis: &nalgebra::Matrix3<f64>,
    ndim: usize,
) -> (nalgebra::Matrix3<f64>, nalgebra::Matrix3<f64>) {
    let mut q = nalgebra::Matrix3::<f64>::zeros();
    let mut r = nalgebra::Matrix3::<f64>::zeros();
    for j in 0..ndim {
        let mut v = basis.column(j).into_owned();
        for i in 0..j {
            r[(i, j)] = q.column(i).dot(&v);
            v -= q.column(i) * r[(i, j)];
        }
        r[(j, j)] = v.norm();
        q.set_column(j, &(v / r[(j, j)]));
    }
    (q, r)
}

impl ImageLattice {
    fn new(matrix: &Matrix3f, axes: &[usize]) -> Self {
        let ndim = axes.len();
        let mut basis = nalgebra::Matrix3::<f64>::zeros();
        for (j, &axis) in axes.iter().enumerate() {
            basis.set_column(j, &matrix.column(axis).into_owned().cast::<f64>());
        }
        let (mut q, mut r) = lattice_qr(&basis, ndim);
        // LLL size reduction and swaps preserve the enabled lattice. Reduction
        // controls search cost; correctness never depends on its completion.
        let mut k = 1;
        for _ in 0..128 {
            if k >= ndim {
                break;
            }
            for j in (0..k).rev() {
                let n = (r[(j, k)] / r[(j, j)]).round();
                if n != 0.0 {
                    let v = basis.column(k) - basis.column(j) * n;
                    basis.set_column(k, &v);
                    (q, r) = lattice_qr(&basis, ndim);
                }
            }
            let mu = r[(k - 1, k)] / r[(k - 1, k - 1)];
            if r[(k, k)].powi(2) >= (0.75 - mu * mu) * r[(k - 1, k - 1)].powi(2) {
                k += 1;
            } else {
                basis.swap_columns(k, k - 1);
                (q, r) = lattice_qr(&basis, ndim);
                k = k.saturating_sub(1).max(1);
            }
        }
        let basis = basis.cast::<Float>();
        let q = q.cast::<Float>();
        let r = r.cast::<Float>();
        let mut inv_diag = Vector3f::zeros();
        let mut min_diag = Float::INFINITY;
        for i in 0..axes.len() {
            inv_diag[i] = r[(i, i)].recip();
            min_diag = min_diag.min(r[(i, i)]);
        }
        let orthogonal = r[(0, 1)] == 0.0 && r[(0, 2)] == 0.0 && r[(1, 2)] == 0.0;
        Self {
            basis,
            q,
            r,
            inv_diag,
            ndim: axes.len(),
            orthogonal,
            corrections: if orthogonal {
                Some(Vec::new())
            } else {
                image_corrections(&r, ndim)
            },
            cartesian: axes.len() == 3 && q == Matrix3f::identity(),
            // Any nonzero lattice vector has a last nonzero coefficient. Its
            // corresponding QR component is at least min_diag in magnitude.
            unique_radius2: 0.25 * min_diag * min_diag * (1.0 - 32.0 * Float::EPSILON),
        }
    }

    #[inline]
    fn shortest(&self, vec: &Vector3f) -> Vector3f {
        let mut y = if self.cartesian {
            *vec
        } else {
            self.q.transpose() * vec
        };
        if y.norm_squared() < self.unique_radius2 {
            return *vec;
        }
        let mut shift = Vector3f::zeros();
        // Nearest-plane reduction gives a short initial image, including for
        // displacements many cells apart. It is not itself an exact solution.
        for i in (0..self.ndim).rev() {
            shift[i] = (y[i] * self.inv_diag[i]).round();
            for j in 0..=i {
                y[j] -= shift[i] * self.r[(j, i)];
            }
        }
        let start = vec - self.basis * shift;
        let mut best2 = y.norm_squared();
        if self.orthogonal || best2 < self.unique_radius2 {
            return start;
        }
        if let Some(corrections) = &self.corrections {
            let start2 = best2;
            let mut correction = Vector3f::zeros();
            for s in corrections {
                let dot = y.dot(&s.projected);
                let d2 = start2 + s.norm2 - 2.0 * dot.abs();
                if d2 < best2 {
                    best2 = d2;
                    correction = s.coefficients * dot.signum();
                    if best2 < self.unique_radius2 {
                        break;
                    }
                }
            }
            return start - self.basis * correction;
        }
        let mut best_shift = Vector3f::zeros();
        // Enumerate only the sphere that could improve the initial image.
        // For each Z and Y coefficient the nearest X coefficient is exact.
        let (zlo, zhi) = if self.ndim == 3 {
            image_bounds(y[2], best2.sqrt(), self.inv_diag[2])
        } else {
            (0, 0)
        };
        for iz in zlo..=zhi {
            let z = iz as Float;
            let dz = y[2] - z * self.r[(2, 2)];
            let z2 = dz * dz;
            if z2 > best2 {
                continue;
            }
            let cy = y[1] - z * self.r[(1, 2)];
            let (ylo, yhi) = image_bounds(cy, (best2 - z2).sqrt(), self.inv_diag[1]);
            for iy in ylo..=yhi {
                let b = iy as Float;
                let dy = cy - b * self.r[(1, 1)];
                let yz2 = z2 + dy * dy;
                if yz2 > best2 {
                    continue;
                }
                let cx = y[0] - b * self.r[(0, 1)] - z * self.r[(0, 2)];
                let a = (cx * self.inv_diag[0]).round();
                let dx = cx - a * self.r[(0, 0)];
                let d2 = yz2 + dx * dx;
                if d2 < best2 {
                    best2 = d2;
                    best_shift = Vector3f::new(a, b, z);
                    if best2 < self.unique_radius2 {
                        return start - self.basis * best_shift;
                    }
                }
            }
        }
        start - self.basis * best_shift
    }
}

/// Complete correction set for the nearest-plane brick. A useful shift obeys
/// |s| <= sqrt(sum(diagonal^2)) and |s|^2 < sum(diagonal_i * |s_i|).
/// Unlike a fixed 27-cell stencil, the coefficient bounds come from the inverse
/// basis. Store one of each +/- pair. Fall back to sphere enumeration if setup
/// or storage would be excessive for a poorly conditioned lattice.
fn image_corrections(r: &Matrix3f, ndim: usize) -> Option<Vec<ImageCorrection>> {
    let mut full = *r;
    for i in ndim..3 {
        full[(i, i)] = 1.0;
    }
    let inverse = full.try_inverse()?;
    let diameter = r.diagonal().norm();
    let mut bounds = [0_i32; 3];
    let mut count = 1_u64;
    for i in 0..ndim {
        let bound = (inverse.row(i).norm() * diameter).ceil();
        if !bound.is_finite() || bound > 32.0 {
            return None;
        }
        bounds[i] = bound as i32;
        count *= (2 * bounds[i] + 1) as u64;
    }
    if count > 4096 {
        return None;
    }
    let mut out = Vec::new();
    for i in -bounds[0]..=bounds[0] {
        for j in -bounds[1]..=bounds[1] {
            for k in -bounds[2]..=bounds[2] {
                if i < 0 || (i == 0 && j < 0) || (i == 0 && j == 0 && k <= 0) {
                    continue;
                }
                let coefficients = Vector3f::new(i as Float, j as Float, k as Float);
                let projected = r * coefficients;
                let norm2 = projected.norm_squared();
                let support = r.diagonal().dot(&projected.abs());
                if norm2 <= support * (1.0 + 32.0 * Float::EPSILON) {
                    out.push(ImageCorrection {
                        projected,
                        coefficients,
                        norm2,
                    });
                    if out.len() > 64 {
                        return None;
                    }
                }
            }
        }
    }
    out.sort_by(|a, b| a.norm2.total_cmp(&b.norm2));
    Some(out)
}

#[inline]
fn image_bounds(center: Float, radius: Float, inv_diag: Float) -> (i64, i64) {
    let c = center * inv_diag;
    let r = radius * inv_diag;
    // Keep candidates on sphere boundaries despite floating-point rounding.
    let margin = 16.0 * Float::EPSILON * (1.0 + c.abs() + r);
    (
        (c - r - margin).ceil() as i64,
        (c + r + margin).floor() as i64,
    )
}

/// Periodic lattice directions. Bits X, Y and Z select box vectors a, b and c,
/// respectively; for a triclinic box these are not Cartesian-axis constraints.
#[derive(Debug, PartialEq, Clone, Copy)]
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
        if n > 2 {
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
        if n > 2 {
            panic!("pbc has only 3 dimentions")
        }
        (self.0 & (1 << n)) != 0
    }

    pub fn any(&self) -> bool {
        (self.0 & (1 << 0)) != 0 || (self.0 & (1 << 1)) != 0 || (self.0 & (1 << 2)) != 0
    }
}

/// All dimentions are periodic
pub const PBC_FULL: PbcDims = PbcDims(0b0000_0111);
/// All dimentions are non-periodic
pub const PBC_NONE: PbcDims = PbcDims(0b0000_0000);

/// Errors related to periodic boxes and periodicity
#[derive(Error, Debug)]
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
    pub fn from_matrix<S>(
        matrix: Matrix<Float, Const<3>, Const<3>, S>,
    ) -> Result<Self, PeriodicBoxError>
    where
        S: Storage<Float, Const<3>, Const<3>>,
    {
        if !matrix.iter().all(|v| v.is_finite()) {
            return Err(PeriodicBoxError::InverseFailed);
        }
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
        let diagonal = (0..3).all(|i| (0..3).all(|j| i == j || matrix[(i, j)] == 0.0));
        let lattices = if diagonal {
            None
        } else {
            Some(Box::new([
                ImageLattice::new(&matrix, &[0, 1]),
                ImageLattice::new(&matrix, &[0, 2]),
                ImageLattice::new(&matrix, &[1, 2]),
                ImageLattice::new(&matrix, &[0, 1, 2]),
            ]))
        };
        Ok(Self {
            matrix,
            inv,
            lattices,
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
    pub fn shortest_vector<S>(&self, vec: &nalgebra::Vector<Float, Const<3>, S>) -> Vector3f
    where
        S: Storage<Float, Const<3>>,
    {
        self.shortest_vector_dims(vec, PBC_FULL)
    }

    /// Computes the shortest vector between two points considering periodicity only in specified dimensions.
    #[inline(always)]
    pub fn shortest_vector_dims<S>(
        &self,
        vec: &nalgebra::Vector<Float, Const<3>, S>,
        pbc_dims: PbcDims,
    ) -> Vector3f
    where
        S: Storage<Float, Const<3>>,
    {
        let vec = vec.clone_owned();
        if !pbc_dims.any() || !vec.iter().all(|v| v.is_finite()) {
            return vec;
        }
        let Some(lattices) = &self.lattices else {
            let mut result = vec;
            for i in 0..3 {
                if pbc_dims.get_dim(i) {
                    result[i] -= (vec[i] * self.inv[(i, i)]).round() * self.matrix[(i, i)];
                }
            }
            return result;
        };
        let index = match pbc_dims.0 {
            0b011 => 0,
            0b101 => 1,
            0b110 => 2,
            0b111 => 3,
            // One periodic vector: orthogonal projection gives the exact shift.
            _ => {
                let axis = pbc_dims.0.trailing_zeros() as usize;
                let v = self.matrix.column(axis);
                return vec - v * (vec.dot(&v) / v.norm_squared()).round();
            }
        };
        lattices[index].shortest(&vec)
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
    pub fn to_box_coords<S>(&self, vec: &nalgebra::Vector<Float, Const<3>, S>) -> Vector3f
    where
        S: Storage<Float, Const<3>>,
    {
        self.inv * vec
    }

    /// Checks if a point is inside the box (coordinates between 0 and 1 in box frame).
    #[inline(always)]
    pub fn is_inside(&self, point: &Pos) -> bool {
        let v = self.inv * point;
        v[0] < 1.0 && v[1] < 1.0 && v[2] < 1.0 && v[0] >= 0.0 && v[1] >= 0.0 && v[2] >= 0.0
    }

    /// Converts coordinates from box frame to lab frame.
    #[inline(always)]
    pub fn to_lab_coords<S>(&self, vec: &nalgebra::Vector<Float, Const<3>, S>) -> Vector3f
    where
        S: Storage<Float, Const<3>>,
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
        Vector3f::from_fn(|i, _| self.matrix.row(i).abs().sum())
    }

    /// Perpendicular distances between opposite unit-cell faces.
    /// A fractional search grid must use these spacings, not vector lengths.
    pub(crate) fn face_spacings(&self) -> Vector3f {
        Vector3f::from_fn(|i, _| self.inv.row(i).norm().recip())
    }

    /// Computes squared distance between two points considering periodic boundary conditions.
    #[inline(always)]
    pub fn distance_squared(&self, p1: &Pos, p2: &Pos, pbc_dims: PbcDims) -> Float {
        self.shortest_vector_dims(&(p2 - p1), pbc_dims)
            .norm_squared()
    }

    /// Computes distance between two points considering periodic boundary conditions.
    #[inline(always)]
    pub fn distance(&self, p1: &Pos, p2: &Pos, pbc: PbcDims) -> Float {
        self.distance_squared(p1, p2, pbc).sqrt()
    }

    /// Checks if the box is triclinic (has non-orthogonal vectors).
    pub fn is_triclinic(&self) -> bool {
        self.matrix[(0, 1)] != 0.0
            || self.matrix[(0, 2)] != 0.0
            || self.matrix[(1, 0)] != 0.0
            || self.matrix[(1, 2)] != 0.0
            || self.matrix[(2, 0)] != 0.0
            || self.matrix[(2, 1)] != 0.0
    }

    /// Scales box vectors by given factors.
    pub(crate) fn scale_vectors(
        &mut self,
        scale_factors: [Float; 3],
    ) -> Result<(), PeriodicBoxError> {
        let mut matrix = self.matrix;
        for c in 0..3 {
            matrix.column_mut(c).scale_mut(scale_factors[c]);
        }
        // Rebuild all derived data, and leave self unchanged if scaling fails.
        *self = Self::from_matrix(matrix)?;
        Ok(())
    }

    /// Wraps a point into the primary unit cell.
    #[inline(always)]
    pub fn wrap_point(&self, p: &Pos) -> Pos {
        Pos::from(self.wrap_vec(&p.coords))
    }

    #[inline(always)]
    pub fn wrap_vec<S>(&self, vec: &nalgebra::Vector<Float, Const<3>, S>) -> Vector3f
    where
        S: Storage<Float, Const<3>>,
    {
        // Get vector in box fractional coordinates
        let mut bv = self.inv * vec;
        for i in 0..3 {
            bv[i] -= bv[i].floor();
        }
        return self.matrix * bv;
    }
}

#[cfg(test)]
mod tests {
    use super::PeriodicBox;
    use crate::prelude::*;

    const EPSILON: Float = 1e-6;

    fn assert_vec_eq(v1: &Vector3f, v2: &Vector3f) {
        assert!(
            (v1 - v2).norm() < EPSILON,
            "Vectors not equal: {:?} != {:?}",
            v1,
            v2
        );
    }

    #[test]
    #[should_panic]
    fn invalid_from_vec_ang() {
        let _b = PeriodicBox::from_vectors_angles(10.0, 0.2, 15.0, 90.0, 9.0, 90.0).unwrap();
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
    fn test_orthogonal_needs_no_image_search() {
        let box_matrix = Matrix3f::from_diagonal(&Vector3f::new(10.0, 20.0, 30.0));
        let pbox = PeriodicBox::from_matrix(box_matrix).unwrap();
        assert!(pbox.lattices.is_none());
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
        let box_matrix = Matrix3f::new(10.0, 4.0, -4.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0);
        let pbox = PeriodicBox::from_matrix(box_matrix).unwrap();
        let p1 = Pos::new(38.9214, 40.0078, -34.0795);
        let p2 = Pos::new(-26.6187, 40.8926, 30.9709);
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
        let box_matrix = Matrix3f::new(6.0, 0.0, 3.0, 0.0, 6.0, 3.0, 0.0, 0.0, 6.0);
        let pbox = PeriodicBox::from_matrix(box_matrix).unwrap();
        let dx = Vector3f::new(2.9, 2.9, 2.9);

        // Reference: exhaustive 5x5x5 brute-force minimum over lattice images.
        let a: Vector3f = box_matrix.column(0).into();
        let b: Vector3f = box_matrix.column(1).into();
        let c: Vector3f = box_matrix.column(2).into();
        let mut best = Float::INFINITY;
        for i in -2..=2i32 {
            for j in -2..=2i32 {
                for k in -2..=2i32 {
                    let cand = dx + (i as Float) * a + (j as Float) * b + (k as Float) * c;
                    best = best.min(cand.norm());
                }
            }
        }

        let got = pbox.shortest_vector(&dx).norm();
        assert!((got - best).abs() < 1e-5, "expected {best}, got {got}");
    }

    // Far-apart points (|dx| >> box) must still be handled correctly via the
    // initial fractional-coord reduction step.
    #[test]
    fn test_triclinic_far_apart_reduction() {
        let box_matrix = Matrix3f::new(10.0, 4.0, -4.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0);
        let pbox = PeriodicBox::from_matrix(box_matrix).unwrap();
        // p1 and p2 ~6 box-lengths apart.
        let p1 = Pos::new(0.1, 0.2, 0.3);
        let p2 = Pos::new(60.1, 0.2, 0.3);
        let d = pbox.distance(&p1, &p2, PBC_FULL);
        assert!(
            d < 1e-4,
            "pure a-vector shift should collapse to ~0, got {d}"
        );
    }

    fn exhaustive(m: &Matrix3f, v: &Vector3f, pbc: PbcDims) -> Float {
        let range = |d| if pbc.get_dim(d) { -8..=8 } else { 0..=0 };
        let mut best = f64::INFINITY;
        for i in range(0) {
            for j in range(1) {
                for k in range(2) {
                    let n = [i as f64, j as f64, k as f64];
                    let mut d2 = 0.0;
                    for d in 0..3 {
                        let x = v[d] as f64 - (0..3).map(|a| m[(d, a)] as f64 * n[a]).sum::<f64>();
                        d2 += x * x;
                    }
                    best = best.min(d2);
                }
            }
        }
        best as Float
    }

    #[test]
    fn triclinic_requires_two_vector_shift() {
        let m = Matrix3f::new(10., 0., -2., 0., 10., 0., 0., 0., 1.);
        let b = PeriodicBox::from_matrix(m).unwrap();
        assert_vec_eq(
            &b.shortest_vector(&Vector3f::new(4.9, 0., 0.)),
            &Vector3f::new(0.9, 0., 2.),
        );
    }

    #[test]
    fn image_search_matches_exhaustive_all_masks() {
        let matrices = [
            Matrix3f::new(10., 4., -4., 0., 10., 2., 0., 0., 10.),
            Matrix3f::new(10., 0., -2., 0., 10., 0., 0., 0., 1.),
            // A non-triangular, rotated box exercises the QR path.
            Matrix3f::new(0., -10., -2., 10., 4., -4., 0., 0., 10.),
        ];
        let mut seed = 2026_u64;
        for m in matrices {
            let b = PeriodicBox::from_matrix(m).unwrap();
            for mask in 0..8 {
                let pbc = PbcDims::new(mask & 1 != 0, mask & 2 != 0, mask & 4 != 0);
                for _ in 0..32 {
                    let f = Vector3f::from_fn(|_, _| {
                        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                        ((seed >> 32) as u32 as Float / u32::MAX as Float - 0.5) * 3.0
                    });
                    let v = m * f;
                    let got = b.shortest_vector_dims(&v, pbc);
                    let expected = exhaustive(&m, &v, pbc);
                    assert!(
                        (got.norm_squared() - expected).abs() < 3e-5 * (1.0 + expected),
                        "mask {mask}, matrix {m:?}, displacement {v:?}: {got:?}, expected squared {expected}"
                    );
                    let shift = b.to_box_coords(&(v - got));
                    for d in 0..3 {
                        let expected = if pbc.get_dim(d) {
                            shift[d].round()
                        } else {
                            0.0
                        };
                        assert!((shift[d] - expected).abs() < 2e-5);
                    }
                }
            }
        }
    }

    #[test]
    fn partial_xy_uses_shortest_permitted_image() {
        let m = Matrix3f::new(10., 4., 0., 0., 10., 0., 0., 0., 10.);
        let b = PeriodicBox::from_matrix(m).unwrap();
        let v = m * Vector3f::new(0.49, 0.49, 0.);
        let got = b.shortest_vector_dims(&v, PbcDims::new(true, true, false));
        assert!((got - (v - m.column(0))).norm() < 2e-6);
        // The nonperiodic distance must not enlarge the lattice search sphere.
        let v = v + Vector3f::new(0., 0., 1e6);
        let got = b.shortest_vector_dims(&v, PbcDims::new(true, true, false));
        assert!((got - (v - m.column(0))).norm() < 2e-6);
    }

    #[test]
    fn one_periodic_vector_uses_projection() {
        let b =
            PeriodicBox::from_matrix(Matrix3f::new(10., 4., 0., 0., 10., 0., 0., 0., 10.)).unwrap();
        let v = Vector3f::new(100., 0., 0.);
        // Fractional y is zero, but the nearest image requires three b shifts.
        assert_vec_eq(
            &b.shortest_vector_dims(&v, PbcDims::new(false, true, false)),
            &Vector3f::new(88., -30., 0.),
        );
    }

    #[test]
    fn scaling_rebuilds_search_and_failure_is_atomic() {
        let m = Matrix3f::new(10., 4., 0., 0., 10., 0., 0., 0., 10.);
        let mut b = PeriodicBox::from_matrix(m).unwrap();
        b.scale_vectors([2., 2., 2.]).unwrap();
        let fresh = PeriodicBox::from_matrix(m * 2.).unwrap();
        let v = Vector3f::new(9., 9., 0.);
        assert_vec_eq(&b.shortest_vector(&v), &fresh.shortest_vector(&v));
        assert!(b.scale_vectors([1., 0., 1.]).is_err());
        assert_eq!(b.get_matrix(), fresh.get_matrix());
        assert_vec_eq(&b.shortest_vector(&v), &fresh.shortest_vector(&v));
    }

    #[test]
    fn wrapping_negative_and_multicell_coordinates() {
        let m = Matrix3f::new(10., 4., -4., 0., 10., 0., 0., 0., 10.);
        let b = PeriodicBox::from_matrix(m).unwrap();
        let v = m * Vector3f::new(-0.2, -2.7, 3.4);
        let want = m * Vector3f::new(0.8, 0.3, 0.4);
        assert!((b.wrap_vec(&v) - want).norm() < 1e-5);
        let p = b.wrap_point(&Pos::from(v));
        assert!(b.is_inside(&p));
        assert!((p.coords - want).norm() < 1e-5);
    }

    #[test]
    fn triclinic_face_spacings_and_signed_extents() {
        let m = Matrix3f::new(10., 4., -4., 0., 10., 0., 0., 0., 10.);
        let b = PeriodicBox::from_matrix(m).unwrap();
        assert_vec_eq(&b.get_lab_extents(), &Vector3f::new(18., 10., 10.));
        assert!((b.face_spacings()[0] - 10.0 / (1.32 as Float).sqrt()).abs() < 1e-6);
    }

    #[test]
    fn bounded_fallback_matches_cached_search_at_faces() {
        let m = Matrix3f::new(10., 4., -4., 0., 10., 2., 0., 0., 10.);
        let cached = PeriodicBox::from_matrix(m).unwrap();
        let mut bounded = cached.clone();
        for lattice in bounded.lattices.as_mut().unwrap().iter_mut() {
            lattice.corrections = None;
        }
        for mask in 0..8 {
            let pbc = PbcDims::new(mask & 1 != 0, mask & 2 != 0, mask & 4 != 0);
            for x in [-0.50001, -0.5, -0.49999, 0.49999, 0.5, 0.50001] {
                for y in [-0.5, 0.2, 0.5] {
                    let v = m * Vector3f::new(x, y, 0.5);
                    let expected = exhaustive(&m, &v, pbc);
                    for b in [&cached, &bounded] {
                        let got = b.shortest_vector_dims(&v, pbc).norm_squared();
                        assert!((got - expected).abs() < 3e-5 * (1.0 + expected));
                    }
                }
            }
        }
    }
}
