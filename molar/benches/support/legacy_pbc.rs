// Audit baseline: the previous fractional reduction and 26-shift search.
// It is NOT an exact oracle for thin boxes or partial PBC.
use molar::prelude::*;

pub struct LegacyBox {
    matrix: Matrix3f,
    inv: Matrix3f,
    shifts: Vec<Vector3f>,
    pairs: Vec<(Vector3f, Float)>,
}
impl LegacyBox {
    pub fn new(matrix: Matrix3f) -> Self {
        Self::build(matrix, true)
    }

    pub fn original(matrix: Matrix3f) -> Self {
        Self::build(matrix, false)
    }

    fn build(matrix: Matrix3f, paired: bool) -> Self {
        let inv = matrix.try_inverse().unwrap();
        let a = matrix.column(0).into_owned();
        let b = matrix.column(1).into_owned();
        let c = matrix.column(2).into_owned();
        let half_diag = 0.5
            * (a + b + c)
                .norm()
                .max((a + b - c).norm())
                .max((a - b + c).norm())
                .max((-a + b + c).norm());
        let bound2 = (2.0 * half_diag).powi(2);
        let mut shifts = Vec::new();
        if (0..3).any(|i| (0..3).any(|j| i != j && matrix[(i, j)] != 0.0)) {
            for i in -1..=1 {
                for j in -1..=1 {
                    for k in -1..=1 {
                        let s = a * i as Float + b * j as Float + c * k as Float;
                        if (i != 0 || j != 0 || k != 0) && s.norm_squared() < bound2 {
                            shifts.push(s);
                        }
                    }
                }
            }
        }
        let pairs = shifts
            .iter()
            .take(if paired { shifts.len() } else { 0 })
            .filter(|s| s.iter().find(|x| **x != 0.0).is_some_and(|x| *x > 0.0))
            .filter(|s| s.norm_squared() < (matrix.transpose() * (*s)).abs().sum())
            .map(|s| (*s, s.norm_squared()))
            .collect();
        Self {
            matrix,
            inv,
            shifts,
            pairs,
        }
    }
    #[inline]
    pub fn shortest(&self, v: &Vector3f, pbc: PbcDims) -> Vector3f {
        let mut f = self.inv * v;
        for i in 0..3 {
            if pbc.get_dim(i) {
                f[i] -= f[i].round();
            }
        }
        let start = self.matrix * f;
        if pbc != PBC_FULL {
            return start;
        }
        let mut best = start;
        let mut best2 = start.norm_squared();
        for s in &self.shifts {
            let cand = start + s;
            let n2 = cand.norm_squared();
            if n2 < best2 {
                best2 = n2;
                best = cand;
            }
        }
        best
    }
    #[inline]
    pub fn paired(&self, v: &Vector3f) -> Vector3f {
        let mut f = self.inv * v;
        f.apply(|x| *x -= x.round());
        let start = self.matrix * f;
        let start2 = start.norm_squared();
        let mut best = start;
        let mut best2 = start2;
        for (s, s2) in &self.pairs {
            let dot = start.dot(s);
            let n2 = start2 + s2 - 2.0 * dot.abs();
            if n2 < best2 {
                best2 = n2;
                best = start - s * dot.signum();
            }
        }
        best
    }
}
