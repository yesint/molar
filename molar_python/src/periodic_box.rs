use crate::utils::*;
use molar::prelude::*;
use numpy::{
    nalgebra::{Const, Dyn, MatrixView, VectorView},
    AllowTypeChange, PyArray1, PyArray2, PyArrayLike1, PyArrayLike2, ToPyArray,
};
use pyo3::{
    exceptions::{PyTypeError, PyValueError},
    prelude::*,
    types::PyTuple,
};
/// Periodic simulation box geometry and minimum-image/PBC utilities.
///
/// Construct either from a `3x3` box matrix or from vectors and angles.
///
/// **Example**
///
/// .. code-block:: python
///
///    import pymolar as molar
///    import numpy as np
///
///    box = molar.PeriodicBox(np.eye(3, dtype=np.float32) * 10.0)
///    d = box.distance(
///        np.array([0.0, 0.0, 0.0], dtype=np.float32),
///        np.array([9.0, 0.0, 0.0], dtype=np.float32),
///        [True, True, True],
///    )

#[pyclass(name = "PeriodicBox")]
pub(super) struct PeriodicBoxPy(pub(super) PeriodicBox);

#[pymethods]
impl PeriodicBoxPy {
    #[new]
    #[pyo3(signature = (*py_args))]
    /// Create a periodic box.
    ///
    /// :param py_args:
    ///     Positional arguments:
    ///
    ///     - ``(matrix,)`` where ``matrix`` is a ``3x3`` array-like
    ///     - ``(vectors, angles)`` where each is length-3 array-like
    /// :returns: Initialized periodic box.
    /// :rtype: PeriodicBox
    /// :raises TypeError: If argument count is not 1 or 2.
    /// :raises ValueError: If input shapes are incompatible.
    fn new<'py>(py_args: &Bound<'py, PyTuple>) -> PyResult<Self> {
        match py_args.len() {
            1 => {
                // From matrix
                let arr: PyArrayLike2<'py, f32, AllowTypeChange> =
                    py_args.get_item(0)?.extract()?;

                let m: MatrixView<f32, Const<3>, Const<3>, Dyn, Dyn> = arr
                    .try_as_matrix()
                    .ok_or_else(|| PyValueError::new_err("conversion to 3x3 matrix has failed"))?;

                let pb = PeriodicBox::from_matrix(m)
                    .map_err(|e| PyValueError::new_err(e.to_string()))?;

                Ok(PeriodicBoxPy(pb))
            }
            2 => {
                // From vectors and angles
                let v_arr: PyArrayLike1<'py, f32, AllowTypeChange> =
                    py_args.get_item(0)?.extract()?;
                let a_arr: PyArrayLike1<'py, f32, AllowTypeChange> =
                    py_args.get_item(1)?.extract()?;

                let v: VectorView<f32, Const<3>, Dyn> = v_arr.try_as_matrix().ok_or_else(|| {
                    PyValueError::new_err("conversion of vectors to Vector3 has failed")
                })?;

                let a: VectorView<f32, Const<3>, Dyn> = a_arr.try_as_matrix().ok_or_else(|| {
                    PyValueError::new_err("conversion of angles to Vector3 has failed")
                })?;

                let pb = PeriodicBox::from_vectors_angles(v[0], v[1], v[2], a[0], a[1], a[2])
                    .map_err(|e| PyValueError::new_err(e.to_string()))?;

                Ok(PeriodicBoxPy(pb))
            }
            _ => Err(PyTypeError::new_err(
                "wrong number of arguments: 1 or 2 required",
            )),
        }
    }

    /// Return vector-length and angle representation.
    ///
    /// :returns: Tuple ``(vectors, angles)`` as NumPy arrays of length 3.
    /// :rtype: tuple[numpy.ndarray, numpy.ndarray]
    fn to_vectors_angles<'py>(
        slf: Bound<'py, Self>,
    ) -> (Bound<'py, PyArray1<f32>>, Bound<'py, PyArray1<f32>>) {
        let (v, a) = slf.borrow().0.to_vectors_angles();
        let v_arr = clone_vec_to_pyarray1(&v, slf.py());
        let a_arr = clone_vec_to_pyarray1(&a, slf.py());
        (v_arr, a_arr)
    }

    #[pyo3(signature = (arr, dims=None))]
    #[pyo3(text_signature = "($self, arr, dims=None)")]
    /// Compute minimum-image displacement vector.
    ///
    /// :param arr: Input displacement vector (length 3).
    /// :param dims: Periodic dimensions as ``[x, y, z]`` booleans.
    /// :returns: Minimum-image displacement vector.
    /// :rtype: numpy.ndarray
    /// :raises ValueError: If vector conversion fails.
    fn shortest_vector<'py>(
        &self,
        py: Python<'py>,
        arr: PyArrayLike1<'py, f32, AllowTypeChange>,
        dims: Option<[bool; 3]>,
    ) -> PyResult<Bound<'py, PyArray1<f32>>> {
        let dims = dims.unwrap_or([true, true, true]);
        let v: VectorView<f32, Const<3>, Dyn> = arr
            .try_as_matrix()
            .ok_or_else(|| PyValueError::new_err("conversion to Vector3 has failed"))?;

        let pbc = PbcDims::new(dims[0], dims[1], dims[2]);
        let out_v = self.0.shortest_vector_dims(&v, pbc);
        Ok(clone_vec_to_pyarray1(&out_v, py))
    }

    #[pyo3(signature = (point, target, dims=None))]
    #[pyo3(text_signature = "($self, point, target, dims=None)")]
    /// Return periodic image of ``point`` closest to ``target``.
    ///
    /// :param point: Point to image (length 3).
    /// :param target: Reference target point (length 3).
    /// :param dims: Periodic dimensions as ``[x, y, z]`` booleans.
    /// :returns: Imaged point closest to ``target``.
    /// :rtype: numpy.ndarray
    /// :raises ValueError: If vector conversion fails.
    fn closest_image<'py>(
        &self,
        py: Python<'py>,
        point: PyArrayLike1<'py, f32, AllowTypeChange>,
        target: PyArrayLike1<'py, f32, AllowTypeChange>,
        dims: Option<[bool; 3]>,
    ) -> PyResult<Bound<'py, PyArray1<f32>>> {
        let dims = dims.unwrap_or([true, true, true]);
        let p: VectorView<f32, Const<3>, Dyn> = point
            .try_as_matrix()
            .ok_or_else(|| PyValueError::new_err("conversion of point to Vector3 has failed"))?;
        let t: VectorView<f32, Const<3>, Dyn> = target
            .try_as_matrix()
            .ok_or_else(|| PyValueError::new_err("conversion of target to Vector3 has failed"))?;

        let pbc = PbcDims::new(dims[0], dims[1], dims[2]);
        let out = t + self.0.shortest_vector_dims(&(p - t), pbc);
        Ok(clone_vec_to_pyarray1(&out, py))
    }

    /// Return full box matrix.
    ///
    /// :returns: ``3x3`` box matrix.
    /// :rtype: numpy.ndarray
    fn get_matrix<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f32>> {
        self.0.get_matrix().to_pyarray(py)
    }

    /// Convert Cartesian coordinates to box coordinates.
    ///
    /// :param point: Cartesian coordinate vector (length 3).
    /// :returns: Coordinate in box basis.
    /// :rtype: numpy.ndarray
    /// :raises ValueError: If vector conversion fails.
    fn to_box_coords<'py>(
        &self,
        py: Python<'py>,
        point: PyArrayLike1<'py, f32, AllowTypeChange>,
    ) -> PyResult<Bound<'py, PyArray1<f32>>> {
        let p: VectorView<f32, Const<3>, Dyn> = point
            .try_as_matrix()
            .ok_or_else(|| PyValueError::new_err("conversion of point to Vector3 has failed"))?;
        let v = self.0.to_box_coords(&p);
        Ok(clone_vec_to_pyarray1(&v, py))
    }

    /// Convert box coordinates to Cartesian coordinates.
    ///
    /// :param point: Coordinate in box basis (length 3).
    /// :returns: Coordinate in Cartesian/lab basis.
    /// :rtype: numpy.ndarray
    /// :raises ValueError: If vector conversion fails.
    fn to_lab_coords<'py>(
        &self,
        py: Python<'py>,
        point: PyArrayLike1<'py, f32, AllowTypeChange>,
    ) -> PyResult<Bound<'py, PyArray1<f32>>> {
        let p: VectorView<f32, Const<3>, Dyn> = point
            .try_as_matrix()
            .ok_or_else(|| PyValueError::new_err("conversion of point to Vector3 has failed"))?;
        let v = self.0.to_lab_coords(&p);
        Ok(clone_vec_to_pyarray1(&v, py))
    }

    /// Return box extents in box basis.
    ///
    /// :returns: Length-3 vector of extents.
    /// :rtype: numpy.ndarray
    fn get_box_extents<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f32>> {
        clone_vec_to_pyarray1(&self.0.get_box_extents(), py)
    }

    /// Return box extents in lab basis.
    ///
    /// :returns: Length-3 vector of extents in lab basis.
    /// :rtype: numpy.ndarray
    fn get_lab_extents<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f32>> {
        clone_vec_to_pyarray1(&self.0.get_box_extents(), py)
    }

    /// Check whether the box is triclinic.
    ///
    /// :returns: ``True`` for triclinic geometry, else ``False``.
    /// :rtype: bool
    fn is_triclinic(&self) -> bool {
        self.0.is_triclinic()
    }

    /// Compute squared distance with optional periodic dimensions.
    ///
    /// :param p1: First point (length 3).
    /// :param p2: Second point (length 3).
    /// :param dims: Periodic dimensions as ``[x, y, z]`` booleans.
    /// :returns: Squared minimum-image distance.
    /// :rtype: float
    /// :raises ValueError: If vector conversion fails.
    fn distance_squared<'py>(
        &self,
        p1: PyArrayLike1<'py, f32, AllowTypeChange>,
        p2: PyArrayLike1<'py, f32, AllowTypeChange>,
        dims: [bool; 3],
    ) -> PyResult<f32> {
        let p1: VectorView<f32, Const<3>, Dyn> = p1
            .try_as_matrix()
            .ok_or_else(|| PyValueError::new_err("conversion of point to Vector3 has failed"))?;
        let p2: VectorView<f32, Const<3>, Dyn> = p2
            .try_as_matrix()
            .ok_or_else(|| PyValueError::new_err("conversion of target to Vector3 has failed"))?;

        let pbc = PbcDims::new(dims[0], dims[1], dims[2]);
        Ok(self.0.shortest_vector_dims(&(p2 - p1), pbc).norm_squared())
    }

    /// Compute distance with optional periodic dimensions.
    ///
    /// :param p1: First point (length 3).
    /// :param p2: Second point (length 3).
    /// :param dims: Periodic dimensions as ``[x, y, z]`` booleans.
    /// :returns: Minimum-image distance.
    /// :rtype: float
    fn distance<'py>(
        &self,
        p1: PyArrayLike1<'py, f32, AllowTypeChange>,
        p2: PyArrayLike1<'py, f32, AllowTypeChange>,
        dims: [bool; 3],
    ) -> PyResult<f32> {
        Ok(self.distance_squared(p1, p2, dims)?.sqrt())
    }

    /// Wrap point into the primary unit cell.
    ///
    /// :param p: Input point (length 3).
    /// :returns: Wrapped point.
    /// :rtype: numpy.ndarray
    /// :raises ValueError: If vector conversion fails.
    fn wrap_point<'py>(
        &self,
        py: Python<'py>,
        p: PyArrayLike1<'py, f32, AllowTypeChange>,
    ) -> PyResult<Bound<'py, PyArray1<f32>>> {
        let v: VectorView<f32, Const<3>, Dyn> = p
            .try_as_matrix()
            .ok_or_else(|| PyValueError::new_err("conversion of point to Vector3 has failed"))?;
        Ok(clone_vec_to_pyarray1(&self.0.wrap_vec(&v), py))
    }
}
