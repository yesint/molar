use crate::utils::*;
use anyhow::{anyhow, bail};
use molar::prelude::*;
use numpy::{
    nalgebra::{Const, Dyn, MatrixView, VectorView},
    AllowTypeChange, PyArray1, PyArray2, PyArrayLike1, PyArrayLike2, ToPyArray,
};
use pyo3::{prelude::*, types::PyTuple};

#[pyclass]
pub(super) struct PeriodicBox(pub(super) molar::core::PeriodicBox);

#[pymethods]
impl PeriodicBox {
    #[new]
    #[pyo3(signature = (*py_args))]
    fn new<'py>(py_args: &Bound<'py, PyTuple>) -> anyhow::Result<Self> {
        if py_args.len() == 1 {
            // From matrix
            let arr: PyArrayLike2<'py, f32, AllowTypeChange> = py_args.get_item(0)?.extract()?;
            let m: MatrixView<f32, Const<3>, Const<3>, Dyn, Dyn> = arr
                .try_as_matrix()
                .ok_or_else(|| anyhow!("conversion to 3x3 matrix has failed"))?;
            Ok(PeriodicBox(molar::core::PeriodicBox::from_matrix(m)?))
        } else if py_args.len() == 2 {
            // From vectors and angles
            let v_arr: PyArrayLike1<'py, f32, AllowTypeChange> = py_args.get_item(0)?.extract()?;
            let a_arr: PyArrayLike1<'py, f32, AllowTypeChange> = py_args.get_item(1)?.extract()?;
            let v: VectorView<f32, Const<3>, Dyn> = v_arr
                .try_as_matrix()
                .ok_or_else(|| anyhow!("conversion of vectors to Vector3 has failed"))?;
            let a: VectorView<f32, Const<3>, Dyn> = a_arr
                .try_as_matrix()
                .ok_or_else(|| anyhow!("conversion of angles to Vector3 has failed"))?;
            Ok(PeriodicBox(molar::core::PeriodicBox::from_vectors_angles(
                v[0], v[1], v[2], a[0], a[1], a[2],
            )?))
        } else {
            bail!("wrong number of arguments: 1 or 2 reqired")
        }
    }

    fn to_vectors_angles<'py>(
        slf: Bound<'py, Self>,
    ) -> (Bound<'py, PyArray1<f32>>, Bound<'py, PyArray1<f32>>) {
        let (v, a) = slf.borrow().0.to_vectors_angles();
        let v_arr = clone_vec_to_pyarray1(&v, slf.py());
        let a_arr = clone_vec_to_pyarray1(&a, slf.py());
        (v_arr, a_arr)
    }

    #[pyo3(signature = (arr, dims=[true,true,true]))]
    fn shortest_vector<'py>(
        &self,
        py: Python<'py>,
        arr: PyArrayLike1<'py, f32, AllowTypeChange>,
        dims: [bool; 3],
    ) -> anyhow::Result<Bound<'py, PyArray1<f32>>> {
        let v: VectorView<f32, Const<3>, Dyn> = arr
            .try_as_matrix()
            .ok_or_else(|| anyhow!("conversion to Vector3 has failed"))?;

        let pbc = PbcDims::new(dims[0], dims[1], dims[2]);
        let out_v = self.0.shortest_vector_dims(&v, pbc);
        Ok(clone_vec_to_pyarray1(&out_v, py))
    }

    #[pyo3(signature = (point, target, dims=[true,true,true]))]
    fn closest_image<'py>(
        &self,
        py: Python<'py>,
        point: PyArrayLike1<'py, f32, AllowTypeChange>,
        target: PyArrayLike1<'py, f32, AllowTypeChange>,
        dims: [bool; 3],
    ) -> anyhow::Result<Bound<'py, PyArray1<f32>>> {
        let p: VectorView<f32, Const<3>, Dyn> = point
            .try_as_matrix()
            .ok_or_else(|| anyhow!("conversion of point to Vector3 has failed"))?;
        let t: VectorView<f32, Const<3>, Dyn> = target
            .try_as_matrix()
            .ok_or_else(|| anyhow!("conversion of target to Vector3 has failed"))?;

        let pbc = PbcDims::new(dims[0], dims[1], dims[2]);
        let out = t + self.0.shortest_vector_dims(&(p - t), pbc);
        Ok(clone_vec_to_pyarray1(&out, py))
    }

    fn get_matrix<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f32>> {
        self.0.get_matrix().to_pyarray(py)
    }

    fn to_box_coords<'py>(
        &self,
        py: Python<'py>,
        point: PyArrayLike1<'py, f32, AllowTypeChange>,
    ) -> anyhow::Result<Bound<'py, PyArray1<f32>>> {
        let p: VectorView<f32, Const<3>, Dyn> = point
            .try_as_matrix()
            .ok_or_else(|| anyhow!("conversion of point to Vector3 has failed"))?;
        let v = self.0.to_box_coords(&p);
        Ok(clone_vec_to_pyarray1(&v, py))
    }

    fn to_lab_coords<'py>(
        &self,
        py: Python<'py>,
        point: PyArrayLike1<'py, f32, AllowTypeChange>,
    ) -> anyhow::Result<Bound<'py, PyArray1<f32>>> {
        let p: VectorView<f32, Const<3>, Dyn> = point
            .try_as_matrix()
            .ok_or_else(|| anyhow!("conversion of point to Vector3 has failed"))?;
        let v = self.0.to_lab_coords(&p);
        Ok(clone_vec_to_pyarray1(&v, py))
    }

    fn get_box_extents<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f32>> {
        clone_vec_to_pyarray1(&self.0.get_box_extents(), py)
    }

    fn get_lab_extents<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f32>> {
        clone_vec_to_pyarray1(&self.0.get_box_extents(), py)
    }

    fn is_triclinic(&self) -> bool {
        self.0.is_triclinic()
    }

    fn distance_squared<'py>(
        &self,
        p1: PyArrayLike1<'py, f32, AllowTypeChange>,
        p2: PyArrayLike1<'py, f32, AllowTypeChange>,
        dims: [bool; 3],
    ) -> anyhow::Result<f32> {
        let p1: VectorView<f32, Const<3>, Dyn> = p1
            .try_as_matrix()
            .ok_or_else(|| anyhow!("conversion of point to Vector3 has failed"))?;
        let p2: VectorView<f32, Const<3>, Dyn> = p2
            .try_as_matrix()
            .ok_or_else(|| anyhow!("conversion of target to Vector3 has failed"))?;

        let pbc = PbcDims::new(dims[0], dims[1], dims[2]);
        Ok(self.0.shortest_vector_dims(&(p2 - p1), pbc).norm_squared())
    }

    fn distance<'py>(
        &self,
        p1: PyArrayLike1<'py, f32, AllowTypeChange>,
        p2: PyArrayLike1<'py, f32, AllowTypeChange>,
        dims: [bool; 3],
    ) -> anyhow::Result<f32> {
        Ok(self.distance_squared(p1, p2, dims)?.sqrt())
    }
}
