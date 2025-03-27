use molar::prelude::*;
use nalgebra::DVector;

use crate::lipid_species::LipidSpecies;

pub struct LipidProperties {
    pub(super) normal: Vector3f,
    pub(super) area: f32,
    pub(super) tilt: f32,
    pub(super) order: Vec<DVector<f32>>,
    // mean_curv: f32,
    // gauss_curv: f32,
    // princ_curv: [f32; 2],
    // princ_curv_axes: [Vector3f; 2],
}