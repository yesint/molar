use molar::prelude::*;
use nalgebra::DVector;
use std::sync::Arc;

use crate::lipid_species::LipidSpecies;

#[derive(Debug)]
pub struct SingleLipidProperties {
    pub(super) normal: Vector3f,
    pub(super) area: f32,
    pub(super) tilt: f32,
    pub(super) order: Vec<DVector<f32>>,
}

impl SingleLipidProperties {
    pub fn new(species: &LipidSpecies) -> Self {
        let mut order = Vec::with_capacity(species.tails.len());
        for t in &species.tails {
            order.push(DVector::from_element(t.bond_orders.len() - 1, 0.0));
        }
        Self {
            normal: Default::default(),
            area: Default::default(),
            tilt: Default::default(),
            order,
        }
    }
}

pub struct LipidMolecule {
    pub(super) sel: Sel<MutableSerial>,
    pub(super) species: Arc<LipidSpecies>,
    pub(super) head_sel: Sel<MutableSerial>,
    pub(super) mid_sel: Sel<MutableSerial>,
    pub(super) tail_end_sel: Sel<MutableSerial>,
    pub(super) tail_sels: Vec<Sel<MutableSerial>>,
    pub(super) head_marker: Pos,
    pub(super) mid_marker: Pos,
    pub(super) tail_marker: Pos,
    pub(super) stats: SingleLipidProperties,
}

impl LipidMolecule {
    pub fn update_markers(&mut self) -> anyhow::Result<()> {
        self.head_marker = self.head_sel.center_of_mass_pbc()?;
        self.mid_marker = self.mid_sel.center_of_mass_pbc()?;
        self.tail_marker = self.tail_end_sel.center_of_mass_pbc()?;
        Ok(())
    }

    pub fn compute_order(&mut self, order_type: OrderType, normal: &Vector3f) {
        for i in 0..self.tail_sels.len() {
            self.stats.order[i] = self.tail_sels[i]
                .lipid_tail_order(
                    order_type.clone(),
                    &vec![normal.clone()],
                    &self.species.tails[i].bond_orders,
                )
                .unwrap();
        }
    }

    pub fn num_tails(&self) -> usize {
        self.tail_sels.len()
    }
}
