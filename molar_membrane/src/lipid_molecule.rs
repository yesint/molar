use molar::prelude::*;
use nalgebra::DVector;
use std::sync::Arc;

use crate::lipid_species::LipidSpecies;

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

    //pub(super) props: SingleLipidProperties,
    pub(super) order: Vec<DVector<f32>>,
    pub(super) tail_head_vec: Vector3f,
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
            self.order[i] = self.tail_sels[i]
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

    pub fn set_state(&mut self, st: Holder<State, MutableSerial>) -> anyhow::Result<()> {
        self.sel.set_state(st.clone())?;
        self.head_sel.set_state(st.clone())?;
        self.mid_sel.set_state(st.clone())?;
        for t in &mut self.tail_sels {
            t.set_state(st.clone())?;
        }
        self.update_markers()?;
        Ok(())
    }
}
