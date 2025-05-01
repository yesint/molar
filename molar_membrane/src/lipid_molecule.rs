use molar::prelude::*;
use nalgebra::DVector;
use std::rc::Rc;

use crate::lipid_species::LipidSpecies;

pub struct LipidMolecule<K> {
    pub id: usize,
    pub sel: Sel<K>,
    pub species: Rc<LipidSpecies>,
    pub head_sel: Sel<K>,
    pub mid_sel: Sel<K>,
    pub tail_end_sel: Sel<K>,
    pub tail_sels: Vec<Sel<K>>,
    pub head_marker: Pos,
    pub mid_marker: Pos,
    pub tail_marker: Pos,

    //pub(super) props: SingleLipidProperties,
    pub(super) order: Vec<DVector<f32>>,
    pub(super) tail_head_vec: Vector3f,

    pub(super) valid: bool,
}

impl<K: UserCreatableKind+MutableKind> LipidMolecule<K> {
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

    pub fn set_state(&mut self, st: Holder<State, K>) -> anyhow::Result<()> {
        self.sel.set_state(st.new_ref())?;
        self.head_sel.set_state(st.new_ref())?;
        self.mid_sel.set_state(st.new_ref())?;
        for t in &mut self.tail_sels {
            t.set_state(st.new_ref())?;
        }
        self.update_markers()?;
        Ok(())
    }
}
