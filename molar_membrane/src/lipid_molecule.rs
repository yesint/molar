use molar::prelude::*;
use nalgebra::DVector;
use std::rc::Rc;

use crate::lipid_species::LipidSpecies;

pub struct LipidMolecule {
    pub id: usize,
    pub valid: bool,
    pub sel: Sel<ImmutableParallel>,
    pub species: Rc<LipidSpecies>,
    pub head_sel: Sel<ImmutableParallel>,
    pub mid_sel: Sel<ImmutableParallel>,
    pub tail_end_sel: Sel<ImmutableParallel>,
    pub tail_sels: Vec<Sel<ImmutableParallel>>,
    pub head_marker: Pos,
    pub mid_marker: Pos,
    pub tail_marker: Pos,

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

    pub fn set_state(&mut self, st: State) -> anyhow::Result<()> {
        let st: Holder<State,ImmutableParallel> = st.into();
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
