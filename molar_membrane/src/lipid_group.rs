use molar::core::{MutableKind, UserCreatableKind};

use crate::{stats::GroupProperties, LipidMolecule, Surface};

#[derive(Default)]
pub struct LipidGroup {
    pub(crate) lipid_ids: Vec<usize>,
    pub(crate) stats: GroupProperties,
}

impl LipidGroup {
    pub(crate) fn frame_update(
        &mut self,
        lipids: &Vec<LipidMolecule>,
        surf: &Surface,
    ) -> anyhow::Result<()> {
        // Init update
        for stat in self.stats.per_species.values_mut() {
            stat.init_frame_update();
        }
        // Update group by adding individual lipid stats
        for lip_id in &self.lipid_ids {
            // If lipid is not valid skip
            if !lipids[*lip_id].valid {continue}
            // Get species name
            let sp_name = &lipids[*lip_id].species.name;
            // Update stats for this species
            self.stats
                .per_species
                .get_mut(sp_name)
                .unwrap()
                .add_single_lipid_stats(*lip_id, &lipids, &surf)?;
        }
        // Finish update
        for stat in self.stats.per_species.values_mut() {
            stat.finish_frame_update();
        }
        Ok(())
    }
}
