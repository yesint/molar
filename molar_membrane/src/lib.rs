use anyhow::{bail, Context};
use molar::prelude::*;
use std::{
    any, collections::HashMap, path::{Path, PathBuf}, sync::Arc
};

#[cfg(not(test))]
use log::{info, warn}; // Use log crate when building application

#[cfg(test)]
use std::{println as info, println as warn};

mod stats;
use stats::{GroupProperties, StatProperties};

mod lipid_molecule;
use lipid_molecule::{LipidMolecule, SingleLipidProperties};

mod lipid_species;
use lipid_species::{LipidSpecies, LipidSpeciesDescr};

pub struct Membrane {
    // source: Source<MutableSerial>,
    //species: Vec<(String, LipidSpecies)>,
    lipids: Vec<LipidMolecule>,
    groups: HashMap<String, LipidGroup>,
    global_normal: Option<Vector3f>,
    order_type: OrderType,
    output_dir: PathBuf,

    // Local patches
    patches: Vec<usize>,
}

impl Membrane {
    pub fn new(source: Source<MutableSerial>, defstr: &str) -> anyhow::Result<Self> {
        // Load species descriptions
        let species_descr: HashMap<String, LipidSpeciesDescr> = toml::from_str(defstr)?;

        // Load lipids from provided source
        //let mut species: HashMap<String, Arc<LipidSpecies>> = Default::default();
        let mut lipids = vec![];
        for (name, descr) in species_descr.iter() {
            if let Ok(lips) = source.select(&descr.whole) {
                let lips = lips.split_resindex::<Vec<_>>()?;
                // Use first lipid to create lipid species object
                info!("Creating {} '{}' lipids", lips.len(), name);
                let sp = Arc::new(LipidSpecies::new(name.clone(), descr.clone(), &lips[0])?);

                // Now create individual lipids
                for lip in lips {
                    //let sp = &species[name];
                    let head_sel = lip.subsel(&sp.head_marker_offsets)?;
                    let mid_sel = lip.subsel(&sp.mid_marker_offsets)?;

                    let mut tail_sels = vec![];
                    let mut tail_ends = vec![];
                    for t in sp.tails.iter() {
                        tail_sels.push(lip.subsel(&t.offsets)?);
                        tail_ends.push(*t.offsets.last().unwrap());
                    }
                    let tail_end_sel = lip.subsel(tail_ends)?;

                    // Unwrap the lipid
                    lip.unwrap_simple()?;
                    // Compute markers
                    let head_marker = head_sel.center_of_mass()?;
                    let mid_marker = mid_sel.center_of_mass()?;
                    let tail_marker = tail_end_sel.center_of_mass()?;

                    lipids.push(LipidMolecule {
                        sel: lip, // Selection is moved to the lipid
                        species: Arc::clone(&sp),
                        head_sel,
                        mid_sel,
                        tail_end_sel,
                        tail_sels,
                        head_marker,
                        mid_marker,
                        tail_marker,
                        stats: SingleLipidProperties::new(&sp),
                    });
                }
            } else {
                // No such lipid species found
                warn!("No '{name}' lipids found");
            }
        }
        Ok(Self {
            // source,
            // species,
            lipids,
            groups: Default::default(),
            global_normal: None,
            order_type: OrderType::ScdCorr,
            output_dir: PathBuf::from("membrane_results"),
            patches: vec![],
        })
    }

    pub fn with_global_normal(mut self, global_normal: impl Into<Option<Vector3f>>) -> Self {
        self.global_normal = global_normal.into();
        self
    }

    pub fn with_order_type(mut self, order_type: OrderType) -> Self {
        self.order_type = order_type;
        self
    }

    pub fn with_output_dir(mut self, output_dir: impl AsRef<Path>) -> anyhow::Result<Self> {
        // Check if we can access this directory or create it if needed
        let path = output_dir.as_ref();
        if !path.exists() {
            std::fs::create_dir_all(path).with_context(|| {
                format!("Failed to create output directory '{}'", path.display())
            })?;
        }

        self.output_dir = path.to_path_buf();
        info!("will use output directory '{}'", path.display());
        Ok(self)
    }

    pub fn iter_lipids(&self) -> impl Iterator<Item = (usize, &LipidMolecule)> {
        self.lipids.iter().enumerate()
    }

    pub fn add_group(&mut self, gr_name: impl AsRef<str>, ids: Vec<usize>) -> anyhow::Result<()> {
        let gr_name = gr_name.as_ref();
        for id in &ids {
            if *id >= self.lipids.len() {
                bail!(
                    "lipid id {} is out of bounds 0:{}",
                    id,
                    self.lipids.len() - 1
                );
            }
        }
        if self.groups.contains_key(gr_name) {
            bail!("group '{}' already exists", gr_name);
        } else {
            info!("Creating group '{}' with {} lipids", gr_name, ids.len());
            self.groups
                .insert(gr_name.to_owned(), LipidGroup::new(&self, ids));
        }
        Ok(())
    }

    pub fn compute(&mut self) -> anyhow::Result<()> {
        // Compute order for all lipids
        for lip in &mut self.lipids {
            lip.compute_order(self.order_type.clone(), &self.global_normal.unwrap());
        }
        // Add stats to groups
        for gr in self.groups.values_mut() {
            gr.add_lipid_stats(&self.lipids)?;
        }
        Ok(())
    }

    pub fn finalize(&self) -> anyhow::Result<()> {
        info!(
            "Writing results to directory '{}'",
            self.output_dir.display()
        );
        // Write results for groups
        for (gr_name, gr) in &self.groups {
            info!("\tGroup '{gr_name}'...");
            gr.stats
                .save_order_to_file(self.output_dir.as_path(), gr_name)?;
        }
        Ok(())
    }

    pub fn set_state(&mut self, st: Holder<State, MutableSerial>) -> anyhow::Result<()> {
        for lip in &mut self.lipids {
            lip.set_state(st.clone())?;
            lip.update_markers()?;
        }
        Ok(())
    }

    fn compute_patches(&mut self) -> anyhow::Result<()>{
        let ind: Vec<(usize,usize)> = distance_search_single_pbc(
            2.0,
            self.lipids.iter().map(|l| &l.mid_marker),
            0..self.lipids.len(),
            self.lipids[0].sel.require_box()?,
            PBC_FULL,
        );
        let conn = LocalConnectivity::from_iter(ind, self.lipids.len());
        println!("{:?}",conn);
        Ok(())
    }
}

pub struct LipidGroup {
    lipid_ids: Vec<usize>,
    stats: GroupProperties,
}

impl LipidGroup {
    fn new(memb: &Membrane, lipid_ids: Vec<usize>) -> Self {
        let mut stats = GroupProperties::default();
        for id in lipid_ids.iter() {
            let name = memb.lipids[*id].species.name.clone();

            if let Some(el) = stats.per_species.get_mut(&name) {
                el.num_lip += 1;
            } else {
                let sp = &memb.lipids[*id].species;
                stats.per_species.insert(name, StatProperties::new(sp));
            }
        }

        for (sp, stats) in &stats.per_species {
            info!("\t{}: {}", sp, stats.num_lip);
        }
        Self { lipid_ids, stats }
    }

    fn add_lipid_stats(&mut self, all_lipids: &Vec<LipidMolecule>) -> anyhow::Result<()> {
        self.stats
            .add_lipid_stats(self.lipid_ids.iter().map(|id| &all_lipids[*id]))?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::{io::Read, path::PathBuf};

    use molar::prelude::*;

    use crate::{lipid_species::PredefinedLipidSpecies, LipidSpecies, LipidSpeciesDescr, Membrane};

    #[test]
    fn test_descr() -> anyhow::Result<()> {
        let src = Source::serial_from_file("tests/membr.gro")?;
        let pope = src.select("resid 60")?;
        let descr = LipidSpeciesDescr {
            whole: "resname POPE".into(),
            head: "name P N".into(),
            mid: "name C21 C22".into(),
            tails: vec![
                "C21-C22-C23-C24-C25-C26-C27-C28-C29=C210-C211-C212-C213-C214-C215-C216-C217-C218"
                    .into(),
                "C31-C32-C33-C34-C35-C36-C37-C38-C39=C310-C311-C312-C313-C314-C315-C316".into(),
            ],
        };
        let lip_sp = LipidSpecies::new("DOPE".to_owned(), descr, &pope)?;
        println!("{lip_sp:?}");
        Ok(())
    }

    // #[test]
    // fn test_descr_from_itp() -> anyhow::Result<()> {
    //     let top = FileHandler::open("../../tests/POPE.itp")?.read_topology()?;
    //     let n = top.num_atoms();
    //     let src = Source::new_serial(top.into(), State::new_fake(n).into())?;
    //     let pope = src.select_all()?;
    //     let resname = &pope.first_atom().resname;
    //     let descr = LipidSpeciesDescr::pope();
    //     let lip_sp = LipidSpecies::new(descr, &pope)?;
    //     println!("{lip_sp:?}");
    //     Ok(())
    // }

    #[test]
    fn test_descr_serde() -> anyhow::Result<()> {
        let top = FileHandler::open("tests/POPE.itp")?.read_topology()?;
        let n = top.num_atoms();
        let src = Source::new_serial(top.into(), State::new_fake(n).into())?;

        let descr: LipidSpeciesDescr = toml::from_str(
            r#"
            whole = "resname POPE"
            head = "name P N"
            mid = "name C21 C22"
            tails = [
                "C21-C22-C23-C24-C25-C26-C27-C28-C29=C210-C211-C212-C213-C214-C215-C216-C217-C218",
                "C31-C32-C33-C34-C35-C36-C37-C38-C39=C310-C311-C312-C313-C314-C315-C316"
            ]
        "#,
        )?;

        let lip_sp = LipidSpecies::new("DOPE".to_owned(), descr, &src.select_all()?)?;
        println!("{lip_sp:?}");
        Ok(())
    }

    #[test]
    fn test_toml_descr() -> anyhow::Result<()> {
        let mut toml = String::new();
        std::fs::File::open("data/lipid_species.toml")?.read_to_string(&mut toml)?;
        let descr: PredefinedLipidSpecies = toml::from_str(&toml)?;
        println!("{descr:?}");
        Ok(())
    }

    #[test]
    fn test_toml2() -> anyhow::Result<()> {
        let toml: PredefinedLipidSpecies = toml::from_str(
            r#"
            [POPC]
            whole = "resname POPE"
            head = "name P N"
            mid = "name C21 C22"
            tails = [
                "C21-C22-C23-C24-C25-C26-C27-C28-C29=C210-C211-C212-C213-C214-C215-C216-C217-C218",
                "C31-C32-C33-C34-C35-C36-C37-C38-C39=C310-C311-C312-C313-C314-C315-C316"
            ]

            [POPE]
            whole = "resname POPE"
            head = "name P N"
            mid = "name C21 C22"
            tails = [
                "C21-C22-C23-C24-C25-C26-C27-C28-C29=C210-C211-C212-C213-C214-C215-C216-C217-C218",
                "C31-C32-C33-C34-C35-C36-C37-C38-C39=C310-C311-C312-C313-C314-C315-C316"
            ]
        "#,
        )?;
        println!("{toml:?}");
        Ok(())
    }

    #[test]
    fn test_whole() -> anyhow::Result<()> {
        let path = PathBuf::from("/home/semen/work/Projects/Misha/balanced/not_depleted");
        let src = Source::serial_from_file(path.join("inp_7.gro"))?;
        let z0 = src.select("not resname TIP3 POT CLA")?.center_of_mass()?.z;
        let mut toml = String::new();
        std::fs::File::open("data/lipid_species.toml")?.read_to_string(&mut toml)?;
        let mut memb = Membrane::new(src, &toml)?
            .with_global_normal(Vector3f::z())
            .with_order_type(OrderType::ScdCorr)
            .with_output_dir("../target/membr_test_results")?;

        let mut upper = vec![];
        let mut lower = vec![];
        for (id, lip) in memb.iter_lipids() {
            if lip.head_marker.z > z0 {
                upper.push(id);
            } else {
                lower.push(id);
            }
        }

        memb.add_group("upper", upper)?;
        memb.add_group("lower", lower)?;

        let traj = FileHandler::open(path.join("traj_comp.xtc"))?;
        for st in traj.into_iter().take(10) {
            //println!(">> {}",st.get_time());
            memb.set_state(st.into())?;
            memb.compute()?;
        }

        memb.finalize()?;

        Ok(())
    }
}
