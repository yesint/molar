use anyhow::{bail, Context};
use molar::prelude::*;
use nalgebra::DVector;
use serde::Deserialize;
use std::{
    collections::HashMap,
    f32::consts::FRAC_PI_2,
    path::{Path, PathBuf},
    rc::Rc,
};

#[cfg(not(test))]
use log::{info, warn}; // Use log crate when building application

#[cfg(test)]
use std::{println as info, println as warn};

mod stats;
use stats::SpeciesStats;

mod lipid_molecule;
pub use lipid_molecule::LipidMolecule;

mod lipid_species;
use lipid_species::{LipidSpecies, LipidSpeciesDescr};

mod surface;
use surface::*;

mod vmd_visual;
use vmd_visual::VmdVisual;

mod lipid_group;
use lipid_group::LipidGroup;

pub struct Membrane<K> {
    lipids: Vec<LipidMolecule<K>>,
    surface: Surface,
    groups: HashMap<String, LipidGroup>,
    species: Vec<Rc<LipidSpecies>>,
    // Options
    options: MembraneOptions,
}

#[derive(Deserialize)]
#[serde(default)]
struct MembraneOptions {
    global_normal: Option<Vector3f>,
    order_type: OrderType,
    output_dir: PathBuf,
    max_smooth_iter: usize,
    cutoff: f32,
    lipids: HashMap<String, LipidSpeciesDescr>,
    groups: Vec<String>,
    sel: String,
}

impl Default for MembraneOptions {
    fn default() -> Self {
        MembraneOptions {
            global_normal: None,
            cutoff: 2.5,
            max_smooth_iter: 1,
            order_type: OrderType::ScdCorr,
            output_dir: ".".into(),
            lipids: Default::default(),
            groups: vec![],
            sel: "all".into(),
        }
    }
}

impl<K: MutableKind+UserCreatableKind> Membrane<K> {
    pub fn new(source: &Source<K>, optstr: &str) -> anyhow::Result<Self> {
        // Load options
        info!("Processing membrane options...");
        let options: MembraneOptions = toml::from_str(optstr)?;

        // Create working selection
        let source_sel = source.select(&options.sel)?;

        // Load lipids from provided source
        let mut lipids = vec![];
        let mut species = vec![];
        for (name, descr) in options.lipids.iter() {
            if let Ok(lips) = source_sel.subsel(&descr.whole) {
                let lips = lips.split_resindex::<Vec<_>>()?;
                // Use first lipid to create lipid species object
                info!("Creating {} '{}' lipids", lips.len(), name);
                let sp = Rc::new(LipidSpecies::new(name.clone(), descr.clone(), &lips[0])?);
                species.push(Rc::clone(&sp));

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

                    let mut order = Vec::with_capacity(sp.tails.len());
                    for t in &sp.tails {
                        order.push(DVector::from_element(t.bond_orders.len() - 1, 0.0));
                    }

                    let tail_head_vec = tail_marker - head_marker;

                    lipids.push(LipidMolecule {
                        sel: lip, // Selection is moved to the lipid
                        species: Rc::clone(&sp),
                        head_sel,
                        mid_sel,
                        tail_end_sel,
                        tail_sels,
                        head_marker,
                        mid_marker,
                        tail_marker,
                        order,
                        tail_head_vec,
                    });
                }
            } else {
                // No such lipid species found
                warn!("No '{name}' lipids found");
            }
        }

        // We need parallel-friendly data structure for surface computations
        let surface = Surface {
            pbox: lipids[0].sel.require_box()?.clone(),
            nodes: lipids
                .iter()
                .map(|l| SurfNode {
                    marker: l.head_marker,
                    ..Default::default()
                })
                .collect(),
        };

        // If there are any groups in the input toml, create them
        let mut groups: HashMap<String, LipidGroup> = Default::default();
        for gr in &options.groups {
            groups.insert(gr.to_string(), Default::default());
        }

        Ok(Self {
            lipids,
            groups,
            surface,
            species,
            options,
        })
    }

    pub fn with_groups(mut self, group_names: &[impl AsRef<str>]) -> Self {
        for gr in group_names {
            self.groups
                .insert(gr.as_ref().to_string(), Default::default());
        }
        self
    }

    pub fn with_global_normal(mut self, global_normal: impl Into<Option<Vector3f>>) -> Self {
        self.options.global_normal = global_normal.into();
        self
    }

    pub fn with_order_type(mut self, order_type: OrderType) -> Self {
        self.options.order_type = order_type;
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

        self.options.output_dir = path.to_path_buf();
        info!("will use output directory '{}'", path.display());
        Ok(self)
    }

    pub fn with_cutoff(mut self, cutoff: f32) -> Self {
        self.options.cutoff = cutoff;
        self
    }

    pub fn with_max_iter(mut self, max_iter: usize) -> Self {
        self.options.max_smooth_iter = max_iter;
        self
    }

    pub fn iter_lipids(&self) -> impl Iterator<Item = (usize, &LipidMolecule<K>)> {
        self.lipids.iter().enumerate()
    }

    pub fn reset_groups(&mut self) {
        for gr in self.groups.values_mut() {
            gr.lipid_ids.clear();
            gr.stats.per_species.clear();
        }
    }

    pub fn add_lipids_to_group(
        &mut self,
        gr_name: impl AsRef<str>,
        ids: &Vec<usize>,
    ) -> anyhow::Result<()> {
        let gr_name = gr_name.as_ref();
        if let Some(gr) = self.groups.get_mut(gr_name) {
            for id in ids {
                // Check if lipid id is in range
                if *id >= self.lipids.len() {
                    bail!(
                        "lipid id {} is out of bounds 0:{}",
                        id,
                        self.lipids.len() - 1
                    );
                }
                // Add this lipid id
                gr.lipid_ids.push(*id);
                // Initialize statistics for this lipid species if not yet done
                let sp = &self.lipids[*id].species.name;
                if !gr.stats.per_species.contains_key(sp) {
                    gr.stats.per_species.insert(
                        sp.to_string(),
                        SpeciesStats::new(
                            &self.lipids[*id].species,
                            self.species.iter().map(|sp| sp.name.to_owned()),
                        ),
                    );
                }
            }
        } else {
            bail!("No such group: {gr_name}");
        }

        Ok(())
    }

    pub fn compute(&mut self) -> anyhow::Result<()> {
        // Compute patches
        self.compute_patches(self.options.cutoff)?;

        // Get initial normals
        self.compute_initial_normals();

        // Initialize markers and set all points as valid
        for i in 0..self.lipids.len() {
            self.surface.nodes[i].marker = self.lipids[i].head_marker;
            self.surface.nodes[i].valid = true;
        }

        // Do smoothing
        let mut iter = 0;
        loop {
            self.surface.smooth();

            // self.write_vmd_visualization(format!(
            //     "/home/semen/work/Projects/Misha/PG_flipping/iter_{iter}.tcl"
            // ))
            // .unwrap();

            iter += 1;

            if iter >= self.options.max_smooth_iter {
                break;
            }
        }

        // Compute order for all lipids
        for i in 0..self.lipids.len() {
            let norm = if let Some(n) = &self.options.global_normal {
                n
            } else {
                &self.surface.nodes[i].normal
            };
            self.lipids[i].compute_order(self.options.order_type.clone(), norm);
        }

        // Update group stats
        for gr in self.groups.values_mut() {
            gr.frame_update(&self.lipids, &self.surface)?;
        }

        Ok(())
    }

    fn compute_initial_normals(&mut self) {
        // Compute tail-to-head vectors for all lipids
        // This is unwrapped with the same lipid
        for lip in &mut self.lipids {
            lip.tail_head_vec = (lip.head_marker - lip.tail_marker).normalize();
        }

        // First pass - average of tail-to-head distances in each patch
        for i in 0..self.lipids.len() {
            self.surface.nodes[i].normal = self.surface.nodes[i]
                .patch_ids
                .iter()
                .map(|l| &self.lipids[*l].tail_head_vec)
                .sum::<Vector3f>()
                .normalize();
        }

        // Second pass - average of vectors from the first pass in eac patch
        for i in 0..self.lipids.len() {
            self.surface.nodes[i].normal = self.surface.nodes[i]
                .patch_ids
                .iter()
                .map(|l| &self.surface.nodes[*l].normal)
                .sum::<Vector3f>()
                .normalize();
        }

        // Correct normal orientations if needed
        for i in 0..self.lipids.len() {
            if self.surface.nodes[i]
                .normal
                .angle(&self.lipids[i].tail_head_vec)
                > FRAC_PI_2
            {
                self.surface.nodes[i].normal *= -1.0;
            }
        }
    }

    pub fn finalize(&self) -> anyhow::Result<()> {
        info!(
            "Writing results to directory '{}'",
            self.options.output_dir.display()
        );
        // Write results for groups
        for (gr_name, gr) in &self.groups {
            info!("\tGroup '{gr_name}'...");
            gr.stats
                .save_group_stats(self.options.output_dir.as_path(), gr_name)?;
            gr.stats
                .save_order_to_file(self.options.output_dir.as_path(), gr_name)?;
        }
        Ok(())
    }

    pub fn set_state(&mut self, st: impl Into<Holder<State,K>>) -> anyhow::Result<()> {
        let st: Holder<State,K> = st.into();
        // Go over all lipids and set their "shallow" states
        for lip in &mut self.lipids {
            lip.sel.set_state(st.new_ref())?;
            lip.head_sel.set_state(st.new_ref())?;
            lip.mid_sel.set_state(st.new_ref())?;
            lip.tail_end_sel.set_state(st.new_ref())?;
            for t in &mut lip.tail_sels {
                t.set_state(st.new_ref())?;
            }
        }
        Ok(())
    }

    pub fn get_state(&self) -> Holder<State, K> {
        self.lipids.first().unwrap().sel.get_state()
    }

    fn compute_patches(&mut self, d: f32) -> anyhow::Result<()> {
        let ind: Vec<(usize, usize)> = distance_search_single_pbc(
            d,
            self.lipids.iter().map(|l| &l.head_marker),
            0..self.lipids.len(),
            self.lipids[0].sel.require_box()?,
            PBC_FULL,
        );

        // Clear all patches first!
        for node in self.surface.nodes.iter_mut() {
            node.patch_ids.clear();
        }
        // Add ids to patches
        for (i, j) in ind {
            self.surface.nodes[i].patch_ids.push(j);
            self.surface.nodes[j].patch_ids.push(i);
        }
        Ok(())
    }

    pub fn write_vmd_visualization(&self, fname: impl AsRef<Path>) -> anyhow::Result<()> {
        let mut vis = VmdVisual::new();

        for i in 0..self.lipids.len() {
            let lip = &self.surface.nodes[i];
            // Initial lipid marker
            vis.sphere(&self.lipids[i].head_marker, 0.8, "white");
            // tail-head vector
            vis.arrow(&lip.marker, &self.lipids[i].tail_head_vec, "yellow");

            // Fitted lipid marker
            vis.sphere(&lip.marker, 0.8, "red");
            // fitted normal
            vis.arrow(&lip.marker, &lip.normal, "orange");

            // Voronoi cell
            let n = lip.voro_vertexes.len();
            for i in 0..n {
                let p1 = lip.voro_vertexes[i];
                let p2 = lip.voro_vertexes[(i + 1) % n];
                vis.cylinder(&p1, &p2, "green");
            }
            //vis.cylinder(&lip.voro_vertexes[n - 1], &lip.voro_vertexes[0], "green");

            // Fitted patch
            for p in &lip.fitted_patch_points {
                vis.sphere(p, 0.3, "green");
            }
        }

        vis.save_to_file(fname)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::{collections::HashMap, io::Read, path::PathBuf};

    use molar::prelude::*;

    use crate::{LipidSpecies, LipidSpeciesDescr, Membrane};

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

    #[allow(dead_code)]
    #[derive(Clone, Debug, serde::Deserialize)]
    struct TestInp {
        cutoff: f32,
        lipids: HashMap<String, LipidSpeciesDescr>,
    }

    #[test]
    fn test_toml2() -> anyhow::Result<()> {
        let toml: TestInp = toml::from_str(
            r#"
            cutoff = 2.5
            [lipids.POPC]
            whole = "resname POPE"
            head = "name P N"
            mid = "name C21 C22"
            tails = [
                "C21-C22-C23-C24-C25-C26-C27-C28-C29=C210-C211-C212-C213-C214-C215-C216-C217-C218",
                "C31-C32-C33-C34-C35-C36-C37-C38-C39=C310-C311-C312-C313-C314-C315-C316"
            ]

            [lipids.POPE]
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
        let path = PathBuf::from("tests");
        let src = Source::serial_from_file(path.join("membr.gro"))?;

        // let path = PathBuf::from("/home/semen/work/Projects/Misha/PG_flipping");
        // let src = Source::serial_from_file(path.join("inp_7.gro"))?;

        let z0 = 5.6; //src.select("not resname TIP3 POT CLA /A+/ /B+/")?.center_of_mass()?.z;
        let mut toml = String::new();
        std::fs::File::open("data/inp.toml")?.read_to_string(&mut toml)?;
        let mut memb = Membrane::new(&src, &toml)?;
        //.with_groups(&["upper", "lower"]);
        // .with_global_normal(Vector3f::z())
        // .with_order_type(OrderType::ScdCorr)
        // .with_output_dir("../target/membr_test_results")?
        // .with_cutoff(2.5)
        // .with_max_iter(1);

        let mut upper = vec![];
        let mut lower = vec![];
        for (id, lip) in memb.iter_lipids() {
            if lip.head_marker.z > z0 + 1.0 {
                upper.push(id);
            } else if lip.head_marker.z < z0 - 1.0 {
                lower.push(id);
            }
        }

        memb.add_lipids_to_group("upper", &upper)?;
        memb.add_lipids_to_group("lower", &lower)?;

        let traj = FileHandler::open(path.join("traj_comp.xtc"))?;
        for st in traj.into_iter().take(10) {
            println!(">> {}", st.get_time());
            memb.set_state(st)?;
            memb.compute()?;
        }

        memb.finalize()?;

        Ok(())
    }

    #[test]
    fn test_vmd_vis() -> anyhow::Result<()> {
        let path = PathBuf::from("tests");
        let src = Source::serial_from_file(path.join("membr.gro"))?;
        let z0 = 5.6; //src.select("not resname TIP3 POT CLA")?.center_of_mass()?.z;
        let mut toml = String::new();
        std::fs::File::open("data/inp.toml")?.read_to_string(&mut toml)?;

        let mut memb = Membrane::new(&src, &toml)?.with_max_iter(1);
        // .with_global_normal(Vector3f::z())
        // .with_order_type(OrderType::ScdCorr)
        // .with_output_dir("../target/membr_test_results")?;

        let mut upper = vec![];
        let mut lower = vec![];
        for (id, lip) in memb.iter_lipids() {
            if lip.head_marker.z > z0 + 1.0 {
                upper.push(id);
            } else if lip.head_marker.z < z0 - 1.0 {
                lower.push(id);
            }
            if lip.sel.first_atom().resname == "POGL" {
                println!("{}", lip.head_marker);
            }
        }

        memb.add_lipids_to_group("upper", &upper)?;
        memb.add_lipids_to_group("lower", &lower)?;

        memb.compute()?;
        memb.finalize()?;
        memb.write_vmd_visualization("../target/vis.tcl")?;

        Ok(())
    }
}
