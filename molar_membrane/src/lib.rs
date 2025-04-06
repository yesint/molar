use anyhow::{bail, Context};
use molar::prelude::*;
use nalgebra::DVector;
use std::{
    collections::HashMap,
    f32::consts::FRAC_PI_2,
    path::{Path, PathBuf},
    sync::Arc,
};

#[cfg(not(test))]
use log::{info, warn}; // Use log crate when building application

#[cfg(test)]
use std::{println as info, println as warn};

mod stats;
use stats::{GroupProperties, StatProperties};

mod lipid_molecule;
use lipid_molecule::LipidMolecule;

mod lipid_species;
use lipid_species::{LipidSpecies, LipidSpeciesDescr};

mod surface;
use surface::*;

pub struct Membrane {
    source: Source<MutableSerial>,
    lipids: Vec<LipidMolecule>,
    surface: Surface,
    groups: HashMap<String, LipidGroup>,
    // Options
    global_normal: Option<Vector3f>,
    order_type: OrderType,
    output_dir: PathBuf,
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

                    let mut order = Vec::with_capacity(sp.tails.len());
                    for t in &sp.tails {
                        order.push(DVector::from_element(t.bond_orders.len() - 1, 0.0));
                    }

                    let tail_head_vec = tail_marker - head_marker;

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

        Ok(Self {
            source,
            lipids,
            groups: Default::default(),
            global_normal: None,
            order_type: OrderType::ScdCorr,
            output_dir: PathBuf::from("membrane_results"),
            surface,
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
        // Compute patches
        self.compute_patches(2.5)?;

        // Get initial normals
        self.compute_initial_normals();

        // Initialize markers and set all points as valid
        for i in 0..self.lipids.len() {
            self.surface.nodes[i].marker = self.lipids[i].head_marker;
            self.surface.nodes[i].valid = true;
        }

        // Do smoothing
        let mut iter = 0;
        //const TOL: f32 = 1e-3;
        const MAX_ITER: u8 = 5;

        loop {
            self.surface.smooth();

            self.write_vmd_visualization(format!(
                "/home/semen/work/Projects/Misha/PG_flipping/iter_{iter}.tcl"
            ))
            .unwrap();

            iter += 1;

            if iter >= MAX_ITER {
                break;
            }
        }

        // Compute order for all lipids
        for i in 0..self.lipids.len() {
            let norm = if let Some(n) = &self.global_normal {
                n
            } else {
                &self.surface.nodes[i].normal
            };
            self.lipids[i].compute_order(self.order_type.clone(), norm);
        }
        // Add stats to groups
        for gr in self.groups.values_mut() {
            for lip_id in &gr.lipid_ids {
                gr.stats.add_single_lipid_stats(&self.lipids[*lip_id], &self.surface.nodes[*lip_id])?;
            }
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

    pub fn set_state(&mut self, st: State) -> anyhow::Result<()> {
        self.source.set_state(st)?;
        // for lip in &mut self.lipids {
        //     lip.set_state(st.clone())?;
        //     lip.sel.unwrap_simple()?;
        //     lip.update_markers()?;
        // }
        Ok(())
    }

    fn compute_patches(&mut self, d: f32) -> anyhow::Result<()> {
        let ind: Vec<(usize, usize)> = distance_search_single_pbc(
            d,
            self.lipids.iter().map(|l| &l.head_marker),
            0..self.lipids.len(),
            self.lipids[0].sel.require_box()?,
            PBC_FULL,
        );

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

//===========================================
struct VmdVisual {
    buf: String,
}

impl VmdVisual {
    fn new() -> Self {
        Self { buf: String::new() }
    }

    fn save_to_file(&self, fname: impl AsRef<Path>) -> anyhow::Result<()> {
        use std::io::Write;
        let mut f = std::fs::File::create(fname)?;
        writeln!(f, "{}", self.buf)?;
        Ok(())
    }

    fn sphere(&mut self, point: &Pos, radius: f32, color: &str) {
        use std::fmt::Write;
        writeln!(self.buf, "draw color {color}").unwrap();
        writeln!(
            self.buf,
            "draw sphere \"{} {} {}\" radius {radius} resolution 12",
            point.x * 10.0,
            point.y * 10.0,
            point.z * 10.0
        )
        .unwrap();
    }

    fn arrow(&mut self, point: &Pos, dir: &Vector3f, color: &str) {
        use std::fmt::Write;

        const LENGTH: f32 = 5.0;

        let p1 = point * 10.0;
        let p2 = p1 + dir * 0.5 * LENGTH;
        let p3 = p1 + dir * 0.7 * LENGTH;

        writeln!(self.buf, "draw color {color}").unwrap();

        writeln!(
            self.buf,
            "draw cylinder \"{} {} {}\" \"{} {} {}\" radius 0.2 resolution 12",
            p1.x, p1.y, p1.z, p2.x, p2.y, p2.z
        )
        .unwrap();

        writeln!(
            self.buf,
            "draw cone \"{} {} {}\" \"{} {} {}\" radius 0.4 resolution 12\n",
            p2.x, p2.y, p2.z, p3.x, p3.y, p3.z
        )
        .unwrap();
    }

    fn cylinder(&mut self, point1: &Pos, point2: &Pos, color: &str) {
        use std::fmt::Write;

        let p1 = point1 * 10.0;
        let p2 = point2 * 10.0;

        writeln!(self.buf, "draw color {color}").unwrap();

        writeln!(
            self.buf,
            "draw cylinder \"{} {} {}\" \"{} {} {}\" radius 0.2 resolution 12",
            p1.x, p1.y, p1.z, p2.x, p2.y, p2.z
        )
        .unwrap();
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
                stats.per_species.get_mut(&memb.lipids[*id].species.name).unwrap().num_lip += 1;
            }
        }

        for (sp, stats) in &stats.per_species {
            info!("\t{}: {}", sp, stats.num_lip);
        }
        Self { lipid_ids, stats }
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
        // let path = PathBuf::from("tests");
        // let src = Source::serial_from_file(path.join("membr.gro"))?;

        let path = PathBuf::from("/home/semen/work/Projects/Misha/PG_flipping");
        let src = Source::serial_from_file(path.join("inp_7.gro"))?;
        
        let z0 = 5.6;//src.select("not resname TIP3 POT CLA /A+/ /B+/")?.center_of_mass()?.z;
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
            println!(">> {}", st.get_time());
            memb.set_state(st.into())?;
            memb.compute()?;
        }

        memb.finalize()?;

        Ok(())
    }

    #[test]
    fn test_vmd_vis() -> anyhow::Result<()> {
        let path = PathBuf::from("tests");
        let src = Source::serial_from_file(path.join("membr.gro"))?;
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

        memb.compute()?;
        memb.finalize()?;
        memb.write_vmd_visualization("../target/vis.tcl")?;

        Ok(())
    }

    #[test]
    fn periodic_sel_test() {
        let path = PathBuf::from("tests");
        let src = Source::serial_from_file(path.join("membr.gro")).unwrap();
        let pbox = src.get_box().unwrap();

        let p = src.select("name P").unwrap();
        let pairs: Vec<(usize, usize)> =
            distance_search_single_pbc(2.5, p.iter_pos(), p.iter_index(), pbox, PBC_FULL);

        let r162 = src.select("resid 162 and name P").unwrap();
        let p1 = r162.first_particle().pos;
        let r123 = src.select("resid 123 and name P").unwrap();
        let p2 = r123.first_particle().pos;

        let d = pbox.distance(p1, p2, PBC_FULL);
        println!("d={d}");
    }
}
