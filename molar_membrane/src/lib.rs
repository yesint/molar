use anyhow::bail;
use molar::prelude::*;
use nalgebra::{DVector, OVector, VectorN};
use serde::Deserialize;
use sorted_vec::SortedSet;
use std::{
    collections::{hash_map, HashMap},
    sync::Arc,
};

#[cfg(not(test))]
use log::{info, warn}; // Use log crate when building application

#[cfg(test)]
use std::{println as info, println as warn};

pub struct Membrane {
    // source: Source<MutableSerial>,
    //species: Vec<(String, LipidSpecies)>,
    lipids: Vec<LipidMolecule>,
    groups: HashMap<String, LipidGroup>,
    global_normal: Option<Vector3f>,
    order_type: OrderType,
}

impl Membrane {
    pub fn new(source: Source<MutableSerial>, defstr: &str) -> anyhow::Result<Self> {
        // Load species descriptions
        let species_descr: HashMap<String, LipidSpeciesDescr> = toml::from_str(defstr)?;

        // Load lipids from provided source
        //let mut species: HashMap<String, Arc<LipidSpecies>> = Default::default();
        let mut lipids = vec![];
        for (name, descr) in species_descr.iter() {
            println!(">>> {}",descr.whole);
            if let Ok(lips) = source.select(&descr.whole) {
                println!(">>> {}",lips.len());
                let lips = lips.split_resindex::<Vec<_>>()?;
                // Use first lipid to create lipid species object
                info!("Creating {} '{}' lipids", lips.len(), name);
                let sp = Arc::new(LipidSpecies::new(name.clone(),descr.clone(), &lips[0])?);

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
                        stats: LipidProperties::new(&sp),
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
            self.groups.insert(
                gr_name.to_owned(),
                LipidGroup {
                    lipid_ids: ids,
                    stats: Default::default(),
                },
            );
        }
        Ok(())
    }

    pub fn compute(&mut self) -> anyhow::Result<()> {
        // Compute order for all lipids
        for lip in &mut self.lipids {
            lip.compute_order(self.order_type.clone(), &self.global_normal.unwrap());
        }
        // Add values to groups
        for gr in self.groups.values() {
            for i in &gr.lipid_ids {
                let sp_name = &self.lipids[*i].species.name;
                let order = &gr.stats.per_species[sp_name].order;
            }
        }
        Ok(())
    }
}

#[derive(Debug)]
pub struct LipidProperties {
    pub normal: Vector3f,
    pub area: f32,
    pub tilt: f32,
    pub order: Vec<Vec<f32>>,
    // mean_curv: f32,
    // gauss_curv: f32,
    // princ_curv: [f32; 2],
    // princ_curv_axes: [Vector3f; 2],
}

impl LipidProperties {
    pub fn new(species: &LipidSpecies) -> Self {
        // Allocate needed number of tails
        let mut order = Vec::with_capacity(species.tails.len());
        // Initialize each tail
        for t in &species.tails {
            order.push(vec![0.0; t.bond_orders.len()]);
        }
        Self {
            normal: Default::default(),
            area: Default::default(),
            tilt: Default::default(),
            order,
        }
    }
}

#[derive(Default, Debug)]
pub struct GroupProperties {
    pub per_species: HashMap<String, StatProperties>,
}

#[derive(Default, Debug)]
pub struct StatProperties {
    pub area: MeanStd,
    pub tilt: MeanStd,
    pub order: Vec<MeanStdVec>,
}

impl StatProperties {
    pub fn new(species: &LipidSpecies) -> Self {
        // Allocate needed number of tails
        let mut order = Vec::with_capacity(species.tails.len());
        // Initialize each tail
        for t in &species.tails {
            order.push(MeanStdVec::new(t.bond_orders.len()));
        }
        Self {
            area: Default::default(),
            tilt: Default::default(),
            order,
        }
    }

    pub fn add_lipid_stats(&mut self, lip_stats: LipidProperties) {

    }
}

#[derive(Default, Debug)]
pub struct MeanStd {
    x: f32,
    x2: f32,
    n: f32,
}

impl MeanStd {
    pub fn add(&mut self, val: f32) {
        self.x = self.x + val;
        self.x2 = self.x2 + val * val;
        self.n += 1.0;
    }

    pub fn compute(&self) -> anyhow::Result<(f32, f32)> {
        if self.n == 0.0 {
            bail!("no values accumulated in MeanStd");
        }
        let mean = self.x / self.n;
        let stddev = ((self.x2 / self.n) - (mean * mean)).sqrt();
        Ok((mean, stddev))
    }
}

#[derive(Debug)]
pub struct MeanStdVec {
    x: DVector<f32>,
    x2: DVector<f32>,
    n: f32,
}

impl MeanStdVec {
    pub fn new(sz: usize) -> Self {
        Self {
            x: DVector::zeros(sz),
            x2: DVector::zeros(sz),
            n: 0.0,
        }
    }

    pub fn add(&mut self, val: &DVector<f32>) -> anyhow::Result<()> {
        if val.len() != self.x.len() {
            bail!(
                "incomtible vector size {} in MeanStdVec::add, {} expected",
                val.len(),
                self.x.len()
            );
        }
        self.x = &self.x + val;
        self.x2 = &self.x2 + val.component_mul(val);
        self.n += 1.0;
        Ok(())
    }

    pub fn compute(&self) -> anyhow::Result<(DVector<f32>, DVector<f32>)> {
        if self.n == 0.0 {
            bail!("no values accumulated in MeanStd");
        }
        let mean = &self.x / self.n;
        let mut stddev = (&self.x2 / self.n) - mean.component_mul(&mean);
        stddev.apply(|v| *v = v.sqrt());
        Ok((mean, stddev))
    }
}

pub struct LipidGroup {
    lipid_ids: Vec<usize>,
    stats: GroupProperties,
}

#[derive(Clone, Debug, Deserialize)]
pub struct PredefinedLipidSpecies(HashMap<String, LipidSpeciesDescr>);

#[derive(Clone, Debug, Deserialize)]
pub struct LipidSpeciesDescr {
    pub whole: String,
    pub head: String,
    pub mid: String,
    pub tails: Vec<String>,
}

#[derive(Debug)]
pub struct LipidSpecies {
    name: String,
    descr: LipidSpeciesDescr,
    head_marker_offsets: SortedSet<usize>,
    mid_marker_offsets: SortedSet<usize>,
    tails: Vec<LipidTail>,
}

impl LipidSpecies {
    pub fn new<K: UserCreatableKind>(
        name: String,
        descr: LipidSpeciesDescr,
        lipid: &Sel<K>,
    ) -> anyhow::Result<Self> {
        let first_index = lipid.first_index();
        let mut tails = vec![];
        // Parse tails
        for t in &descr.tails {
            let mut names = vec![];
            let mut bond_orders = vec![];
            let mut cur = &t[..];

            while let Some(e) = cur.find(['-', '=']) {
                if cur[..e].is_empty() {
                    bail!("missing carbon atom name");
                }
                names.push(&cur[..e]);
                //println!("name: {} order: {}",&cur[..e], &cur[e..=e]);
                bond_orders.push(match &cur[e..=e] {
                    "-" => 1,
                    "=" => 2,
                    _ => unreachable!(),
                });
                cur = &cur[e + 1..];
            }

            if cur.is_empty() {
                bail!("missing last carbon atom name");
            }
            names.push(&cur);

            // Now find offsets for each atom name
            let mut offsets = vec![];
            for name in names {
                let atom = lipid.subsel(format!("name {name}"))?;
                if atom.len() > 1 {
                    bail!("more than one tail atom {name} in lipid");
                }
                offsets.push(atom.first_index() - first_index);
            }

            tails.push(LipidTail {
                descr: t.to_owned(),
                offsets,
                bond_orders,
            });
        }

        Ok(Self {
            name,
            descr: descr.clone(),
            head_marker_offsets: SortedSet::from_unsorted(
                lipid
                    .subsel(descr.head)?
                    .iter_index()
                    .map(|i| i - first_index)
                    .collect(),
            ),
            mid_marker_offsets: SortedSet::from_unsorted(
                lipid
                    .subsel(descr.mid)?
                    .iter_index()
                    .map(|i| i - first_index)
                    .collect(),
            ),
            tails,
        })
    }
}

#[derive(Debug)]
pub struct LipidTail {
    descr: String,
    offsets: Vec<usize>,
    bond_orders: Vec<u8>,
}

pub struct LipidMolecule {
    sel: Sel<MutableSerial>,
    species: Arc<LipidSpecies>,
    head_sel: Sel<MutableSerial>,
    mid_sel: Sel<MutableSerial>,
    tail_end_sel: Sel<MutableSerial>,
    tail_sels: Vec<Sel<MutableSerial>>,
    head_marker: Pos,
    mid_marker: Pos,
    tail_marker: Pos,
    stats: LipidProperties,
}

impl LipidMolecule {
    pub fn update_markers<K: UserCreatableKind>(&mut self) -> anyhow::Result<()> {
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
}

#[cfg(test)]
mod tests {
    use std::io::Read;

    use molar::prelude::*;

    use crate::{LipidSpecies, LipidSpeciesDescr, Membrane, PredefinedLipidSpecies};

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
        let lip_sp = LipidSpecies::new("DOPE".to_owned(),descr, &pope)?;
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

        let lip_sp = LipidSpecies::new("DOPE".to_owned(),descr, &src.select_all()?)?;
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
        let src = Source::serial_from_file("/home/semen/work/Projects/Misha/PG_flipping/last.gro")?;
        let z0 = src
            .select("not resname TIP3 POT CLA")?
            .center_of_mass()?
            .z;
        let mut toml = String::new();
        std::fs::File::open("data/lipid_species.toml")?.read_to_string(&mut toml)?;
        let mut memb = Membrane::new(src, &toml)?
            .with_global_normal(Vector3f::z())
            .with_order_type(OrderType::ScdCorr);

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

        Ok(())
    }
}
