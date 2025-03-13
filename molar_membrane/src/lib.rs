use anyhow::bail;
use log::warn;
use molar::prelude::*;
use serde::Deserialize;
use sorted_vec::SortedSet;
use std::{collections::HashMap, sync::Arc};

pub struct Membrane {
    source: Source<MutableSerial>,
    species: HashMap<String, Arc<LipidSpecies>>,
    lipids: Vec<LipidMolecule>,
    groups: Vec<LipidGroup>,
}

impl Membrane {
    pub fn new(source: Source<MutableSerial>) -> Self {
        Self {
            source,
            species: Default::default(),
            lipids: vec![],
            groups: vec![],
        }
    }

    pub fn species_from_str(source: Source<MutableSerial>, defstr: &str) -> anyhow::Result<Self> {
        // Load species descriptions
        let species_descr: HashMap<String, LipidSpeciesDescr> = toml::from_str(defstr)?;
        // Load lipids from provided source
        let mut species: HashMap<String, Arc<LipidSpecies>> = Default::default();
        let mut lipids = vec![];
        for (name, descr) in species_descr.iter() {
            let lips = source
                .select_str(&descr.whole)?
                .split_resindex::<Vec<_>>()?;
            if !lips.is_empty() {
                // Use first lipid to create lipid species object
                log::info!("Creating {} '{}' lipids", lips.len(), name);
                species.insert(
                    name.to_owned(),
                    Arc::new(LipidSpecies::new(descr.clone(), &lips[0])?),
                );
                // Now create individual lipids
                
                for lip in lips {
                    let sp = &species[name];
                    let head_sel = lip.subsel_iter(sp.head_marker_offsets.iter().cloned())?;
                    let mid_sel = lip.subsel_iter(sp.mid_marker_offsets.iter().cloned())?;
                    
                    let mut tail_sels = vec![];
                    let mut tail_ends = vec![];
                    for t in sp.tails.iter() {
                        tail_sels.push(lip.subsel_iter(t.offsets.iter().cloned())?);
                        tail_ends.push(*t.offsets.last().unwrap());
                    }
                    let tail_end_sel = lip.subsel_iter(tail_ends.iter().cloned())?;
                    
                    lipids.push(LipidMolecule {
                        sel: lip, // Selection is moved to the lipid
                        species: Arc::clone(sp),
                        head_sel,
                        mid_sel,
                        tail_end_sel,
                        tail_sels,
                        head_marker: Default::default(),
                        mid_marker: Default::default(),
                        tail_marker: Default::default(),
                    });
                }
            }
        }
        Ok(Self {
            source,
            species,
            lipids,
            groups: Default::default(),
        })
    }
}

pub struct LipidGroup {
    name: String,
    lipids: Vec<usize>,
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
    descr: LipidSpeciesDescr,
    head_marker_offsets: SortedSet<usize>,
    mid_marker_offsets: SortedSet<usize>,
    tails: Vec<LipidTail>,
}

impl LipidSpecies {
    pub fn new<K: UserCreatableKind>(
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
                let atom = lipid.subsel_str(format!("name {name}"))?;
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
            descr: descr.clone(),
            head_marker_offsets: SortedSet::from_unsorted(
                lipid
                    .subsel_str(descr.head)?
                    .iter_index()
                    .map(|i| i - first_index)
                    .collect(),
            ),
            mid_marker_offsets: SortedSet::from_unsorted(
                lipid
                    .subsel_str(descr.mid)?
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

#[derive(Debug)]
pub struct LipidMolecule<K: UserCreatableKind> {
    sel: Sel<K>,
    species: Arc<LipidSpecies>,
    markers: LipidMarkers,
    cutoff_neib: Vec<usize>,
    voronoi_neib: Vec<usize>,
    properties: LipidProperties,
}

#[derive(Debug,Default)]
struct LipidMarkers {
    // Markers
    head_marker: Pos,
    mid_marker: Pos,
    tail_marker: Pos,
}

/// Instanteneous lipid properties
#[derive(Default,Debug)]
pub struct LipidProperties {
    // Normals and orientations
    pub normal: Vector3f,
    pub head_tail_vec: Vector3f,
    pub tilt_angle: f32,
    // Curvatures
    pub mean_curvature: f32,
    pub gaussian_curvature: f32,
    pub prinical_curvature_vecs: [Vector3f;2],
    pub prinical_curvatures: [f32;2],
    // Tails properties
    pub tail_orders: Vec<Vec<f32>>,
    pub tail_dihedrals: Vec<Vec<f32>>,
}

impl<K: UserCreatableKind> LipidMolecule<K> {
    pub fn new(species: Arc<LipidSpecies>, sel: Sel<K>) -> anyhow::Result<Self> {
        let head_marker = sel
            .subsel_iter(species.head_marker_offsets.iter().cloned())?
            .center_of_mass_pbc()?;
        let mid_marker = sel
            .subsel_iter(species.mid_marker_offsets.iter().cloned())?
            .center_of_mass_pbc()?;

        let tail_marker;
        if !species.tails.is_empty() {
            // Collect ends of all tails
            let mut end_ind = Vec::with_capacity(species.tails.len());
            for tail in species.tails.iter() {
                end_ind.push(tail.offsets.last().unwrap());
            }
            tail_marker = sel
                .subsel_iter(end_ind.into_iter().cloned())?
                .center_of_mass_pbc()?;
        } else {
            tail_marker = mid_marker;
        }
        Ok(Self {
            sel,
            species,
            markers: LipidMarkers{
                head_marker,
                mid_marker,
                tail_marker,
            },
            cutoff_neib: vec![],
            voronoi_neib: vec![],
            properties: Default::default(),
        })
    }

    pub fn update(&mut self, sel: &Sel<K>) -> anyhow::Result<()> {
        self.markers.head_marker = sel
            .subsel_iter(self.species.head_marker_offsets.iter().cloned())?
            .center_of_mass_pbc()?;
        self.markers.mid_marker = sel
            .subsel_iter(self.species.mid_marker_offsets.iter().cloned())?
            .center_of_mass_pbc()?;

        if !self.species.tails.is_empty() {
            // Collect ends of all tails
            let mut end_ind = Vec::with_capacity(self.species.tails.len());
            for tail in self.species.tails.iter() {
                end_ind.push(tail.offsets.last().unwrap());
            }
            self.markers.tail_marker = sel
                .subsel_iter(end_ind.into_iter().cloned())?
                .center_of_mass_pbc()?;
        } else {
            self.markers.tail_marker = self.markers.mid_marker;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::io::Read;

    use molar::prelude::*;

    use crate::{LipidSpecies, LipidSpeciesDescr, PredefinedLipidSpecies};

    #[test]
    fn test_descr() -> anyhow::Result<()> {
        let src = Source::serial_from_file("../../tests/membr.gro")?;
        let pope = src.select_str("resid 60")?;
        let descr = LipidSpeciesDescr {
            whole: "resname POPE".into(),
            head: "name P N".into(),
            mid: "name C21 C22".into(),
            tails: vec![
                "C21-C22-C23-C24-C25-C26-C27-C28-C29=C210-C211-C212-C213-C214-C215-C216-C217-C218"
                    .into(),
                "C21-C22-C23-C24-C25-C26-C27-C28-C29=C210-C211-C212-C213-C214-C215-C216-C217-C218"
                    .into(),
                "C31-C32-C33-C34-C35-C36-C37-C38-C39=C310-C311-C312-C313-C314-C315-C316".into(),
            ],
        };
        let lip_sp = LipidSpecies::new(descr, &pope)?;
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
        let top = FileHandler::open("../../tests/POPE.itp")?.read_topology()?;
        let n = top.num_atoms();
        let src = Source::new_serial(top.into(), State::new_fake(n).into())?;

        let descr: LipidSpeciesDescr = toml::from_str(
            r#"
            name = "POPE"
            whole_sel_str = "resname POPE"
            head_marker_subsel_str = "name P N"
            mid_marker_subsel_str = "name C21 C22"
            tails_descr = [
                "C21-C22-C23-C24-C25-C26-C27-C28-C29=C210-C211-C212-C213-C214-C215-C216-C217-C218",
                "C31-C32-C33-C34-C35-C36-C37-C38-C39=C310-C311-C312-C313-C314-C315-C316"
            ]
        "#,
        )?;
        "#,
        )?;

        let lip_sp = LipidSpecies::new(descr, &src.select_all()?)?;
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
        "#,
        )?;
        println!("{toml:?}");
        Ok(())
    }
}
