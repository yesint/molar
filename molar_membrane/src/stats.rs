use anyhow::bail;
use log::info;
use molar::core::{MutableKind, UserCreatableKind};
use nalgebra::DVector;
use serde::{ser::SerializeStruct, Serialize, Serializer};
use std::collections::HashMap;
use std::path::Path;

use crate::lipid_species::LipidSpecies;
use crate::{LipidMolecule, Surface};

#[derive(Default, Debug)]
pub struct GroupProperties {
    pub per_species: HashMap<String, SpeciesStats>,
}

impl GroupProperties {
    pub(crate) fn save_order_to_file(
        &self,
        dir: &Path,
        gr_name: impl AsRef<str>,
    ) -> anyhow::Result<()> {
        for (sp, stat) in &self.per_species {
            info!("\tWriting order for '{sp}'");
            stat.save_order_to_file(dir, format!("gr_{}_order_{}.dat", gr_name.as_ref(), sp))?;
        }
        Ok(())
    }

    pub fn save_group_stats(&self, dir: &Path, gr_name: impl AsRef<str>) -> anyhow::Result<()> {
        use std::fmt::Write as FmtWrite;
        use std::io::Write as StdWrite;

        // Write basic statistsics
        info!("\tWriting basic statistics...");
        let mut s = "#species\tnum\tnum_std\tarea\tarea_std\ttilt\ttilt_std\n".to_string();
        for (sp, stat) in &self.per_species {
            let area = stat.area.compute()?;
            let tilt = stat.tilt.compute()?;
            let num = stat.num_lip.compute()?;
            writeln!(
                s,
                "{sp}\t{:>8.3}\t{:>8.3}\t{:>8.3}\t{:>8.3}\t{:>8.3}\t{:>8.3}",
                num.mean, num.stddev, area.mean, area.stddev, tilt.mean, tilt.stddev
            )?;
        }

        let mut f = std::fs::File::create(dir.join(format!("gr_{}_stats.dat", gr_name.as_ref(),)))?;
        write!(f, "{}", s)?;

        // Write neighbours frequencies
        s = "".to_string();
        for (sp, stat) in &self.per_species {
            let nneib = stat.num_neib.compute()?;
            writeln!(s, "{sp}:\t\t{:>8.3}\t{:>8.3}", nneib.mean, nneib.stddev)?;
            for (nsp, counter) in &stat.neib_species {
                let res = counter.compute()?;
                writeln!(s, "\t{nsp}\t{:>8.3}\t{:>8.3}", res.mean, res.stddev)?;
            }
            s.push('\n');
        }

        let mut f =
            std::fs::File::create(dir.join(format!("gr_{}_neib_stats.dat", gr_name.as_ref(),)))?;
        write!(f, "{}", s)?;

        Ok(())
    }
}

#[derive(Debug)]
pub(crate) struct SpeciesStats {
    pub num_lip: MeanStd,
    pub area: MeanStd,
    pub tilt: MeanStd,
    pub order: Vec<MeanStdVec>,
    pub neib_species: HashMap<String, MeanStd>,
    pub num_neib: MeanStd,

    // Accumulated values for a single frame update
    num_lip_cur: usize,
    neib_species_cur: HashMap<String, usize>,
}

impl SpeciesStats {
    pub fn new(
        cur_species: &LipidSpecies,
        all_species_names: impl Iterator<Item = String>,
    ) -> Self {
        let mut order = Vec::with_capacity(cur_species.tails.len());
        for t in &cur_species.tails {
            order.push(MeanStdVec::new(t.bond_orders.len() - 1));
        }

        let neib_species = all_species_names
            .map(|sp| (sp.to_owned(), Default::default()))
            .collect::<HashMap<_, _>>();

        let neib_species_cur = neib_species.keys()
            .map(|sp| (sp.to_owned(), 0))
            .collect::<HashMap<_, _>>();

        Self {
            area: MeanStd::default(),
            tilt: MeanStd::default(),
            order,
            num_lip: MeanStd::default(),
            num_neib: MeanStd::default(),
            neib_species,
            num_lip_cur: 0,
            neib_species_cur,
        }
    }

    pub(crate) fn init_frame_update(&mut self) {
        self.num_lip_cur = 0;
        for v in self.neib_species_cur.values_mut() {
            *v = 0;
        }
    }

    pub fn add_single_lipid_stats<K: UserCreatableKind+MutableKind>(
        &mut self,
        id: usize,
        lipids: &Vec<LipidMolecule<K>>,
        surf: &Surface,
    ) -> anyhow::Result<()> {
        if surf.nodes[id].valid {
            self.area.add(surf.nodes[id].area);

            self.tilt.add(
                surf.nodes[id]
                    .normal
                    .angle(&lipids[id].tail_head_vec)
                    .to_degrees(),
            );

            for tail in 0..lipids[id].order.len() {
                self.order[tail].add(&lipids[id].order[tail])?;
            }

            self.num_neib.add(surf.nodes[id].neib_ids.len() as f32);

            // Update lipid counter
            self.num_lip_cur += 1;

            // Update neighbours
            //let cur_sp_name = &lipids[id].species.name;
            for nid in &surf.nodes[id].neib_ids {
                let neib_sp_name = &lipids[*nid].species.name;
                // Initialize map entry if not yet exists
                // if !self.neib_species.contains_key(neib_sp_name) {
                //     self.neib_species
                //         .insert(neib_sp_name.to_string(), Default::default());
                // }
                self.neib_species_cur.get_mut(neib_sp_name).map(|v| *v+=1);
            }      
        }

        Ok(())
    }

    pub fn finish_frame_update(&mut self) {
        // Number of lipids have to be counted once per frame
        self.num_lip.add(self.num_lip_cur as f32);
        // Take into account the number of lipids of the current type
        for (sp,val) in self.neib_species_cur.iter() {
            let v = *val as f32 / self.num_lip_cur as f32;
            //println!("{} {} {}",self.num_lip_cur,sp,val);
            self.neib_species.get_mut(sp).unwrap().add(v);
        }
    }

    pub fn save_order_to_file(&self, dir: &Path, fname: impl AsRef<str>) -> anyhow::Result<()> {
        // We need Write traits from both fmt and io, so rewrite the imports
        // to avoid name clash
        use std::fmt::Write as FmtWrite;
        use std::io::Write as StdWrite;

        let max_len = self.order.iter().map(|v| v.x.len()).max().unwrap();

        let mstd = self
            .order
            .iter()
            .map(|mstd| mstd.compute())
            .collect::<anyhow::Result<Vec<_>>>()?;

        // Array for average order
        let mut ave_order = DVector::from_element(max_len, 0.0f32);
        // Average order is computed for each i only from the tails that have atom i
        for i in 0..max_len {
            let mut ave = 0.0;
            let mut nave = 0;
            for m in mstd.iter() {
                if let Some(v) = m.mean.get(i) {
                    ave += v;
                    nave += 1;
                }
            }
            ave_order[i] = ave / nave as f32;
        }

        // Write to file. Missing value are written as "--"
        let mut s = "# time\taver\t".to_string();
        // Headers for indivisual tails
        for t in 0..mstd.len() {
            write!(s, "tail{}", t + 1)?;
            if t < mstd.len() - 1 {
                s.push('\t');
            } else {
                s.push('\n');
            }
        }
        // Data lines
        for i in 0..max_len {
            // Carbon number and average
            write!(s, "{:.3}\t{:.3}\t", i + 1, ave_order[i])?;
            // Individual tails
            for t in 0..mstd.len() {
                let sv = mstd[t]
                    .mean
                    .get(i)
                    .map(|v| format!("{:.3}", v))
                    .unwrap_or("--".into());
                write!(s, "{sv}")?;
                if t < mstd.len() - 1 {
                    s.push('\t');
                } else {
                    s.push('\n');
                }
            }
        }
        // Write to file
        let mut f = std::fs::File::create(dir.join(fname.as_ref()))?;
        write!(f, "{}", s)?;
        Ok(())
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
        self.x += val;
        self.x2 += val * val;
        self.n += 1.0;
    }

    pub fn compute(&self) -> anyhow::Result<MeanStdResult> {
        if self.n == 0.0 {
            bail!("no values accumulated in MeanStd");
        }

        let mean = self.x / self.n;
        let x2_n = self.x2 / self.n;
        let m2 = mean * mean;
        let stddev = if x2_n > m2 {
            (x2_n - m2).sqrt()
        } else {
            0.0
        };
        Ok(MeanStdResult { mean, stddev })
    }
}

#[derive(Serialize, Debug, Clone)]
pub struct MeanStdResult {
    pub mean: f32,
    pub stddev: f32,
}

impl Serialize for MeanStd {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("MeanStd", 2)?;
        let res = self.compute().unwrap();
        state.serialize_field("mean", &res.mean)?;
        state.serialize_field("stddev", &res.stddev)?;
        state.end()
    }
}

#[derive(Debug)]
pub struct MeanStdVec {
    x: DVector<f32>,
    x2: DVector<f32>,
    n: f32,
}

#[derive(Serialize, Debug, Clone)]
pub struct MeanStdVecResult {
    mean: DVector<f32>,
    stddev: DVector<f32>,
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
                "incompatible vector size in MeanStdVec::add: {} provided, {} expected",
                val.len(),
                self.x.len()
            );
        }
        self.x = &self.x + val;
        self.x2 = &self.x2 + val.component_mul(val);
        self.n += 1.0;
        Ok(())
    }

    pub fn compute(&self) -> anyhow::Result<MeanStdVecResult> {
        if self.n == 0.0 {
            bail!("no values accumulated in MeanStd");
        }
        let mean = &self.x / self.n;
        let mut stddev = (&self.x2 / self.n) - mean.component_mul(&mean);
        stddev.apply(|v| *v = v.sqrt());
        Ok(MeanStdVecResult { mean, stddev })
    }
}

impl Serialize for MeanStdVec {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("MeanStdVec", 2)?;
        let res = self.compute().unwrap();
        state.serialize_field("mean", &res.mean)?;
        state.serialize_field("stddev", &res.stddev)?;
        state.end()
    }
}
