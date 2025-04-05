use anyhow::bail;
use log::info;
use nalgebra::DVector;
use serde::{ser::SerializeStruct, Serialize, Serializer};
use std::collections::HashMap;
use std::path::Path;

use crate::lipid_species::LipidSpecies;
use crate::LipidMolecule;

#[derive(Default, Debug)]
pub struct GroupProperties {
    pub per_species: HashMap<String, StatProperties>,
}

impl GroupProperties {
    pub fn add_lipid_stats<'a>(
        &'a mut self,
        gr_lipids: impl Iterator<Item = &'a LipidMolecule>,
    ) -> anyhow::Result<()> {
        for lip in gr_lipids {
            let sp_name = &lip.species.name;
            self.per_species
                .get_mut(sp_name)
                .unwrap()
                .add_single_lipid_stats(lip)?;
        }
        Ok(())
    }

    pub fn save_order_to_file(
        &self,
        dir: &Path,
        gr_name: impl AsRef<str>,
    ) -> anyhow::Result<()> {
        for (sp, stat) in &self.per_species {
            info!("\tWriting order for '{sp}'");
            stat.save_order_to_file(
                dir,
                format!("gr_{}_order_{}.dat", gr_name.as_ref(), sp),
            )?;
        }
        Ok(())
    }
}

#[derive(Default, Debug)]
pub struct StatProperties {
    pub num_lip: usize,
    pub area: MeanStd,
    pub tilt: MeanStd,
    pub order: Vec<MeanStdVec>,
}

impl StatProperties {
    pub fn new(species: &LipidSpecies) -> Self {
        let mut order = Vec::with_capacity(species.tails.len());
        for t in &species.tails {
            order.push(MeanStdVec::new(t.bond_orders.len() - 1));
        }
        Self {
            area: Default::default(),
            tilt: Default::default(),
            order,
            num_lip: 0,
        }
    }

    pub fn add_single_lipid_stats(&mut self, lip: &LipidMolecule) -> anyhow::Result<()> {
        // self.area.add(lip.props.area);
        // self.tilt.add(lip.props.tilt);
        for tail in 0..lip.order.len() {
            self.order[tail].add(&lip.order[tail])?;
        }
        Ok(())
    }

    pub fn save_order_to_file(
        &self,
        dir: &Path,
        fname: impl AsRef<str>,
    ) -> anyhow::Result<()> {
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
        let stddev = ((self.x2 / self.n) - (mean * mean)).sqrt();
        Ok(MeanStdResult { mean, stddev })
    }
}

#[derive(Serialize, Debug, Clone)]
pub struct MeanStdResult {
    mean: f32,
    stddev: f32,
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
