use std::sync::Arc;

use anyhow::bail;
use molar::prelude::*;
use sorted_vec::SortedSet;

#[derive(Clone, Debug)]
pub struct LipidSpeciesDescr {
    pub name: String,
    pub whole_sel_str: String,
    pub head_marker_subsel_str: String,
    pub mid_marker_subsel_str: String,
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
                cur = &cur[e+1..];
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
                    .subsel_str(descr.head_marker_subsel_str)?
                    .iter_index()
                    .map(|i| i - first_index)
                    .collect(),
            ),
            mid_marker_offsets: SortedSet::from_unsorted(
                lipid
                    .subsel_str(descr.mid_marker_subsel_str)?
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
    species: Arc<LipidSpecies>,
    head_marker: Pos,
    mid_marker: Pos,
    tail_marker: Pos,
    tails: Vec<LipidTail>,
}

impl LipidMolecule {
    pub fn update<K: UserCreatableKind>(&mut self, sel: &Sel<K>) -> anyhow::Result<()> {
        self.head_marker = sel
            .subsel_iter(self.species.head_marker_offsets.iter().cloned())?
            .center_of_mass_pbc()?;
        self.mid_marker = sel
            .subsel_iter(self.species.mid_marker_offsets.iter().cloned())?
            .center_of_mass_pbc()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use molar::core::Source;

    use crate::{LipidSpecies, LipidSpeciesDescr};

    #[test]
    fn test_descr() -> anyhow::Result<()> {
        let src = Source::serial_from_file("../../tests/membr.gro")?;
        let pope = src.select_str("resid 60")?;
        let descr = LipidSpeciesDescr {
            name: "POPE".into(),
            whole_sel_str: "resname POPE".into(),
            head_marker_subsel_str: "name P N".into(),
            mid_marker_subsel_str: "name C21 C22".into(),
            tails: vec![
                "C21-C22-C23-C24-C25-C26-C27-C28-C29=C210-C211-C212-C213-C214-C215-C216-C217-C218".into(),
                "C31-C32-C33-C34-C35-C36-C37-C38-C39=C310-C311-C312-C313-C314-C315-C316".into(),
            ],
        };
        let lip_sp = LipidSpecies::new(descr, &pope)?;
        println!("{lip_sp:?}");
        Ok(())
    }
}