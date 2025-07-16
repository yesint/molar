use anyhow::bail;
use molar::prelude::*;
use serde::Deserialize;
use sorted_vec::SortedSet;


#[derive(Clone, Debug, Deserialize)]
pub struct LipidSpeciesDescr {
    pub(super) whole: String,
    pub(super) head: String,
    pub(super) mid: String,
    pub(super) tails: Vec<String>,
    #[serde(default)]
    pub(super) max_area: f32,
}

#[derive(Debug)]
pub struct LipidSpecies {
    pub(super) name: String,
    //pub(super) descr: LipidSpeciesDescr,
    pub(super) head_marker_offsets: SortedSet<usize>,
    pub(super) mid_marker_offsets: SortedSet<usize>,
    pub(super) tails: Vec<LipidTail>,
    pub(super) max_area: f32,
}

#[derive(Debug)]
pub struct LipidTail {
    //pub(super) descr: String,
    pub(super) offsets: Vec<usize>,
    pub(super) bond_orders: Vec<u8>,
}

impl LipidSpecies {
    pub fn new<K: UserCreatableKind>(
        name: String,
        descr: LipidSpeciesDescr,
        lipid: &Sel<K>,
    ) -> anyhow::Result<Self> {
        let first_index = lipid.first_index();
        let mut tails = vec![];

        for t in &descr.tails {
            let mut names = vec![];
            let mut bond_orders = vec![];
            let mut cur = &t[..];

            while let Some(e) = cur.find(['-', '=']) {
                if cur[..e].is_empty() {
                    bail!("missing carbon atom name");
                }
                names.push(&cur[..e]);
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

            // Convert names to offsets
            let offsets = names.iter()
                .map(|name| {
                    let atom = lipid.subsel(format!("name {name}"))?;
                    if atom.len() > 1 {
                        bail!("more than one tail atom {name} in lipid");
                    }
                    Ok(atom.first_index() - first_index)
                })
                .collect::<anyhow::Result<Vec<_>>>()?;

            tails.push(LipidTail {
                //descr: t.to_owned(),
                offsets,
                bond_orders,
            });
        }

        Ok(Self {
            name,
            //descr: descr.clone(),
            head_marker_offsets: SortedSet::from_unsorted(
                lipid.subsel(descr.head)?
                    .iter_index()
                    .map(|i| i - first_index)
                    .collect(),
            ),
            mid_marker_offsets: SortedSet::from_unsorted(
                lipid.subsel(descr.mid)?
                    .iter_index()
                    .map(|i| i - first_index)
                    .collect(),
            ),
            tails,
            max_area: descr.max_area,
        })
    }
}
