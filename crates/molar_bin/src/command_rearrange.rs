use std::collections::HashSet;

use anyhow::{anyhow, bail, Context, Result};
use log::info;
use molar::{
    core::{MutableSerial, Sel, Source},
    io::{FileHandler, IndexProvider, TopologyProvider, WritableToFile},
};

pub(super) fn command_rearrange(
    infile: &str,
    outfile: &str,
    begin: &Vec<String>,
    end: &Vec<String>,
) -> Result<()> {
    // Sanity check
    if begin.is_empty() && end.is_empty() {
        bail!("Provide at least one selection for rearranging!");
    }

    info!("Rearranging atoms from file '{infile}'...");
    let mut data = FileHandler::open(infile)?;
    let (top, st) = data.read()?;

    if !begin.is_empty() {
        info!("Selections to put in the beginning:");
    }
    for s in begin.iter() {
        info!("\t{s}");
    }

    if !end.is_empty() {
        info!("Selections to put at the end:");
    }
    for s in end.iter() {
        info!("\t{s}");
    }

    // Make selections
    let in_source = Source::new(top, st)?;
    let begin_sels = begin
        .iter()
        .map(|s| in_source.select_str(s).map_err(|e| anyhow!(e)))
        .collect::<Result<Vec<Sel<MutableSerial>>>>()
        .with_context(|| "can't create begin selections for rearranging")?;

    let end_sels = end
        .iter()
        .map(|s| in_source.select_str(s).map_err(|e| anyhow!(e)))
        .collect::<Result<Vec<Sel<MutableSerial>>>>()
        .with_context(|| "can't create end selections for rearranging")?;

    // Check overlap of selections
    let mut used = HashSet::<usize>::default();
    for sel in begin_sels.iter().chain(end_sels.iter()) {
        for i in sel.iter_index() {
            if !used.insert(i) {
                bail!("selections for rearrangement overlap at atom {i}");
            }
        }
    }

    // Get the rest of indexes, which are not used
    let all_ind = (0..in_source.num_atoms()).collect::<HashSet<usize>>();
    let rest_sel = in_source.select_from_iter(all_ind.difference(&used).cloned())?;

    info!("There are {} untouched atoms",rest_sel.len());

    // Create output builder
    let mut out = Source::empty_builder()?;
    // Add beginning selections
    for sel in begin_sels {
        out.append(&sel);
    }
    // Add the rest
    out.append(&rest_sel);
    // Add ending selections
    for sel in end_sels {
        out.append(&sel);
    }

    info!("Writing rearranged to '{outfile}'...");
    // Write rearranged
    out.save(outfile)?;

    Ok(())
}