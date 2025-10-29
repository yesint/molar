use std::collections::HashSet;

use anyhow::{bail, Context, Result};
use log::info;
use molar::prelude::*;

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

    if !begin.is_empty() {
        info!("Selections to put in the beginning:");
        for s in begin.iter() {
            info!("\t{s}");
        }
    }

    if !end.is_empty() {
        info!("Selections to put at the end:");
        for s in end.iter() {
            info!("\t{s}");
        }
    }

    // Make selections
    info!("Rearranging file '{infile}'...");
    let in_source = System::from_file(infile)?;
    let begin_sels = begin
        .iter()
        .map(|s| in_source.select(s))
        .collect::<Result<Vec<Sel>,_>>()
        .with_context(|| "can't create begin selections for rearranging")?;

    let end_sels = end
        .iter()
        .map(|s| in_source.select(s))
        .collect::<Result<Vec<Sel>,_>>()
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
    let all_ind = (0..in_source.len()).collect::<HashSet<usize>>();

    let rest_sel = in_source
        .select(all_ind.difference(&used).cloned().collect::<Vec<_>>())
        .ok();

    // Create output builder
    let mut out = System::new_empty();
    // Add beginning selections
    for sel in begin_sels {
        out.append(&sel);
    }

    // Add the rest if any
    if let Some(rest) = rest_sel {
        info!("There are {} untouched atoms", rest.len());
        out.append(&rest);
    }

    // Add ending selections
    for sel in end_sels {
        out.append(&sel);
    }

    info!("Writing rearranged to '{outfile}'...");
    // Write rearranged
    out.bind()?.save(outfile)?;

    Ok(())
}
