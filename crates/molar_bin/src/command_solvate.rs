use std::path::PathBuf;

use crate::SolvateMode;
use anyhow::{anyhow, bail, Context, Result};
use log::info;
use molar::{
    core::{MutableSerial, Sel, Source}, distance_search, io::{FileHandler, IndexProvider, TopologyProvider, WritableToFile}, prelude::DistanceSearcherDouble
};

pub(crate) fn command_solvate(
    file: &str,
    outfile: &str,
    solvent: &Option<String>,
    exclude: &Option<String>,
    mode: &SolvateMode,
) -> Result<()> {
    info!("Loading solute from file '{file}'...");
    let mut data = FileHandler::open(file)?;
    let (top, st) = data.read()?;

    // Check periodic box
    if st.get_box().is_none() {
        bail!("can't solvate a system without a periodic box");
    }

    // Solvent
    let solvent_file = match solvent {
        Some(s) => s.clone(),
        None => {
            let gmx_dir = std::env::var("GMXDATA").with_context(|| {
                "GMXDATA env variable is not set, use explicit path to the solvent file"
            })?;
            format!("{gmx_dir}/top/spc216.gro")
        }
    };
    info!("Loading solvent from file '{solvent_file}'...");

    //let searcher = DistanceSearcherDouble::from_sel_vdw(sel1, sel2);

    Ok(())
}
