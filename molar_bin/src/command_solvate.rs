use anyhow::{Context, Result};
use log::info;
use molar::prelude::*;

pub(crate) fn command_solvate(
    file: &str,
    outfile: &str,
    solvent: &Option<String>,
    exclude: &Option<String>,
) -> Result<()> {
    info!("Loading solute from file '{file}'...");
    let mut solute = System::from_file(file)?;

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
    let solvent = System::from_file(&solvent_file)?;

    let n = solute.len();
    let added = solute.solvate(&solvent)?;
    info!("{added} solvent molecules ({} atoms) added", solute.len() - n);

    // If exclude selection is provided remove it
    if let Some(sel_str) = exclude {
        let not_excl_sel = solute.select_bound(format!("not ({})", sel_str))?;
        info!("Excluding atoms by selection '{}'", sel_str);
        not_excl_sel.save(outfile)?;
    } else {
        solute.save(outfile)?;
    }

    Ok(())
}
