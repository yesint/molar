use std::collections::HashSet;

use anyhow::{bail, Context, Result};
use log::info;
use molar::prelude::*;

pub(crate) fn command_solvate(
    file: &str,
    outfile: &str,
    solvent: &Option<String>,
    exclude: &Option<String>,
) -> Result<()> {
    info!("Loading solute from file '{file}'...");
    let solute = System::from_file(file)?;
    let solute_g = solute.bind()?;

    // Check periodic box
    if solute_g.get_box().is_none() {
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
    let solvent = System::from_file(&solvent_file)?;
    let solvent_g = solvent.bind()?;
    if solvent_g.get_box().is_none() {
        bail!("solvent lacks a periodic box");
    }
    if solvent_g.get_box().unwrap().is_triclinic() {
        bail!("triclinic solvent boxes are not supported yet");
    }

    // Get solute extents
    let solute_max_ext = solute_g.get_box().unwrap().get_lab_extents();

    // We will fill the rectangular region with solvent
    // even if the actual box is triclinic and then we will
    // remove molecules outside the box
    let solvent_ext = solvent_g.get_box().unwrap().get_box_extents();
    let mut nbox = [0; 3];
    for i in 0..=2 {
        nbox[i] = (solute_max_ext[i] / solvent_ext[i]).ceil() as usize;
    }
    info!(
        "Solute box extents [{},{},{}]",
        solute_max_ext[0], solute_max_ext[1], solute_max_ext[2]
    );
    info!(
        "Solvent box extents [{},{},{}]",
        solvent_ext[0], solvent_ext[1], solvent_ext[2]
    );
    info!(
        "Duplicating solvent box [{},{},{}]...",
        nbox[0], nbox[1], nbox[2]
    );

    // Duplicating the solvent
    solvent.multiply_periodically(nbox)?;

    info!("Removing solvent molecules outside the box...");
    //solvent.save("./target/duplicated_solvent.gro")?;

    // Removing solvent molecules not in the box
    let mut inside_ind = vec![];
    let b = solute.get_box().unwrap();
    let all = solvent.select_all()?;
    'outer: for res in all.split_resindex_iter() {
        for p in res.iter_pos() {
            if !b.is_inside(p) {
                // Break without adding this residue to good list
                continue 'outer;
            }
        }
        // If we are here then all atoms of this residue are inside the box
        inside_ind.extend(res.iter_index());
    }
    // It is safe to call here since indexes are guaranteed to be properly sorted
    let inside_sel = solvent.select(inside_ind)?;

    //inside_sel.save("target/inside.gro")?;

    info!("Searching for overlapping solvent molecules...");
    // Do the distance search
    let vdw1 = inside_sel.iter_atoms().map(|a| a.vdw()).collect();
    let vdw2 = solute.iter_atoms().map(|a| a.vdw()).collect();
    let local_overlap_ind: Vec<usize> = distance_search_double_vdw_pbc(
        inside_sel.iter_pos(),
        solute.iter_pos(),
        &vdw1,
        &vdw2,
        b,
        PBC_FULL,
    );
    //distance_search_double_pbc(0.3,&inside_sel, &solute, b, PBC_FULL);
    info!("{} overlapping atoms", local_overlap_ind.len());

    // Find all resindexes for atoms to be removed
    let mut resind_to_remove = HashSet::new();
    for i in local_overlap_ind {
        unsafe {
            resind_to_remove.insert(inside_sel.get_atom_unchecked(i).resindex);
        }
    }
    info!(
        "{} overlapping water molecules to remove",
        resind_to_remove.len()
    );

    let mut good_ind = vec![];
    for res in inside_sel.split_resindex_iter() {
        if !resind_to_remove.contains(&res.first_atom().resindex) {
            good_ind.extend(res.iter_index());
        }
    }
    info!("{} atoms to add", good_ind.len());

    let good_sel = solvent.select(good_ind)?;

    // Add solvent
    solute.append(&good_sel);

    // If exclude selection is provided remove it
    if exclude.is_some() {
        let sel_str = exclude.as_ref().unwrap();
        let not_excl_sel = solute.select(format!("not ({})", sel_str))?;
        info!("Excluding atoms by selection '{}'", sel_str);
        not_excl_sel.save(outfile)?;
    } else {
        solute.save(outfile)?;
    }

    Ok(())
}
