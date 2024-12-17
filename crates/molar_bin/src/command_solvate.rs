use anyhow::{bail, Context, Result};
use log::{debug, info};
use molar::prelude::*;

pub(crate) fn command_solvate(
    file: &str,
    outfile: &str,
    solvent: &Option<String>,
    exclude: &Option<String>,
) -> Result<()> {
    info!("Loading solute from file '{file}'...");
    let solute = Source::builder_from_file(file)?;

    // Check periodic box
    if solute.get_box().is_none() {
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
    let solvent = Source::builder_from_file(&solvent_file)?;
    if solvent.get_box().is_none() {
        bail!("solvent lacks a periodic box");
    }
    if solvent.get_box().unwrap().is_triclinic() {
        bail!("triclinic solvent boxes are not supported yet");
    }

    // Get solute extents
    let solute_max_ext = solute.get_box().unwrap().get_lab_extents();

    // We will fill the rectangular region with solvent
    // even if the actual box is triclinic and then we will
    // remove molecules outside the box
    let solvent_ext = solvent.get_box().unwrap().get_box_extents();
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
    //solvent.multiply_periodically(nbox)?;

    let solvent_sel = solvent.select_all()?;
    for x in 0..=nbox[0] {
        for y in 0..=nbox[1] {
            for z in 0..=nbox[2] {
                if x == 0 && y == 0 && z == 0 {
                    continue;
                }
                let added = solvent.append(&solvent_sel);
                let shift = Vector3f::new(
                    solvent_ext[0] * x as f32,
                    solvent_ext[1] * y as f32,
                    solvent_ext[2] * z as f32,
                );
                added.translate(&shift);
            }
        }
    }

    info!("Removing solvent molecules outside the box...");
    //solvent.save("./target/duplicated_solvent.gro")?;

    // Removing solvent molecules not in the box
    let mut inside_ind = vec![];
    let b = solute.get_box().unwrap();
    let all = solvent.select_all()?;
    'outer: for res in all.into_iter_fragments_resindex() {
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
    let inside_sel = unsafe { solvent.select_vec_unchecked(inside_ind)? };

    //inside_sel.save("target/inside.gro")?;

    info!("Removing overlapping solvent molecules...");
    // Do the distance search
    let vdw1 = inside_sel.iter_atoms().map(|a| a.vdw()).collect();
    let vdw2 = solute.iter_atoms().map(|a| a.vdw()).collect();
    let to_remove_ind: Vec<usize> =
        distance_search_double_vdw_pbc(&inside_sel, &solute, &vdw1, &vdw2, b, PBC_FULL);
        //distance_search_double_pbc(0.3,&inside_sel, &solute, b, PBC_FULL);

    // Remove overlapping
    solvent.remove(&to_remove_ind)?;
    // Add solvent
    solute.append(&solvent);
    // If exclude selection is provided remove it
    if exclude.is_some() {
        let sel_str = exclude.as_ref().unwrap();
        let excl_sel = solute.select_str(sel_str)?;
        info!(
            "Excluding {} atoms by selection '{}'",
            excl_sel.len(),
            sel_str
        );
        solute.remove(&excl_sel)?;
    }

    // Save
    solute.save(outfile)?;

    Ok(())
}
