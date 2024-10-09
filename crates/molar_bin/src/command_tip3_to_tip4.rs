
use std::iter;

use anyhow::Result;
use log::info;
use molar::{
    core::{Atom, AtomsMutProvider, AtomsProvider, BoxProvider, PosProvider, Source},
    io::{TopologyProvider, WritableToFile},
};

pub(crate) fn command_tip3_to_tip4(
    file: &str,
    outfile: &str
) -> Result<()> {
    info!("Loading file '{file}'...");
    let mut inp = Source::serial_from_file(file)?;

    // Make output system
    let mut out = Source::empty_builder();

    // Select water
    let water = inp.select_str("resname TIP3")?;

    // For correct re-assembly of the system
    // select what is before and what is after water
    let w_first = water.first_index();
    let w_last = water.last_index();

    let sel_before = inp.select_range(&(0..w_first))?;
    let sel_after = inp.select_range(&(w_last+1..inp.num_atoms()))?;

    // Add before selection
    out.append(&sel_before);

    // Now go over water molecules one by one                   
    for mol in water.into_fragments_resindex() {
        // TIP3 is arranged as O->H->H
        // so atom 0 is O, atoms 1 and 2 are H
	    // Get cooridnates
        let o_pos = *mol.nth_pos(0).unwrap();
        let h1_pos = *mol.nth_pos(1).unwrap();
        let h2_pos = *mol.nth_pos(2).unwrap();
	    // Get center of masses of H
	    let hc = 0.5*(h1_pos.coords + h2_pos.coords);
	    // Unit vector from o to hc
	    let v = (hc-o_pos.coords).normalize();	
	    // Position of the M dummy particle in TIP4
	    let m_pos = o_pos + v*0.01546;
        // Dummy atom M
        let m_at = Atom {   
            resname: "TIP4".into(),
            name: "M".into(),
            ..mol.first_particle().atom.clone()
        };
        println!("{:?}",m_at);

        // Change resname for existing atoms
        for a in mol.iter_atoms_mut() {
            a.resname = "TIP4".into();
        }
        
        // Add new converted water molecule
        // We assume that the dummy is the last atom.
        out.append_atoms(
            mol.iter_atoms().cloned().chain(iter::once(m_at)),
            mol.iter_pos().cloned().chain(iter::once(m_pos)),
        );
    }

    // Add after selection
    out.append(&sel_after);

    // Transfer box
    out.set_box(inp.get_box().cloned());

    // Write out new system
    out.save(outfile)?;

    Ok(())
}
