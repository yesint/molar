use std::iter;

use anyhow::Result;
use log::info;
use molar::prelude::*;

pub(crate) fn command_tip3_to_tip4(file: &str, outfile: &str) -> Result<()> {
    info!("Loading file '{file}'...");
    let inp = System::from_file(file)?;

    // Make output system
    let mut out = System::default();

    // Select water
    let water = inp.select("resname TIP3")?;

    // For correct re-assembly of the system
    // select what is before and what is after water
    let w_first = water.bind(&inp).first_index();
    let w_last = water.bind(&inp).last_index();

    let sel_before = inp.select(0..w_first)?;
    let sel_after = inp.select(w_last + 1..inp.len())?;

    // Add before selection
    out.append(&sel_before.bind(&inp))?;

    // Now go over water molecules one by one
    for mol in water.bind(&inp).split_resindex_iter() {
        let mol = mol.bind(&inp);
        // TIP3 is arranged as O->H->H
        // so atom 0 is O, atoms 1 and 2 are H
        // Get cooridnates
        let o_pos = *mol.get_pos(0).unwrap();
        let h1_pos = *mol.get_pos(1).unwrap();
        let h2_pos = *mol.get_pos(2).unwrap();
        // Get center of masses of H
        let hc = 0.5 * (h1_pos.coords + h2_pos.coords);
        // Unit vector from o to hc
        let v = (hc - o_pos.coords).normalize();
        // Position of the M dummy particle in TIP4
        let m_pos = o_pos + v * 0.01546;
        // Dummy atom M
        let m_at = Atom {
            name: "M".into(),
            ..mol.first_particle().atom.clone()
        };
        println!("{:?}", m_at);

        // Add new converted water molecule
        // We assume that the dummy is the last atom.
        let added = out.append_atoms_pos(
            mol.iter_atoms().chain(iter::once(&m_at)),
            mol.iter_pos().chain(iter::once(&m_pos)),
        )?;
        
        // Change resname for added atoms
        added.bind_mut(&mut out).set_same_resname("TIP4");
    }

    // Add after selection
    out.append(&sel_after.bind(&inp))?;

    // Transfer box
    out.set_box_from(&inp);

    // Write out new system
    out.save(outfile)?;

    Ok(())
}
