use anyhow::{anyhow, Context, Result};
use log::info;
use molar::prelude::*;

pub(super) fn command_last(files: &Vec<String>, outfile: &str) -> Result<()> {
    let mut trj; // Trajectory file
    // Reading topology
    let top = if files.len() == 1 {
        trj = FileHandler::open(&files[0])?;
        let top = trj
            .read_topology()
            .with_context(|| "loading both topology and state")?;
        info!("Topology loaded from {}", files[0]);
        info!("Trajectory file {}", files[0]);
        top
    } else {
        let mut fh = FileHandler::open(&files[0])?;
        if let Ok(top) = fh.read_topology() {
            // File 1 is topology, file 2 is trajectory
            trj = FileHandler::open(&files[1])?;
            info!("Topology loaded from {}", files[0]);
            info!("Trajectory file {}", files[1]);
            top
        } else {
            // File 2 is topology, file 1 is trajectory
            trj = fh;
            let top = FileHandler::open(&files[1])?.read_topology()?;
            info!("Topology loaded from {}", files[1]);
            info!("Trajectory file {}", files[0]);
            top
        }
    };

    let mut slow_forward = false;

    if let Err(_) = trj.seek_last() {slow_forward = true};

    let st = if slow_forward {
        info!("Fast-forward is not possible, reading the whole trajectory...");
        trj.into_iter().last()
    } else {
        info!("Fast-forward to the last frame...");
        Some(trj.read_state()?)
    }
    .ok_or_else(|| anyhow!("Last frame can't be read"))?;

    info!(
        "Writing last frame time={} to '{}'...",
        st.get_time(), outfile
    );
    let mut out = FileHandler::create(&outfile)?;
    out.write(&System::new(top, st)?)?;

    Ok(())
}
