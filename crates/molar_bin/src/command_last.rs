use anyhow::{anyhow, Context, Result};
use log::info;
use molar::io::{FileHandler, StateProvider};

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
    let mut last_fr = 0;

    match trj.tell_last() {
        Ok((fr, t)) => {
            info!("Fast-forward to last frame {last_fr}...");
            if let Err(_) = trj.seek_frame(fr) {
                slow_forward = true;
            } else {
                last_fr = fr;
            }
        }
        Err(_) => {
            slow_forward = true;
        }
    }

    let st = if slow_forward {
        info!("Fast-forward is not possible, reading the whole trajectory...");
        trj.into_iter().last()
    } else {
        trj.read_state()?
    }
    .ok_or_else(|| anyhow!("Last frame can't be read"))?;

    info!(
        "Writing last frame #{last_fr}, time {} to '{outfile}'...",
        st.get_time()
    );
    let mut out = FileHandler::create(&outfile)?;
    out.write(&(top, st))?;

    Ok(())
}
