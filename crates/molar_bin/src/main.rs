use anyhow::{anyhow, bail, Result};
use clap::{Parser, Subcommand};
use log::info;
use molar::prelude::*;

/// MolAR binary utility
#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cmd {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Extracts last frame from the trajectory and writes it to file
    Last {
        /// Two files for structure and trajectory in any order,
        /// or a single file with both structure and trajectory
        #[arg(short, num_args=1..=2)]
        files: Vec<String>,
        /// Output file (structure+trajectory)
        #[arg(short, default_value = "last.gro")]
        outfile: String,
    },
}

fn main() -> Result<()> {
    env_logger::builder()
        .format_timestamp(None).format_indent(Some(8))
        .init();

    let cmd = Cmd::parse();
    match &cmd.command {
        Commands::Last { files, outfile } => {
            let mut trj; // Trajectory file
                         // Reading topology
            let top = if files.len() == 1 {
                trj = FileHandler::open(&files[0])?;
                let top = trj.read_topology()?;
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

            let (last_fr, last_time) = trj.tell_last()?;
            trj.seek_frame(last_fr)?;

            let st = trj
                .read_state()?
                .ok_or_else(|| anyhow!("Frame {last_fr} can't be read"))?;

            info!("Writing last frame #{last_fr}, time {last_time} to '{outfile}'...");
            let mut out = FileHandler::create(&outfile)?;
            out.write(&(top, st))?;
        }
    }
    Ok(())
}
