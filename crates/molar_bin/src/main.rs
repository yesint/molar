use anyhow::{anyhow, bail, Result};
use clap::{Parser, Subcommand};
use log::info;
use molar::prelude::*;

mod command_last;
use command_last::command_last;

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
        #[arg(short, num_args=1..=2, required=true)]
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
            command_last(files, outfile)?;
        }
    }
    Ok(())
}
