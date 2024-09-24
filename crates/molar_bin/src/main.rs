use anyhow::Result;
use clap::{Parser, Subcommand, ValueEnum};

mod command_last;
mod command_rearrange;
mod command_solvate;
mod command_tip3_to_tip4;

use command_last::command_last;
use command_rearrange::command_rearrange;
use command_solvate::command_solvate;
use command_tip3_to_tip4::command_tip3_to_tip4;

/// MolAR binary utility
#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cmd {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum SolvateMode {
    /// Plain distance cutoff
    Distance,
    /// Use VdW radii
    Vdw,
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

    /// Rearranges the atoms in a specified order
    Rearrange {
        /// Input file containing structure and coordinates
        #[arg(short, required = true)]
        file: String,
        /// Output file (structure+coordinates)
        #[arg(short, default_value = "rearranged.gro")]
        outfile: String,
        /// Selections to be put in the beginning
        #[arg(short, long)]
        begin: Vec<String>,
        /// Selections to be put at the end
        #[arg(short, long)]
        end: Vec<String>,
    },

    /// Performs solvation
    Solvate {
        /// Input structure to be solvated
        /// Must contain structure and coordinates
        #[arg(short, required = true)]
        file: String,

        /// Output file (structure+coordinates)
        #[arg(short, default_value = "rearranged.gro")]
        outfile: String,

        /// File with a periodic solvent box
        /// If not provided tries to use Gromacs SPC water box
        #[arg(short, long)]
        solvent: Option<String>,

        /// Selection for solvent molecules to be deleted
        #[arg(long)]
        exclude: Option<String>,

        /// Distance search mode
        #[arg(long,value_enum,default_value_t=SolvateMode::Distance)]
        mode: SolvateMode,
    },

    /// Converts TIP3 water to TIP4 water
    Tip3to4 {
        /// Input file (structure+coordinates)
        #[arg(short, required = true)]
        file: String,

        /// Output file (structure+coordinates)
        #[arg(short, default_value = "rearranged.gro")]
        outfile: String,
    }
}

fn main() -> Result<()> {
    env_logger::builder()
        .format_timestamp(None)
        .format_indent(Some(8))
        .filter_level(log::LevelFilter::Info)
        .init();

    let cmd = Cmd::parse();

    // Greeting
    greeting();

    match &cmd.command {
        Commands::Last { files, outfile } => {
            println!("▶ Action: last");
            command_last(files, outfile)?;
        },

        Commands::Rearrange {
            file: infile,
            outfile,
            begin,
            end,
        } => {
            println!("▶ Action: rearrange");
            command_rearrange(infile, outfile, begin, end)?;
        },

        Commands::Solvate {
            file,
            outfile,
            solvent,
            exclude,
            mode
        } => {
            println!("▶ Action: solvate");
            command_solvate(file, outfile, solvent, exclude ,mode)?;
        },

        Commands::Tip3to4 { file, outfile } => {
            println!("▶ Action: Tip3toTip4");
            command_tip3_to_tip4(file, outfile)?;
        },
    }
    Ok(())
}

fn greeting() {
    use comfy_table::modifiers::UTF8_ROUND_CORNERS;
    use comfy_table::presets::UTF8_FULL;
    use comfy_table::{Attribute, Cell, Table};

    let mut table = Table::new();
    table
        .load_preset(UTF8_FULL)
        .apply_modifier(UTF8_ROUND_CORNERS)
        .add_row(vec![
            Cell::new("MolAR - Molecular Analysis for Rust").add_attributes(vec![Attribute::Bold])
        ])
        .add_row(vec![format!(
            "{}\n{}",
            env!("CARGO_PKG_HOMEPAGE"),
            env!("CARGO_PKG_AUTHORS")
        )])
        .add_row(vec![format!("MolAR version: {}", molar::VERSION)])
        .add_row(vec![format!(
            "Utility: {}, Version: {}",
            env!("CARGO_PKG_NAME"),
            env!("CARGO_PKG_VERSION")
        )]);
    println!("{table}");
}
