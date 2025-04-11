use std::{io::Read, path::PathBuf};
use clap::{Args, CommandFactory, FromArgMatches, Parser};
use molar::prelude::*;
use anyhow::Result;
use molar_membrane::Membrane;

#[derive(Args, Debug, Default)]
struct Flags {
    #[arg(short, long, required=true)]
    lipids_file: String,
    
    #[arg(short, long, default_value_t=2.5)]
    cutoff: f32,
    
    #[arg(long, default_value_t=1)]
    max_iter: usize,
    
    #[arg(long, default_value_t=false)]
    use_scd_corr: bool,
    
    #[arg(short,long, default_value=".")]
    output_dir: PathBuf,

    #[arg(short,long, default_value="all")]
    sel_center: String,
}

struct MembraneBilayerTask {
    
}


impl AnalysisTask for MembraneBilayerTask {
    type Options = Flags;

    fn pre_process(&mut self, context: &AnalysisContext<Flags>) -> anyhow::Result<()> {
        Ok(())
    }

    fn process_frame(&mut self, context: &AnalysisContext<Flags>) -> anyhow::Result<()> {
        Ok(())
    }

    fn post_process(&mut self, context: &AnalysisContext<Flags>) -> anyhow::Result<()> {
        Ok(())
    }
}

fn main() -> Result<()> {
    env_logger::builder()
        .format_timestamp(None)
        .format_indent(Some(8))
        .filter_level(log::LevelFilter::Info)
        .init();
    
    // Greeting
    //molar::greeting("molar_bin");

    let mut task = MembraneBilayerTask{};
    task.run()?;
    Ok(())
}