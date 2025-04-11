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
    membr: Membrane,
}


impl AnalysisTask<Flags> for MembraneBilayerTask {

    fn new(context: &AnalysisContext<Flags>) -> anyhow::Result<Self> {
        let z0 = context.src.select(&context.args.sel_center)?.center_of_mass()?.z;
        let mut toml = String::new();
        std::fs::File::open(&context.args.lipids_file)?.read_to_string(&mut toml)?;
        let mut membr = Membrane::new(&context.src, &toml)?
            .with_global_normal(Vector3f::z())
            .with_order_type(OrderType::ScdCorr)
            .with_output_dir(&context.args.output_dir)?
            .with_groups(&["upper", "lower"])
            .with_cutoff(context.args.cutoff)
            .with_max_iter(context.args.max_iter);let mut upper = vec![];
        
        let mut lower = vec![];
        for (id, lip) in membr.iter_lipids() {
            if lip.head_marker.z > z0 {
                upper.push(id);
            } else {
                lower.push(id);
            }
        }

        membr.add_lipids_to_group("upper", upper)?;
        membr.add_lipids_to_group("lower", lower)?;

        Ok(Self{membr})
    }

    fn process_frame(&mut self, context: &AnalysisContext<Flags>) -> anyhow::Result<()> {
        self.membr.set_state(context.src.get_state().clone())?;
        self.membr.compute()?;
        Ok(())
    }

    fn post_process(&mut self, _context: &AnalysisContext<Flags>) -> anyhow::Result<()> {
        self.membr.finalize()?;
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

    MembraneBilayerTask::run()?;
    Ok(())
}