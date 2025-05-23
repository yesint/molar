use anyhow::{Context, Result};
use clap::Args;
use molar::prelude::*;
use molar_membrane::Membrane;
use std::io::Read;

#[derive(Args, Debug, Default)]
struct Flags {
    #[arg(short, long, required = true)]
    params_file: String,

    #[arg(short, long, default_value = "all")]
    sel_center: String,
}

struct MembraneBilayerTask {
    membr: Membrane,
}

impl AnalysisTask<Flags> for MembraneBilayerTask {
    fn task_name() -> String {
        "Bilayer analysis".to_owned()
    }

    fn new(context: &AnalysisContext<Flags>) -> anyhow::Result<Self> {
        let z0 = context
            .src
            .select(&context.args.sel_center)?
            .center_of_mass()?
            .z;
        
        let mut toml = String::new();
        std::fs::File::open(&context.args.params_file)
            .context(format!(
                "reading membrane options file {}",
                &context.args.params_file
            ))?
            .read_to_string(&mut toml)?;
        
        let mut membr = Membrane::new(&context.src, &toml)?;

        let mut upper = vec![];
        let mut lower = vec![];
        for lip in membr.iter_valid_lipids() {
            if lip.head_marker.z > z0 {
                upper.push(lip.id);
            } else {
                lower.push(lip.id);
            }
        }

        membr.add_lipids_to_group("upper", &upper)?;
        membr.add_lipids_to_group("lower", &lower)?;

        Ok(Self { membr })
    }

    fn process_frame(&mut self, context: &AnalysisContext<Flags>) -> anyhow::Result<()> {
        // We need to update the state
        self.membr.set_state(context.src.get_state().release()?)?;
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

    MembraneBilayerTask::run()?;
    Ok(())
}
