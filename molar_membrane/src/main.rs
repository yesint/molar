use anyhow::{Context, Result};
use clap::Args;
use molar::prelude::*;
use molar_membrane::{Histogram1D, Membrane};
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
    all_hist: Histogram1D,
}

impl AnalysisTask<Flags> for MembraneBilayerTask {
    fn task_name() -> String {
        "Bilayer analysis".to_owned()
    }

    fn new(context: &AnalysisContext<Flags>) -> anyhow::Result<Self> {
        let mut toml = String::new();
        std::fs::File::open(&context.args.params_file)
            .context(format!(
                "reading membrane options file {}",
                &context.args.params_file
            ))?
            .read_to_string(&mut toml)?;
        
        let membr = Membrane::new(&context.src, &toml)?;
        
        Ok(Self { 
            membr,
            all_hist: Histogram1D::new(-0.15, 0.15, 100),
        })
    }

    fn process_frame(&mut self, context: &AnalysisContext<Flags>) -> anyhow::Result<()> {
        self.membr.set_state_from(&context.src)?;

        // if context.consumed_frames == 500 {
        //     context.src.save("0.gro")?;
        // }
        
        self.membr.reset_groups();
        self.membr.reset_valid_lipids();

        let mut all = vec![];
        for lip in self.membr.iter_valid_lipids() {
            let x = lip.head_marker.x;
            if x > 5.5 && x < 12.5 {
                all.push(lip.id);
            }
        }

        self.membr.add_ids_to_group("all", &all)?;

        self.membr.compute()?;

        // if context.consumed_frames == 500 {
        //     println!("vis");
        //     self.membr.write_vmd_visualization("vis.tcl")?;
        // }

        for lip in self.membr.iter_group_valid("all")? {
            self.all_hist.add_one(lip.mean_curv);
        }

        Ok(())
    }

    fn post_process(&mut self, _context: &AnalysisContext<Flags>) -> anyhow::Result<()> {
        self.membr.finalize()?;
        self.all_hist.normalize_density();
        self.all_hist.save_to_file("results/hist.dat")?;
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
