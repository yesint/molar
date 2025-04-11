use clap::{CommandFactory, FromArgMatches, Parser, Args};
use log::info;
use std::num::{ParseFloatError, ParseIntError};
use std::path::PathBuf;
use thiserror::Error;

use crate::core::{
    Holder, MutableSerial, RandomAtomProvider, SelectionError, Source, State, StateProvider,
    Topology,
};
use crate::io::{FileHandler, FileIoError};

#[derive(Parser)]
#[command(name = "analysis")]
pub struct AnalysisArgs {
    /// Input files (topology and trajectory)
    #[clap(short = 'f', long = "files", required = true, num_args = 2..)]
    pub files: Vec<PathBuf>,

    /// Logging frequency
    #[clap(long = "log", default_value = "100")]
    pub log: usize,

    /// Begin time/frame
    #[clap(short = 'b', long = "begin", default_value = "0")]
    pub begin: String,

    /// End time/frame  
    #[clap(short = 'e', long = "end", default_value = "")]
    pub end: String,

    /// Frame skip interval
    #[clap(long = "skip", default_value = "1")]
    pub skip: usize,
}

#[derive(Error, Debug)]
pub enum AnalysisError {
    #[error(transparent)]
    ParseFloat(#[from] ParseFloatError),

    #[error(transparent)]
    ParseInt(#[from] ParseIntError),

    #[error("invalid time suffix, 'fr', 'ps', 'ns', 'us' allowed")]
    InvalidSuffix,

    #[error(transparent)]
    FileIo(#[from] FileIoError),

    #[error(transparent)]
    Selection(#[from] SelectionError),

    #[error("in task pre_process: {0}")]
    PreProcess(#[source] anyhow::Error),

    #[error("in task process_frame: {0}")]
    ProcessFrame(#[source] anyhow::Error),

    #[error("in task post_process: {0}")]
    PostProcess(#[source] anyhow::Error),

    #[error("argument parsing")]
    Arg(#[from] clap::Error)
}

fn process_suffix(s: &str) -> Result<(Option<usize>, Option<f32>), AnalysisError> {
    if s.is_empty() {
        return Ok((None, None));
    }

    let mut frame = None;
    let mut time = None;

    if s.len() >= 2 {
        let suffix = &s[s.len() - 2..];
        let num_part = &s[..s.len() - 2];

        if let Ok(fr) = s.parse() {
            frame = Some(fr);
        } else {
            match suffix {
                "fr" => frame = Some(num_part.parse::<usize>()?),
                "ps" => time = Some(num_part.parse::<f32>()? * 1.0),
                "ns" => time = Some(num_part.parse::<f32>()? * 1000.0),
                "us" => time = Some(num_part.parse::<f32>()? * 1000_000.0),
                _ => return Err(AnalysisError::InvalidSuffix),
            }
        }
    }

    Ok((frame, time))
}

pub trait AnalysisTask {
    type Options: clap::Args;

    fn pre_process(&mut self, context: &AnalysisContext<Self::Options>) -> anyhow::Result<()>;

    fn process_frame(&mut self, context: &AnalysisContext<Self::Options>) -> anyhow::Result<()>;

    fn post_process(&mut self, context: &AnalysisContext<Self::Options>) -> anyhow::Result<()>;
    
    fn run(&mut self) -> Result<(), AnalysisError> {
        // Get the generic command line arguments
        let mut cmd = AnalysisArgs::command();
        // Add custom arguments from the implementor
        cmd = Self::Options::augment_args(cmd);
        
        let matches = cmd.get_matches();
        // Trajectory processing arguments
        let traj_args = AnalysisArgs::from_arg_matches(&matches)?;
        // Custom arguments
        let custom_args = Self::Options::from_arg_matches(&matches)?;
        
        // Greeting
        crate::greeting("molar_bin");
    
        if traj_args.files.len() < 2 {
            panic!("At least one trajectory file is required");
        }
    
        // Read topology
        let top: Holder<Topology, MutableSerial> =
            FileHandler::open(&traj_args.files[0])?.read_topology()?.into();
    
        let (begin_frame, begin_time) = process_suffix(&traj_args.begin)?;
        let (end_frame, end_time) = process_suffix(&traj_args.end)?;
    
        let n = top.num_atoms();
        let mut context = AnalysisContext {
            consumed_frames: 0,
            src: Source::new(top, State::new_fake(n))?,
            args: custom_args,
        };
    
        let mut begin_skipped = false;
    
        // Process trajectory files
        for trj_file in &traj_args.files[1..] {
            info!("Processing trajectory '{}'...", trj_file.display());
            let mut trj_handler = FileHandler::open(trj_file)?;
    
            if !begin_skipped {
                if let Some(fr) = begin_frame {
                    trj_handler.skip_to_frame(fr)?;
                } else if let Some(t) = begin_time {
                    trj_handler.skip_to_time(t)?;
                }
                begin_skipped = true;
            }
    
            let mut trj_iter = trj_handler.into_iter();
    
            while let Some(state) = trj_iter.next() {
                // Check if end is reached
                if end_frame.map_or(false, |ef| context.consumed_frames >= ef)
                    || end_time.map_or(false, |et| state.get_time() > et)
                {
                    break;
                }
    
                // Skip frames if needed
                if context.consumed_frames >0 && (context.consumed_frames - 1) % traj_args.skip > 0 {
                    continue;
                }
    
                if context.consumed_frames % traj_args.log == 0 {
                    info!(
                        "At frame {}, time {}",
                        context.consumed_frames,
                        get_log_time(&state)
                    );
                }
    
                context.src.set_state(state)?;
    
                if context.consumed_frames == 0 {
                    self.pre_process(&context).map_err(AnalysisError::PreProcess)?;
                }
    
                context.consumed_frames += 1;
    
                self.process_frame(&context).map_err(AnalysisError::ProcessFrame)?;
            }
        }
    
        self.post_process(&context).map_err(AnalysisError::PostProcess)?;
        Ok(())
    }
}

fn run_analysis_task<T: AnalysisTask>() -> anyhow::Result<()> {
    Ok(())
}

pub struct AnalysisContext<A> {
    pub src: Source<MutableSerial>,
    pub consumed_frames: usize,
    pub args: A,
}


fn get_log_time(state: &State) -> String {
    if state.get_time() < 1000.0 {
        format!("{} ps", state.get_time())
    } else if state.get_time() < 1000_000.0 {
        format!("{} ns", state.get_time() / 1000.0)
    } else {
        format!("{} us", state.get_time() / 1000_000.0)
    }
}
