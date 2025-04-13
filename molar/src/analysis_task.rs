use anyhow::bail;
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
    /// Input files (topology and at least one trajectory)
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

    #[error("no frames consumed")]
    NoFramesConsumed,

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

pub trait AnalysisTask<A: clap::Args> {
    //type Options: clap::Args;
    fn new(context: &AnalysisContext<A>) -> anyhow::Result<Self> where Self: Sized;

    //fn pre_process(&mut self, context: &AnalysisContext<A>) -> anyhow::Result<()>;

    fn process_frame(&mut self, context: &AnalysisContext<A>) -> anyhow::Result<()>;

    fn post_process(&mut self, context: &AnalysisContext<A>) -> anyhow::Result<()>;
    
    fn run() -> Result<(), AnalysisError> where Self: Sized {
        // Get the generic command line arguments
        let mut cmd = AnalysisArgs::command();
        // Add custom arguments from the implementor
        cmd = A::augment_args(cmd);
        
        let matches = cmd.get_matches();
        // Trajectory processing arguments
        let traj_args = AnalysisArgs::from_arg_matches(&matches)?;
        
        // Greeting
        crate::greeting("molar_bin");
    
        // Instance of implementing class (to be created on first valid frame)
        let mut inst: Option<Self> = None;

        // Analysis context (to be created on first valid frame)
        let mut context: Option<AnalysisContext<A>> = None;
        
        let (begin_frame, begin_time) = process_suffix(&traj_args.begin)?;
        let (end_frame, end_time) = process_suffix(&traj_args.end)?;
    
        // Process trajectory files
        let mut consumed_frames = 0;

        // Cycle over trajectory files (file 0 is structure)
        for trj_file in &traj_args.files[1..] {
            info!("Processing trajectory '{}'...", trj_file.display());
            let mut trj_handler = FileHandler::open(trj_file)?;
    
            // Skip frames in the beginning if asked
            if let Some(fr) = begin_frame {
                trj_handler.skip_to_frame(fr)?;
            } else if let Some(t) = begin_time {
                trj_handler.skip_to_time(t)?;
            }
    
            let mut trj_iter = trj_handler.into_iter();
            
            while let Some(state) = trj_iter.next() {
                // Check if end is reached
                if end_frame.map_or(false, |ef| consumed_frames >= ef)
                    || end_time.map_or(false, |et| state.get_time() > et)
                {
                    break;
                }
    
                // Skip frames if needed
                if consumed_frames >0 && (consumed_frames - 1) % traj_args.skip > 0 {
                    continue;
                }
    
                if consumed_frames % traj_args.log == 0 {
                    info!(
                        "At frame {}, time {}",
                        consumed_frames,
                        get_log_time(&state)
                    );
                }
    
                // If we are here than the frame is valid
                if consumed_frames == 0 {
                    // Time for initialization. It only happens ones on the first frame 
                    // of the first trajectory.
                    
                    // Read topology
                    let top = FileHandler::open(&traj_args.files[0])?.read_topology()?;
                    // Custom instance arguments
                    let args = A::from_arg_matches(&matches)?;
                    // Create context
                    context = Some(AnalysisContext {
                        args,
                        consumed_frames,
                        src: Source::new(top, state)?,
                    });
                    
                    // Create analysis object instance
                    inst = Some(Self::new(&context.as_ref().unwrap()).map_err(AnalysisError::PreProcess)?);
                } else if let Some(ref ctx) = context {
                    // Not the first frame. Update the context and call process_frame
                    // of the analysis instance
                    ctx.src.set_state(state)?;
                    inst.as_mut().unwrap().process_frame(ctx).map_err(AnalysisError::ProcessFrame)?;
                }
                // Update frame counter
                consumed_frames += 1;
            }
        } // trajectories
    
        if let Some(i) = inst.as_mut() {
            i.post_process(&context.as_ref().unwrap()).map_err(AnalysisError::PostProcess)?;
        } else {
            return Err(AnalysisError::NoFramesConsumed);
        }
        Ok(())
    }
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
