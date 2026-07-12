//! Support for custom analysis tasks

use clap::{CommandFactory, FromArgMatches, Parser};
use log::info;
use std::num::{ParseFloatError, ParseIntError};
use std::path::PathBuf;
use thiserror::Error;

use crate::prelude::*;

/// Standard trajectory processing arguments
#[derive(Parser)]
#[command(name = "analysis")]
pub struct TrajAnalysisArgs {
    /// Input files (topology and at least one trajectory)
    #[clap(short = 'f', long = "files", required = true, num_args = 1..)]
    pub files: Vec<PathBuf>,

    /// Logging frequency
    #[clap(long = "log", default_value_t = 100)]
    pub log: usize,

    /// Begin time/frame
    #[clap(short = 'b', long = "begin", default_value = "0")]
    pub begin: String,

    /// End time/frame  
    #[clap(short = 'e', long = "end", default_value = "")]
    pub end: String,

    /// Frame skip interval (process every N-th frame; must be >= 1)
    #[clap(long = "skip", default_value_t = 1, value_parser = clap::builder::RangedU64ValueParser::<usize>::new().range(1..))]
    pub skip: usize,

    /// Use coordinates in structure file for analysis
    #[clap(long = "use_struct_file", default_value_t = false)]
    pub use_struct_file: bool,
}

/// Errors related to analysis task execution and trajectory processing
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
    Arg(#[from] clap::Error),

    #[error("at least one trajectory required if 'use_struct_file' is not set")]
    NoTraj,
}

/// Parses a `-b`/`-e` value into an optional frame index and/or time (in ps).
///
/// An empty string yields `(None, None)` (no limit). A bare number is a frame
/// index; a number with an explicit unit suffix is either a frame (`fr`) or a
/// time (`ps`/`ns`/`us`, all converted to ps).
fn process_suffix(s: &str) -> Result<(Option<usize>, Option<Float>), AnalysisError> {
    let s = s.trim();
    if s.is_empty() {
        return Ok((None, None));
    }

    // A bare number (no unit) is interpreted as a frame index.
    if let Ok(fr) = s.parse::<usize>() {
        return Ok((Some(fr), None));
    }

    // Explicit frame suffix.
    if let Some(num) = s.strip_suffix("fr") {
        return Ok((Some(num.trim().parse::<usize>()?), None));
    }

    // Time suffixes, all converted to ps.
    for (suffix, scale) in [
        ("ps", 1.0 as Float),
        ("ns", 1000.0 as Float),
        ("us", 1_000_000.0 as Float),
    ] {
        if let Some(num) = s.strip_suffix(suffix) {
            return Ok((None, Some(num.trim().parse::<Float>()? * scale)));
        }
    }

    Err(AnalysisError::InvalidSuffix)
}

/// Analysis task trait. Should be implemented by user's types.
pub trait AnalysisTask<A: clap::Args> {
    fn new(context: &mut AnalysisContext<A>) -> anyhow::Result<Self>
    where
        Self: Sized;

    fn process_frame(&mut self, context: &mut AnalysisContext<A>) -> anyhow::Result<()>;

    fn post_process(&mut self, context: &mut AnalysisContext<A>) -> anyhow::Result<()>;

    fn task_name() -> String;

    fn run() -> Result<(), AnalysisError>
    where
        Self: Sized,
    {
        // Get the generic command line arguments
        let mut cmd = TrajAnalysisArgs::command();
        // Add custom arguments from the implementor
        cmd = A::augment_args(cmd);

        let matches = cmd.get_matches();
        // Trajectory processing arguments
        let traj_args = TrajAnalysisArgs::from_arg_matches(&matches)?;

        // Greeting
        crate::greeting(Self::task_name());

        if !traj_args.use_struct_file && traj_args.files.len() < 2 {
            return Err(AnalysisError::NoTraj);
        }

        // Instance of implementing class (to be created on first valid frame)
        let mut inst: Option<Self> = None;
        // Analysis context (to be created on first valid frame)
        let mut context: Option<AnalysisContext<A>> = None;

        let (begin_frame, begin_time) = process_suffix(&traj_args.begin)?;
        let (end_frame, end_time) = process_suffix(&traj_args.end)?;

        // Number of frames actually processed (passed to the task).
        let mut consumed_frames = 0;
        // Number of frames read since `begin`. This advances for *every* frame,
        // including skipped ones. `--skip` is relative to this counter, so the
        // `begin` frame is always the first one processed.
        let mut read_frames = 0;
        // Absolute index of the first read frame, used to make the `-e <frame>`
        // window bound absolute (consistent with `-e <time>` and `-b`). When
        // `begin` is given as a time the absolute index is unknown, so it falls
        // back to 0 (i.e. relative to the begin-time position).
        let begin_offset = begin_frame.unwrap_or(0);

        // Process state from structure file if asked
        if traj_args.use_struct_file {
            info!("Using structure file for task initialization");
            // Read both topology and state
            let (top, state) = FileHandler::open(&traj_args.files[0])?.read()?;
            // Custom instance arguments
            let args = A::from_arg_matches(&matches)?;
            // Do init
            (inst, context) = init_context_and_process_first_frame(args, top, state)?;
            // One frame is consumed
            consumed_frames += 1;
        }

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
                // Absolute trajectory frame index of the current frame.
                let abs_frame = begin_offset + read_frames;

                // Check if the end of the requested window is reached. `-e` is
                // an absolute bound (frame index or time), exclusive.
                if end_frame.map_or(false, |ef| abs_frame >= ef)
                    || end_time.map_or(false, |et| state.get_time() > et)
                {
                    break;
                }

                // Process only every `skip`-th frame. `read_frames` advances for
                // every frame read, so the cadence follows the trajectory
                // position rather than how many frames were kept (which was the
                // bug that capped processing at 2 frames for any `--skip > 1`).
                let frame_index = read_frames;
                read_frames += 1;
                if frame_index % traj_args.skip != 0 {
                    continue;
                }

                if consumed_frames % traj_args.log == 0 {
                    info!("At frame {}, time {}", abs_frame, get_log_time(&state));
                }

                // If we are here than the frame is valid
                if let Some(ref mut ctx) = context {
                    // Context is initialized already, update it with current state
                    ctx.sys.set_state(state)?;
                    // Process frame
                    inst.as_mut()
                        .unwrap()
                        .process_frame(ctx)
                        .map_err(AnalysisError::ProcessFrame)?;
                } else {
                    // Time for context initialization
                    info!("Using first trajectory frame for task initialization");
                    // Read topology, state is read from trajectory already
                    let top = FileHandler::open(&traj_args.files[0])?.read_topology()?;
                    // Custom instance arguments
                    let args = A::from_arg_matches(&matches)?;
                    // Do init
                    (inst, context) = init_context_and_process_first_frame(args, top, state)?;
                }

                // Update frame counter
                consumed_frames += 1;
                context.as_mut().unwrap().consumed_frames += 1;
            }
            info!("Finished with '{}'.", trj_file.display());
        } // trajectories

        if let Some(inst) = inst.as_mut() {
            info!("Post-processing...");
            inst.post_process(&mut context.as_mut().unwrap())
                .map_err(AnalysisError::PostProcess)?;
        } else {
            return Err(AnalysisError::NoFramesConsumed);
        }
        Ok(())
    }
}

fn init_context_and_process_first_frame<A, T>(
    args: A,
    top: Topology,
    state: State,
) -> Result<(Option<T>, Option<AnalysisContext<A>>), AnalysisError>
where
    A: clap::Args,
    T: AnalysisTask<A>,
{
    info!("Initializing analysis task instance");
    // Create context
    let mut context = AnalysisContext {
        args,
        consumed_frames: 0,
        sys: System::new(top, state)?,
    };

    // Create analysis object instance
    let mut inst = T::new(&mut context).map_err(AnalysisError::PreProcess)?;
    // Call process frame
    inst.process_frame(&mut context)
        .map_err(AnalysisError::ProcessFrame)?;

    Ok((Some(inst), Some(context)))
}

/// Context passed to all frame processing methods
pub struct AnalysisContext<A> {
    pub sys: System,
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn suffix_empty_is_no_limit() {
        assert_eq!(process_suffix("").unwrap(), (None, None));
        assert_eq!(process_suffix("   ").unwrap(), (None, None));
    }

    #[test]
    fn suffix_bare_number_is_frame() {
        // Regression: single-digit values used to be silently dropped.
        assert_eq!(process_suffix("0").unwrap(), (Some(0), None));
        assert_eq!(process_suffix("5").unwrap(), (Some(5), None));
        assert_eq!(process_suffix("42").unwrap(), (Some(42), None));
        assert_eq!(process_suffix("100").unwrap(), (Some(100), None));
    }

    #[test]
    fn suffix_explicit_frame() {
        assert_eq!(process_suffix("5fr").unwrap(), (Some(5), None));
        assert_eq!(process_suffix("100fr").unwrap(), (Some(100), None));
    }

    #[test]
    fn suffix_time_units_convert_to_ps() {
        assert_eq!(process_suffix("5ps").unwrap(), (None, Some(5.0)));
        assert_eq!(process_suffix("2ns").unwrap(), (None, Some(2000.0)));
        assert_eq!(process_suffix("1us").unwrap(), (None, Some(1_000_000.0)));
        assert_eq!(process_suffix("1.5ns").unwrap(), (None, Some(1500.0)));
    }

    #[test]
    fn suffix_invalid() {
        assert!(matches!(
            process_suffix("5km"),
            Err(AnalysisError::InvalidSuffix)
        ));
        // A unit with no number is a parse error, not a panic.
        assert!(process_suffix("fr").is_err());
    }
}
