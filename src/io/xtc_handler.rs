use super::{PeriodicBoxError, State, StateProvider};
use molar_xdrfile::xdrfile_bindings::*;
use nalgebra::{Matrix3, Point3};
use thiserror::Error;

use crate::core::{PeriodicBox, StateStorage};

use log::{debug, warn};
use std::ffi::{CString, NulError};
use std::ptr;

pub struct XtcFileHandler {
    handle: *mut XDRFILE,
    natoms: usize,
    steps_per_frame: usize,
    dt: f32,
    step: i32, // time step for C functions
    frame_range: (usize, usize),
    time_range: (f32, f32),
    is_random_access: bool,
}

// Allow sending handler between threads
unsafe impl Send for XtcFileHandler {}

#[derive(Error, Debug)]
pub enum XtcHandlerError {
    #[error("unexpected null characted")]
    CStringNull(#[from] NulError),

    #[error("plugin can't open file in mode {0}")]
    Open(String),

    #[error("failed to read state")]
    ReadState,

    #[error("failed to write state")]
    WriteState,

    #[error("invalid periodic box")]
    Pbc(#[from] PeriodicBoxError),

    #[error("fixed size field is {0} while needed {1}")]
    FixedSizeFieldOverflow(usize, usize),

    #[error("can't read number of atoms")]
    ReadNumAtoms,

    #[error("can't get current frame number")]
    CantGetCurrentFrameNumber,

    #[error("can't get number of steps per frame")]
    CantGetNumberOfStepsPerFrame,

    #[error("can't get current frame time")]
    CantGetCurrentFrameTime,

    #[error("can't get last frame number")]
    CantGetLastFrameNumber,

    #[error("can't get last frame time")]
    CantGetLastFrameTime,

    #[error("this xtc file doesn't allow random access")]
    RandomAccessImpossible,

    #[error("can't seek to frame {0}, allowed range is {1:?}")]
    SeekOutOfBounds(usize, (usize, usize)),

    #[error(transparent)]
    IntegerConvertion(#[from] std::num::TryFromIntError),

    #[error("unknown error seeking to frame {0}")]
    SeekFailed(usize),

    #[error("can't seek to time {0}, allowed range is {1:?}")]
    SeekTimeOutOfBounds(f32, (f32, f32)),

    #[error("unknown error seeking to time {0}")]
    SeekTimeFailed(f32),
}

fn get_xdr_handle(fname: &str, mode: &str) -> Result<*mut XDRFILE, XtcHandlerError> {
    let c_name = CString::new(fname)?;
    let c_mode = CString::new(mode)?;
    let handle = unsafe { xdrfile_open(c_name.as_ptr(), c_mode.as_ptr()) };

    if handle == ptr::null_mut() {
        Err(XtcHandlerError::Open(mode.into()))
    } else {
        Ok(handle)
    }
}

impl Default for XtcFileHandler {
    fn default() -> Self {
        Self {
            handle: ptr::null_mut(),
            natoms: 0,
            steps_per_frame: 0,
            frame_range: (0, 0),
            is_random_access: true,
            dt: 0.0,
            time_range: (0.0, 0.0),
            step: 0,
        }
    }
}

impl XtcFileHandler {
    fn open_read(&mut self, fname: &str) -> Result<(), XtcHandlerError> {
        self.handle = get_xdr_handle(fname, "r")?;

        // Extract number of atoms
        let mut n_at: i32 = 0;
        let ok = unsafe { xdr_xtc_get_natoms(self.handle, &mut n_at) };
        if ok != 1 {
            // 1 is success for this function, f**ing idiots!
            return Err(XtcHandlerError::ReadNumAtoms);
        }

        self.natoms = n_at as usize;

        // XTC file contains step number in terms of simulation steps, not saved frames
        // So we have to extract conversion factor
        let mut ok: bool = false;
        let first_step = unsafe { xtc_get_current_frame_number(self.handle, n_at, &mut ok) };
        if !ok {
            return Err(XtcHandlerError::CantGetCurrentFrameNumber);
        }
        let next_step = unsafe { xtc_get_next_frame_number(self.handle, n_at) };
        if next_step < first_step {
            return Err(XtcHandlerError::CantGetNumberOfStepsPerFrame);
        } else if first_step == next_step {
            self.steps_per_frame = 1;
        } else {
            self.steps_per_frame = (next_step - first_step) as usize;
        }

        // Get first time
        let first_time = unsafe { xtc_get_current_frame_time(self.handle, n_at, &mut ok) };
        if !ok {
            return Err(XtcHandlerError::CantGetCurrentFrameTime);
        }

        // Get last frame
        let last_step = unsafe { xdr_xtc_get_last_frame_number(self.handle, n_at, &mut ok) };
        if !ok {
            return Err(XtcHandlerError::CantGetLastFrameNumber);
        }
        // Get last time
        let last_time = unsafe { xdr_xtc_get_last_frame_time(self.handle, n_at, &mut ok) };
        if !ok {
            return Err(XtcHandlerError::CantGetLastFrameTime);
        }

        if last_step < first_step || last_time < first_time {
            warn!("Last frame seems to be corrupted ({last_step})!");
            // Disable random access
            self.is_random_access = false;
        }

        self.frame_range = (
            first_step as usize / self.steps_per_frame,
            last_step as usize / self.steps_per_frame,
        );

        self.time_range = (first_time, last_time);

        // Get dt
        self.dt = unsafe { xdr_xtc_estimate_dt(self.handle, n_at, &mut ok) };
        if !ok {
            warn!("Can't get dt");
            self.dt = -1.0;
        }

        if !self.is_random_access {
            warn!(target: "XtcHandler", "Random access operations are disabled for this trajectory");
        }

        debug!(
            "Xtc content:\nnatoms={}\nframes={:?}\ntimes={:?}\ndt={}\nsteps_per_frame={}\nrandom_access={}",
            self.natoms,
            self.frame_range,
            self.time_range,
            self.dt,
            self.steps_per_frame,
            self.is_random_access
        );

        Ok(())
    }

    fn open_write(&mut self, fname: &str) -> Result<(), XtcHandlerError> {
        self.handle = get_xdr_handle(fname, "w")?;
        Ok(())
    }

    pub fn open(fname: &str) -> Result<Self, XtcHandlerError> {
        let mut instance = Self::default();
        instance.open_read(fname)?;
        Ok(instance)
    }

    pub fn create(fname: &str) -> Result<Self, XtcHandlerError> {
        let mut instance = Self::default();
        instance.open_write(fname)?;
        Ok(instance)
    }

    #[allow(non_upper_case_globals)]
    pub fn read_state(&mut self) -> Result<Option<State>, XtcHandlerError> {
        let mut st: StateStorage = Default::default();
        // Prepare variables
        let mut prec: f32 = 0.0;
        // Allocate storage for coordinates, but don't initialize them
        // This doesn't waste time for initialization, which will be overwritten anyway
        st.coords = Vec::with_capacity(self.natoms);
        let mut box_matrix = Matrix3::<f32>::zeros();

        let ok = unsafe {
            read_xtc(
                self.handle,
                self.natoms as i32,
                &mut self.step,
                &mut st.time,
                box_matrix.as_mut_ptr().cast::<rvec>(),
                st.coords.as_mut_ptr().cast::<rvec>(),
                &mut prec,
            )
        };

        if ok as u32 == exdrOK {
            // C function populated the coordinates, set the vector size for Rust
            unsafe { st.coords.set_len(self.natoms) };
            // Convert box to column-major form.
            box_matrix.transpose_mut();
            st.pbox = Some(PeriodicBox::from_matrix(box_matrix)?);
        }

        match ok as u32 {
            exdrOK => Ok(Some(st.into())),
            exdrENDOFFILE => Ok(None),
            _ => Err(XtcHandlerError::ReadState),
        }
    }

    pub fn write_state(&mut self, data: &impl StateProvider) -> Result<(), XtcHandlerError> {
        let n = data.num_coords();

        // Box have to be transposed because XTC contains row-major box
        let box_ = match data.get_box() {
            Some(b) => b.get_matrix().transpose(),
            None => Matrix3::<f32>::zeros(),
        };

        // Coordinate buffer
        #[allow(unused_mut)]
        let mut buf = Vec::<Point3<f32>>::from_iter(data.iter_pos().cloned());

        let ok = unsafe {
            write_xtc(
                self.handle,
                n as i32,
                self.step,
                data.get_time(),
                box_.as_ptr().cast::<rvec>(),
                buf.as_ptr().cast::<rvec>(),
                1000.0,
            )
        };

        if ok as u32 != exdrOK {
            return Err(XtcHandlerError::WriteState);
        }

        self.step += 1;

        Ok(())
    }

    pub fn seek_frame(&mut self, fr: usize) -> Result<(), XtcHandlerError> {
        if !self.is_random_access {
            return Err(XtcHandlerError::RandomAccessImpossible);
        }

        if fr < self.frame_range.0 || fr > self.frame_range.1 {
            return Err(XtcHandlerError::SeekOutOfBounds(fr, self.frame_range));
        }

        let ret = unsafe {
            xdr_xtc_seek_frame(
                (fr * self.steps_per_frame).try_into()?,
                self.handle,
                self.natoms.try_into()?,
            )
        };
        if ret < 0 {
            return Err(XtcHandlerError::SeekFailed(fr));
        }

        Ok(())
    }

    pub fn seek_time(&mut self, t: f32) -> Result<(), XtcHandlerError> {
        if !self.is_random_access {
            return Err(XtcHandlerError::RandomAccessImpossible);
        }

        if t < self.time_range.0 || t > self.time_range.1 {
            return Err(XtcHandlerError::SeekTimeOutOfBounds(t, self.time_range));
        }
        // We assume equally spaced frames in the trajectory. It's much faster
        /*
        let ret = unsafe {
            xdr_xtc_seek_frame(
                ((t/self.dt).ceil()*self.steps_per_frame as f32).to_int_unchecked(),
                self.handle,
                self.natoms.try_into()?
            )
        };
        */
        let ret = unsafe { xdr_xtc_seek_time(t, self.handle, self.natoms.try_into()?, false) };
        if ret < 0 {
            return Err(XtcHandlerError::SeekTimeFailed(t));
        }

        Ok(())
    }

    pub fn tell_current(&self) -> Result<(usize, f32), XtcHandlerError> {
        let mut ok = false;
        let ret =
            unsafe { xtc_get_current_frame_number(self.handle, self.natoms.try_into()?, &mut ok) };
        if !ok || ret < 0 {
            return Err(XtcHandlerError::CantGetCurrentFrameNumber);
        }
        let step = ret / self.steps_per_frame as i32;
        let t =
            unsafe { xtc_get_current_frame_time(self.handle, self.natoms.try_into()?, &mut ok) };
        if !ok || t < 0.0 {
            return Err(XtcHandlerError::CantGetCurrentFrameTime);
        }

        Ok((step as usize, t))
    }

    pub fn tell_first(&self) -> Result<(usize, f32), XtcHandlerError> {
        if !self.is_random_access {
            return Err(XtcHandlerError::RandomAccessImpossible);
        }
        Ok((self.frame_range.0, self.time_range.0))
    }

    pub fn tell_last(&self) -> Result<(usize, f32), XtcHandlerError> {
        if !self.is_random_access {
            return Err(XtcHandlerError::RandomAccessImpossible);
        }
        Ok((self.frame_range.1, self.time_range.1))
    }

    // pub fn get_file_name(&self) -> &str {
    //     &self.file_name
    // }
}

impl Drop for XtcFileHandler {
    fn drop(&mut self) {
        if self.handle != ptr::null_mut() {
            let ok = unsafe { xdrfile_close(self.handle) } as u32;
            if ok != exdrOK {
                panic!("Error closing XTC file!")
            }
        }
    }
}
