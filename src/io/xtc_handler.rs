use super::xdrfile_bindings::{matrix, xdrfile_close, xdrfile_open, XDRFILE,
    xdr_xtc_get_natoms,xtc_get_next_frame_number,xtc_get_current_frame_number,
    xdr_xtc_estimate_dt,xdr_xtc_get_last_frame_number, xdr_xtc_get_last_frame_time, exdrOK};

use std::ptr;
use std::ffi::CString;
use anyhow::{Result, bail};

pub struct XtcFileHandler {
    handle: *mut XDRFILE,
    pbox: matrix,
    file_name: String,
    steps_per_frame: u32,
    num_frames: u32,
    is_random_access: bool,
    dt: f32,
    max_t: f32,
}

impl XtcFileHandler {
    pub fn new(fname: &str) -> Self {
        XtcFileHandler {
            handle: ptr::null_mut(),
            pbox: Default::default(),
            file_name: fname.to_owned(),
            steps_per_frame: 0,
            num_frames: 0,
            is_random_access: true,
            dt: 0.0,
            max_t: 0.0,
        }
    }

    fn open_read(&mut self) -> Result<()> {
        let f_name = CString::new(self.file_name.clone()).unwrap();
        let mode = CString::new("r").unwrap();
        self.handle = unsafe{ xdrfile_open(f_name.as_ptr(), mode.as_ptr()) };

        if self.handle == ptr::null_mut() {
            bail!("Can't open file {} for reading!", self.file_name);
        }

        // Extract number of atoms
        let mut natoms: i32 = 0;
        let ok = unsafe{ xdr_xtc_get_natoms(self.handle, &mut natoms) } as u32;
        if ok!=exdrOK {
            bail!("Can't read XTC number of atoms");
        }

        // XTC file contains step number in terms of simulation steps, not saved frames
        // So we have to extract conversion factor
        let next = unsafe{ xtc_get_next_frame_number(self.handle,natoms) };
        let mut b_ok: bool = false;
        let cur = unsafe{ xtc_get_current_frame_number(self.handle,natoms,&mut b_ok) };
        if cur<0 || next<0 || !b_ok {
            bail!("Can't detect number of steps per frame");
        }
        if cur==next {
            println!("It seems that there is only one frame in this trajectory");
            self.steps_per_frame = 1;
        } else {
            self.steps_per_frame = (next-cur) as u32;
        }

        // Get total number of frames in the trajectory
        let n_frames = unsafe {
            xdr_xtc_get_last_frame_number(self.handle,natoms,&mut b_ok)
        };

        if !b_ok {
            bail!("Can't get number of frames");
        }

        if n_frames<0 {
            println!("Weird XTC file: negative number of frames returned ({})!",n_frames);
            // Disable random access
            self.is_random_access = false;
        }

        self.num_frames = n_frames as u32 / self.steps_per_frame;

        // Get time step
        self.dt = unsafe { xdr_xtc_estimate_dt(self.handle,natoms,&mut b_ok) };

        if !b_ok {
            println!("Can't get time step");
            self.dt = -1.0;
        }

        self.max_t = unsafe{ xdr_xtc_get_last_frame_time(self.handle,natoms,&mut b_ok) };
        if !b_ok || self.max_t<0.0 {
            println!("Can't get last frame time");
            self.max_t = -1.0;
        }

        if !self.is_random_access {
            println!("Random access operations disabled for this trajectory");
        }

        //LOG()->debug("There are {} frames, max_t={}, dt={}",
        //             xtc->num_frames,
        //             xtc->max_t ? fmt::format("{}",xtc->max_t) : "N/A",
        //             xtc->dt ? fmt::format("{}",xtc->dt) : "N/A");

        Ok(())
    }


}

impl Drop for XtcFileHandler {
    fn drop(&mut self) {
        if self.handle != ptr::null_mut() {
            let ok = unsafe{ xdrfile_close(self.handle) } as u32;
            if ok != exdrOK {
                panic!("Error closing XTC file!")
            }
        }
    }
}