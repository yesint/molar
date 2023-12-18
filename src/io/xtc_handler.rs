use super::{IoReader, IoStateReader, IoStateWriter, IoWriter, IoRandomAccess};
use molar_xdrfile::xdrfile_bindings::*;
use nalgebra::{Matrix3, Point3};
use num_traits::ToPrimitive;

use crate::core::{State, PeriodicBox};

use anyhow::{bail, Result};
use std::ffi::CString;
use std::ptr;
use std::sync::{Arc, RwLock};

pub struct XtcFileHandler {
    handle: *mut XDRFILE,
    file_name: String,
    natoms: usize,
    steps_per_frame: usize,
    num_frames: usize,
    is_random_access: bool,
    dt: f32,
    max_t: f32,
    step: i32,
}

fn open_xdr_file(fname: &str, mode: &str) -> Result<*mut XDRFILE> {
    let c_name = CString::new(fname).unwrap();
    let c_mode = CString::new(mode).unwrap();
    let handle = unsafe {
        xdrfile_open(c_name.as_ptr(), c_mode.as_ptr())
    };

    if handle == ptr::null_mut() {
        bail!("Can't open file {} in mode {}!", fname, mode);
    }

    Ok(handle)
}

impl XtcFileHandler {
    fn new(fname: &str) -> Self {    
        XtcFileHandler {
            handle: ptr::null_mut(),
            file_name: fname.to_owned(),
            natoms: 0,
            steps_per_frame: 0,
            num_frames: 0,
            is_random_access: true,
            dt: 0.0,
            max_t: 0.0,
            step: 0,
        }
    }

    fn open_read(&mut self) -> Result<()> {
        self.handle = open_xdr_file(&self.file_name,"r")?;

        // Extract number of atoms
        let mut n_at: i32 = 0;
        let ok = unsafe { xdr_xtc_get_natoms(self.handle, &mut n_at) };
        if ok != 1 {
            // 1 is success for this function, f**ing idiots!
            bail!("Can't read XTC number of atoms {} {}", ok, n_at);
        }

        self.natoms = n_at as usize;

        // XTC file contains step number in terms of simulation steps, not saved frames
        // So we have to extract conversion factor
        let next = unsafe { xtc_get_next_frame_number(self.handle, n_at) };
        let mut ok: bool = false;
        let cur = unsafe { xtc_get_current_frame_number(self.handle, n_at, &mut ok) };
        if cur < 0 || next < 0 || !ok {
            bail!("Can't detect number of steps per frame");
        }
        if cur == next {
            println!("It seems that there is only one frame in this trajectory");
            self.steps_per_frame = 1;
        } else {
            self.steps_per_frame = (next - cur) as usize;
        }

        // Get total number of frames in the trajectory
        let n_frames = unsafe {
            xdr_xtc_get_last_frame_number(self.handle, n_at, &mut ok)
        };

        if !ok {
            bail!("Can't get number of frames");
        }

        if n_frames < 0 {
            println!(
                "Weird XTC file: negative number of frames returned ({})!",
                n_frames
            );
            // Disable random access
            self.is_random_access = false;
        }

        self.num_frames = n_frames as usize / self.steps_per_frame;

        // Get time step
        self.dt = unsafe {
            xdr_xtc_estimate_dt(self.handle, n_at, &mut ok)
        };

        if !ok {
            println!("Can't get time step");
            self.dt = -1.0;
        }

        self.max_t = unsafe { xdr_xtc_get_last_frame_time(self.handle, n_at, &mut ok) };
        if !ok || self.max_t < 0.0 {
            println!("Can't get last frame time");
            self.max_t = -1.0;
        }

        if !self.is_random_access {
            println!("Random access operations disabled for this trajectory");
        }

        println!(
            "File: {}, natoms={}, nframes={}, max_t={}, dt={}",
            self.file_name, self.natoms, self.num_frames, self.max_t, self.dt
        );

        Ok(())
    }


    fn open_write(&mut self) -> Result<()> {
        self.handle = open_xdr_file(&self.file_name,"w")?;
        Ok(())
    }
}


impl IoReader for XtcFileHandler {
    fn new_reader(fname: &str) -> Result<Self> {
        let mut instance = Self::new(fname);
        instance.open_read()?;
        Ok(instance)
    }
}

impl IoWriter for XtcFileHandler {
    fn new_writer(fname: &str) -> Result<Self> {
        let mut instance = Self::new(fname);
        instance.open_write()?;
        Ok(instance)
    }
}

impl Drop for XtcFileHandler {
    fn drop(&mut self) {
        if self.handle != ptr::null_mut() {
            let ok = unsafe { xdrfile_close(self.handle) } as u32;
            if ok != exdrOK {
                panic!("Error closing XTC file {}!",self.file_name)
            }
        }
    }
}

impl IoStateReader for XtcFileHandler {
    #[allow(non_upper_case_globals)]
    fn read_next_state(&mut self) -> Result<Option<State>> {
        let mut st: State = Default::default();
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
                box_matrix.as_mut_ptr().cast::<[f32;3]>(),
                st.coords.as_mut_ptr().cast::<rvec>(),
                &mut prec,
            )
        };

        if ok as u32 == exdrOK {
            // C function populated the coordinates, set the vector size for Rust
            unsafe { st.coords.set_len(self.natoms) };
            // Convert box to column-major form.
            box_matrix.transpose_mut();
            st.box_ = PeriodicBox::from_matrix(box_matrix).ok();
        }
        
        match ok as u32 {
            exdrOK => Ok(Some(st.into())),
            exdrENDOFFILE => Ok(None),
            _ => bail!("Error reading timestep!"),
        }
    }
}

impl IoStateWriter for XtcFileHandler {
    fn write_next_state_subset(&mut self, data: &State, 
            subset_indexes: impl ExactSizeIterator<Item=usize>) -> Result<()> 
    {
        let n = subset_indexes.len();

        // Box have to be transposed because XTC contains row-major box
        let box_ = match data.box_.as_ref() {
            Some(b) => b.get_matrix().transpose(),
            None => Matrix3::<f32>::zeros(),
        };

        // Coordinate buffer
        let mut buf = Vec::<Point3<f32>>::new();
        // Pointer to coordinates
        let mut coord_ptr: *const Point3<f32> = data.coords.as_ptr();

        // If not all coordinates are written we have to extract them
        // to a buffer instead of passing the pointer to original coords
        if n != data.coords.len() {
            // Fill the buffer
            buf.reserve(n);
            for ind in subset_indexes {
                buf.push(data.coords[ind]);
            }
            // Reset the pointer to buffer
            coord_ptr = buf.as_ptr();
        }

        let ok = unsafe {
            write_xtc(
                self.handle,
                n as i32,
                self.step,
                data.time,
                box_.as_ptr().cast::<[f32;3]>(),
                coord_ptr.cast::<rvec>(),
                1000.0,
            )
        };

        if ok as u32 != exdrOK {
            bail!("Unable to write time step {} to XTC file! Return code: {}",self.step,ok);
        }

        self.step+=1;

        Ok(())
    }
}

impl IoRandomAccess for XtcFileHandler {
    fn seek_frame(&mut self, fr: usize) -> Result<()> {
        if !self.is_random_access {
            bail!("Random access is not possible for this XTC file!");
        }

        if fr>self.num_frames {
            bail!("Can't seek to frame {}, last frame is {}",fr,self.num_frames);
        }

        let ret = unsafe{
            xdr_xtc_seek_frame(
                (fr*self.steps_per_frame).try_into()?,
                self.handle,
                self.natoms.try_into()?
            )
        };
        if ret<0 {
            bail!("Error seeking to frame {}",fr);
        }

        Ok(())
    }

    fn seek_time(&mut self, t: f32) -> Result<()> {
        if !self.is_random_access {
            bail!("Random access is not possible for this XTC file!");
        }

        if t<0.0 || t>self.max_t {
            bail!("Can't seek to time {}, last time is {}",t,self.max_t);
        }
        // We assume equally spaced frames in the trajectory. It's much faster
        let ret = unsafe {
            xdr_xtc_seek_frame(
                ((t/self.dt).ceil()*self.steps_per_frame as f32).to_int_unchecked(),
                self.handle,
                self.natoms.try_into()?
            )
        };
        //int ret = xdr_xtc_seek_time(t,handle,natoms,false);
        if ret<0 {
            bail!("Error seeking to time {}",t)
        }

        Ok(())
    }

    fn tell_current(&self) -> Result<(usize,f32)> {
        let mut ok = false;
        let ret = unsafe{
            xtc_get_current_frame_number(
                self.handle,
                self.natoms.try_into()?,
                &mut ok)
        };
        if !ok || ret<0 {
            bail!("Can't get current frame number");
        }
        let step = ret/self.steps_per_frame as i32;
        let t = unsafe{
            xtc_get_current_frame_time(
                self.handle,
                self.natoms.try_into()?,
                &mut ok)
        };
        if !ok || t<0.0 {
            bail!("Can't get current frame time");
        }

        Ok((step as usize,t))
    }

    fn tell_last(&self) -> Result<(usize,f32)> {
        if self.max_t<0.0 {
            bail!("This XTC file doesn't allow getting last time!");
        }

        Ok((self.num_frames,self.max_t))
    }
}