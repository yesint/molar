use super::{IoReader, IoStateReader, IoStateWriter, IoWriter, IoRandomAccess, IoIndexAndStateProvider};
use molar_xdrfile::xdrfile_bindings::*;
use nalgebra::{Matrix3, Point3};

use crate::core::{State, PeriodicBox};

use anyhow::{bail, Result};
use std::ffi::CString;
use std::ptr;

pub struct XtcFileHandler {
    handle: *mut XDRFILE,
    file_name: String,
    natoms: usize,
    steps_per_frame: usize,
    dt: f32,
    step: i32, // time step for C functions
    frame_range: (usize,usize),
    time_range: (f32,f32),
    is_random_access: bool,
}

fn get_xdr_handle(fname: &str, mode: &str) -> Result<*mut XDRFILE> {
    let c_name = CString::new(fname)?;
    let c_mode = CString::new(mode)?;
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
            frame_range: (0,0),
            is_random_access: true,
            dt: 0.0,
            time_range: (0.0,0.0),
            step: 0,
        }
    }

    fn open_read(&mut self) -> Result<()> {
        self.handle = get_xdr_handle(&self.file_name,"r")?;

        // Extract number of atoms
        let mut n_at: i32 = 0;
        let ok = unsafe { xdr_xtc_get_natoms(self.handle, &mut n_at) };
        if ok != 1 {
            // 1 is success for this function, f**ing idiots!
            bail!("Can't read XTC number of atoms. Err: {ok}");
        }

        self.natoms = n_at as usize;

        // XTC file contains step number in terms of simulation steps, not saved frames
        // So we have to extract conversion factor
        let mut ok: bool = false;
        let first_step = unsafe { xtc_get_current_frame_number(self.handle, n_at, &mut ok) };
        if !ok { bail!("Can't get first simulation step"); }
        let next_step = unsafe { xtc_get_next_frame_number(self.handle, n_at) };
        if next_step < first_step {
            bail!("Can't detect number of steps per frame");
        } else if first_step == next_step {
            // It seems that there is only one frame in this trajectory
            self.steps_per_frame = 1;
        } else {
            self.steps_per_frame = (next_step - first_step) as usize;
        }

        // Get first time
        let first_time = unsafe { xtc_get_current_frame_time(self.handle, n_at, &mut ok) };
        if !ok { bail!("Can't get first frame time"); }

        // Get last frame
        let last_step = unsafe { xdr_xtc_get_last_frame_number(self.handle, n_at, &mut ok) };
        if !ok { bail!("Can't get last time step"); }
        // Get last time
        let last_time = unsafe { xdr_xtc_get_last_frame_time(self.handle, n_at, &mut ok) };
        if !ok { bail!("Can't get first frame time"); }
        
        if last_step < first_step || last_time < first_time {
            println!("Last frame seems to be corrupted ({last_step})!");
            // Disable random access
            self.is_random_access = false;
        }

        self.frame_range = (
            first_step as usize/ self.steps_per_frame,
            last_step as usize/ self.steps_per_frame
        );

        self.time_range = (first_time,last_time);

        // Get dt
        self.dt = unsafe {
            xdr_xtc_estimate_dt(self.handle, n_at, &mut ok)
        };
        if !ok {
            println!("Can't get dt");
            self.dt = -1.0;
        }

        if !self.is_random_access {
            println!("Random access operations are disabled for this trajectory");
        }

        println!(
            "File: {}, natoms={}, frames={:?}, times={:?}, dt={}",
            self.file_name, self.natoms, self.frame_range, self.time_range, self.dt
        );

        Ok(())
    }


    fn open_write(&mut self) -> Result<()> {
        self.handle = get_xdr_handle(&self.file_name,"w")?;
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
            st.box_ = Some(PeriodicBox::from_matrix(box_matrix)?);
        }
        
        match ok as u32 {
            exdrOK => Ok(Some(st)),
            exdrENDOFFILE => Ok(None),
            _ => bail!("Error reading timestep!"),
        }
    }
}

impl IoStateWriter for XtcFileHandler {
    fn write_next_state(&mut self, data: &impl IoIndexAndStateProvider) -> Result<()> 
    {
        let (index,st) = data.get_index_and_state();
        let n = index.len();

        // Box have to be transposed because XTC contains row-major box
        let box_ = match st.box_.as_ref() {
            Some(b) => b.get_matrix().transpose(),
            None => Matrix3::<f32>::zeros(),
        };

        // Coordinate buffer
        let mut buf = Vec::<Point3<f32>>::new();
        // Pointer to coordinates
        let mut coord_ptr: *const Point3<f32> = st.coords.as_ptr();

        // If not all coordinates are written we have to extract them
        // to a buffer instead of passing the pointer to original coords
        if n != st.coords.len() {
            // Fill the buffer
            buf.reserve(n);
            for ind in index {
                buf.push(st.coords[ind]);
            }
            // Reset the pointer to buffer
            coord_ptr = buf.as_ptr();
        }

        let ok = unsafe {
            write_xtc(
                self.handle,
                n as i32,
                self.step,
                st.time,
                box_.as_ptr().cast::<rvec>(),
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

        if fr<self.frame_range.0 || fr>self.frame_range.1 {
            bail!("Can't seek to frame {}, allowed range is {:?}",fr,self.frame_range);
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

        if t<self.time_range.0 || t>self.time_range.1 {
            bail!("Can't seek to time {}, allowed range is {:?}",t,self.time_range);
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
        let ret = unsafe{
            xdr_xtc_seek_time(t,self.handle,self.natoms.try_into()?,false)
        };
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

    fn tell_first(&self) -> Result<(usize,f32)> {
        if !self.is_random_access {
            bail!("Random access is not possible for this XTC file!");
        }
        Ok((self.frame_range.0,self.time_range.0))
    }

    fn tell_last(&self) -> Result<(usize,f32)> {
        if !self.is_random_access {
            bail!("Random access is not possible for this XTC file!");
        }
        Ok((self.frame_range.1,self.time_range.1))
    }
}