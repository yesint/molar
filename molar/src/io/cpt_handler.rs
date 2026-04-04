use std::{ffi::CString, path::Path, sync::Arc};

use molar_gromacs::{CptHandle, TprPlugin};
use nalgebra::Matrix3;
use thiserror::Error;

use crate::prelude::*;

pub struct CptFileHandler {
    plugin:       Arc<TprPlugin>,
    handle:       *mut CptHandle,
    already_read: bool,
}

// CptHandle is heap-allocated C++ data managed through the plugin functions.
unsafe impl Send for CptFileHandler {}

#[derive(Debug, Error)]
pub enum CptHandlerError {
    #[error("Gromacs plugin not found (is MOLAR_GROMACS_PLUGIN set correctly?): {0}")]
    GromacsNotFound(String),

    #[error("failed to open CPT file: {0}")]
    OpenFailed(String),

    #[error("unexpected null character in path")]
    CStringNull(#[from] std::ffi::NulError),

    #[error("invalid periodic box")]
    Pbc(#[from] PeriodicBoxError),
}

impl Drop for CptFileHandler {
    fn drop(&mut self) {
        if !self.handle.is_null() {
            unsafe { (self.plugin.cpt.cpt_close)(self.handle) };
            self.handle = std::ptr::null_mut();
        }
    }
}

impl FileFormatHandler for CptFileHandler {
    fn open(fname: impl AsRef<Path>) -> Result<Self, FileFormatError>
    where
        Self: Sized,
    {
        let plugin = TprPlugin::get_cached()
            .map_err(|e| CptHandlerError::GromacsNotFound(e.to_string()))?;

        let c_path = CString::new(fname.as_ref().to_str().unwrap())
            .map_err(CptHandlerError::CStringNull)?;
        let handle = unsafe { (plugin.cpt.cpt_open)(c_path.as_ptr()) };

        if handle.is_null() {
            let msg = plugin.last_error();
            return Err(CptHandlerError::OpenFailed(msg).into());
        }

        Ok(CptFileHandler { plugin, handle, already_read: false })
    }

    fn read_state(&mut self) -> Result<State, FileFormatError> {
        if self.already_read {
            return Err(FileFormatError::Eof);
        }
        self.already_read = true;

        let cpt = &self.plugin.cpt;
        let h   = self.handle;

        let natoms = unsafe { (cpt.cpt_natoms)(h) };
        let time   = unsafe { (cpt.cpt_time)(h) };

        let mut coords_buf: Vec<f32>     = Vec::with_capacity(natoms * 3);
        let mut box_buf                  = std::mem::MaybeUninit::<[f32; 9]>::uninit();

        unsafe {
            coords_buf.set_len(natoms * 3);
            (cpt.cpt_fill_coords)(h, coords_buf.as_mut_ptr());
            (cpt.cpt_fill_box)(h, box_buf.as_mut_ptr() as *mut f32);
        }
        let box_buf = unsafe { box_buf.assume_init() };

        let mut st = State::default();
        st.time = time;
        st.coords.resize(natoms, Default::default());
        for i in 0..natoms {
            st.coords[i].coords.copy_from_slice(&coords_buf[i * 3..i * 3 + 3]);
        }

        let m = Matrix3::from_column_slice(&box_buf);
        st.pbox = Some(PeriodicBox::from_matrix(m).map_err(CptHandlerError::Pbc)?);

        Ok(st)
    }
}

#[cfg(test)]
mod tests {
    use crate::io::CptFileHandler;
    use crate::prelude::*;

    #[test]
    fn test_cpt() {
        let mut h = match CptFileHandler::open("tests/state.cpt") {
            Ok(h) => h,
            Err(e) => {
                let is_not_found = std::iter::successors(
                    Some(&e as &dyn std::error::Error),
                    |e| e.source(),
                )
                .any(|e| e.to_string().contains("plugin not found"));
                if is_not_found {
                    eprintln!("Skipping test_cpt: Gromacs plugin not available");
                    return;
                }
                panic!("unexpected error: {e}");
            }
        };

        let st = h.read_state().unwrap();

        // natoms and time from `gmx dump -cp tests/state.cpt`
        assert_eq!(st.len(), 96027);
        assert!((st.time - 100000.0_f32).abs() < 1.0, "time mismatch: {}", st.time);

        // First atom coords from gmx dump: x[0..2] = 7.46414, 4.04902, 8.06754 (nm)
        let p = &st.coords[0];
        assert!((p.x - 7.46414_f32).abs() < 1e-4, "atom0 x: {}", p.x);
        assert!((p.y - 4.04902_f32).abs() < 1e-4, "atom0 y: {}", p.y);
        assert!((p.z - 8.06754_f32).abs() < 1e-4, "atom0 z: {}", p.z);

        // Box diagonal from gmx dump:
        //   box[0] = {9.64104, 0, 0}
        //   box[1] = {4.82052, 8.34932, 0}
        //   box[2] = {0, 0, 11.4521}
        let m = st.pbox.unwrap().get_matrix();
        assert!((m[(0, 0)] - 9.64104_f32).abs() < 1e-4, "box[0][0]: {}", m[(0, 0)]);
        assert!((m[(1, 1)] - 8.34932_f32).abs() < 1e-4, "box[1][1]: {}", m[(1, 1)]);
        assert!((m[(2, 2)] - 11.4521_f32).abs()  < 1e-3, "box[2][2]: {}", m[(2, 2)]);
    }
}
