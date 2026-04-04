use std::{ffi::CString, path::Path, sync::Arc};

use molar_gromacs::{CptHandle, TprPlugin};
use nalgebra::Matrix3;
use thiserror::Error;

use crate::prelude::*;
use super::ReadWhat;

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
        self.read_state_pick(ReadWhat::all())
    }

    fn read_state_pick(&mut self, what: ReadWhat) -> Result<State, FileFormatError> {
        if self.already_read {
            return Err(FileFormatError::Eof);
        }
        self.already_read = true;

        let cpt = &self.plugin.cpt;
        let h   = self.handle;

        let natoms = unsafe { (cpt.cpt_natoms)(h) };
        let time   = unsafe { (cpt.cpt_time)(h) };

        let mut st = State::default();
        st.time = time;

        // Coordinates
        if what.contains(ReadWhat::COORDS) {
            let mut buf: Vec<f32> = vec![0.0; natoms * 3];
            unsafe { (cpt.cpt_fill_coords)(h, buf.as_mut_ptr()) };
            st.coords.resize(natoms, Default::default());
            for i in 0..natoms {
                st.coords[i].coords.copy_from_slice(&buf[i * 3..i * 3 + 3]);
            }
        }

        // Box
        let mut box_buf = std::mem::MaybeUninit::<[f32; 9]>::uninit();
        unsafe { (cpt.cpt_fill_box)(h, box_buf.as_mut_ptr() as *mut f32) };
        let box_buf = unsafe { box_buf.assume_init() };
        let m = Matrix3::from_column_slice(&box_buf);
        st.pbox = Some(PeriodicBox::from_matrix(m).map_err(CptHandlerError::Pbc)?);

        // Velocities
        if what.contains(ReadWhat::VELOCITIES) {
            let has_v = unsafe { (cpt.cpt_has_velocities)(h) };
            if has_v != 0 {
                let mut buf: Vec<f32> = vec![0.0; natoms * 3];
                unsafe { (cpt.cpt_fill_velocities)(h, buf.as_mut_ptr()) };
                let vels: Vec<Vel> = (0..natoms)
                    .map(|i| Vel::new(buf[i * 3], buf[i * 3 + 1], buf[i * 3 + 2]))
                    .collect();
                st.velocities = Some(vels);
            }
        }

        // Forces
        if what.contains(ReadWhat::FORCES) {
            let has_f = unsafe { (cpt.cpt_has_forces)(h) };
            if has_f != 0 {
                let mut buf: Vec<f32> = vec![0.0; natoms * 3];
                unsafe { (cpt.cpt_fill_forces)(h, buf.as_mut_ptr()) };
                let forces: Vec<Force> = (0..natoms)
                    .map(|i| Force::new(buf[i * 3], buf[i * 3 + 1], buf[i * 3 + 2]))
                    .collect();
                st.forces = Some(forces);
            }
        }

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
        let m = st.pbox.as_ref().unwrap().get_matrix();
        assert!((m[(0, 0)] - 9.64104_f32).abs() < 1e-4, "box[0][0]: {}", m[(0, 0)]);
        assert!((m[(1, 1)] - 8.34932_f32).abs() < 1e-4, "box[1][1]: {}", m[(1, 1)]);
        assert!((m[(2, 2)] - 11.4521_f32).abs()  < 1e-3, "box[2][2]: {}", m[(2, 2)]);

        // Velocities from `gmx dump -cp tests/state.cpt`: v[0..2] = 0.816909, -0.184407, 0.448161 (nm/ps)
        let vels = st.velocities.as_ref().expect("CPT should have velocities");
        assert_eq!(vels.len(), 96027, "velocity count mismatch");
        let v0 = &vels[0];
        assert!((v0.x - 0.816909_f32).abs() < 1e-4, "atom0 vx: {}", v0.x);
        assert!((v0.y - (-0.184407_f32)).abs() < 1e-4, "atom0 vy: {}", v0.y);
        assert!((v0.z - 0.448161_f32).abs() < 1e-4, "atom0 vz: {}", v0.z);

        // Forces are not present in this CPT file
        assert!(st.forces.is_none(), "CPT should have no forces");
    }

    #[test]
    fn test_cpt_pick() {
        use crate::io::{FileHandler, ReadWhat};

        let mut h = match FileHandler::open("tests/state.cpt") {
            Ok(h) => h,
            Err(e) => {
                let is_not_found = std::iter::successors(
                    Some(&e as &dyn std::error::Error),
                    |e| e.source(),
                )
                .any(|e| e.to_string().contains("plugin not found"));
                if is_not_found {
                    eprintln!("Skipping test_cpt_pick: Gromacs plugin not available");
                    return;
                }
                panic!("unexpected error: {e}");
            }
        };

        // Read coords only — velocities should be None even though the file has them
        let st = h.read_state_pick(ReadWhat::COORDS).unwrap();
        assert_eq!(st.len(), 96027);
        assert!(st.velocities.is_none(), "velocities should be None when not requested");
        assert!(st.forces.is_none());
    }
}
