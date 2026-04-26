use std::{ffi::CString, path::Path, sync::Arc};

use molar_gromacs::{CptHandle, TprPlugin};
use thiserror::Error;

use crate::prelude::*;

pub struct CptFileHandler {
    plugin: Arc<TprPlugin>,
    handle: *mut CptHandle,
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
        let plugin =
            TprPlugin::get_cached().map_err(|e| CptHandlerError::GromacsNotFound(e.to_string()))?;

        let c_path =
            CString::new(fname.as_ref().to_str().unwrap()).map_err(CptHandlerError::CStringNull)?;
        let handle = unsafe { (plugin.cpt.cpt_open)(c_path.as_ptr()) };

        if handle.is_null() {
            let msg = plugin.last_error();
            return Err(CptHandlerError::OpenFailed(msg).into());
        }

        Ok(CptFileHandler {
            plugin,
            handle,
            already_read: false,
        })
    }

    fn read_state(&mut self) -> Result<State, FileFormatError> {
        self.read_state_pick(true, true, false)
    }

    fn read_state_pick(
        &mut self,
        read_coords: bool,
        read_vels: bool,
        read_forces: bool,
    ) -> Result<State, FileFormatError> {
        if self.already_read {
            return Err(FileFormatError::Eof);
        }
        self.already_read = true;

        let cpt = &self.plugin.cpt;
        let h = self.handle;

        let natoms = unsafe { (cpt.cpt_natoms)(h) };
        let time = unsafe { (cpt.cpt_time)(h) };

        let mut st = State::default();
        // CPT plugin C ABI exposes f32; cast at boundary.
        st.time = time as Float;

        // Coordinates
        if read_coords {
            let mut buf: Vec<f32> = vec![0.0; natoms * 3];
            unsafe { (cpt.cpt_fill_coords)(h, buf.as_mut_ptr()) };
            st.coords.reserve(natoms);
            for i in 0..natoms {
                st.coords.push(Pos::new(
                    buf[i * 3]     as Float,
                    buf[i * 3 + 1] as Float,
                    buf[i * 3 + 2] as Float,
                ));
            }
        }

        // Box
        let mut box_buf = std::mem::MaybeUninit::<[f32; 9]>::uninit();
        unsafe { (cpt.cpt_fill_box)(h, box_buf.as_mut_ptr() as *mut f32) };
        let box_buf = unsafe { box_buf.assume_init() };
        let m = Matrix3f::from_iterator(box_buf.iter().map(|x| *x as Float));
        st.pbox = Some(PeriodicBox::from_matrix(m).map_err(CptHandlerError::Pbc)?);

        // Velocities
        if read_vels {
            let has_v = unsafe { (cpt.cpt_has_velocities)(h) };
            if has_v != 0 {
                let mut buf: Vec<f32> = vec![0.0; natoms * 3];
                unsafe { (cpt.cpt_fill_velocities)(h, buf.as_mut_ptr()) };
                let vels: Vec<Vel> = (0..natoms)
                    .map(|i| Vel::new(
                        buf[i * 3]     as Float,
                        buf[i * 3 + 1] as Float,
                        buf[i * 3 + 2] as Float,
                    ))
                    .collect();
                st.velocities = vels;
            } else {
                return Err(FileFormatError::NoVelocities);
            }
        }

        // Forces
        if read_forces {
            let has_f = unsafe { (cpt.cpt_has_forces)(h) };
            if has_f != 0 {
                let mut buf: Vec<f32> = vec![0.0; natoms * 3];
                unsafe { (cpt.cpt_fill_forces)(h, buf.as_mut_ptr()) };
                let force_data: Vec<Force> = (0..natoms)
                    .map(|i| Force::new(
                        buf[i * 3]     as Float,
                        buf[i * 3 + 1] as Float,
                        buf[i * 3 + 2] as Float,
                    ))
                    .collect();
                st.forces = force_data;
            } else {
                return Err(FileFormatError::NoForces);
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
                let is_not_found =
                    std::iter::successors(Some(&e as &dyn std::error::Error), |e| e.source())
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
        assert!(
            (st.time - 100000.0).abs() < 1.0,
            "time mismatch: {}",
            st.time
        );

        // First atom coords from gmx dump: x[0..2] = 7.46414, 4.04902, 8.06754 (nm)
        let p = &st.coords[0];
        assert!((p.x - 7.46414).abs() < 1e-4, "atom0 x: {}", p.x);
        assert!((p.y - 4.04902).abs() < 1e-4, "atom0 y: {}", p.y);
        assert!((p.z - 8.06754).abs() < 1e-4, "atom0 z: {}", p.z);

        // Box diagonal from gmx dump:
        //   box[0] = {9.64104, 0, 0}
        //   box[1] = {4.82052, 8.34932, 0}
        //   box[2] = {0, 0, 11.4521}
        let m = st.pbox.as_ref().unwrap().get_matrix();
        assert!(
            (m[(0, 0)] - 9.64104).abs() < 1e-4,
            "box[0][0]: {}",
            m[(0, 0)]
        );
        assert!(
            (m[(1, 1)] - 8.34932).abs() < 1e-4,
            "box[1][1]: {}",
            m[(1, 1)]
        );
        assert!(
            (m[(2, 2)] - 11.4521).abs() < 1e-3,
            "box[2][2]: {}",
            m[(2, 2)]
        );

        // Velocities from `gmx dump -cp tests/state.cpt`: v[0..2] = 0.816909, -0.184407, 0.448161 (nm/ps)
        assert!(!st.velocities.is_empty(), "CPT should have velocities");
        assert_eq!(st.velocities.len(), 96027, "velocity count mismatch");
        let v0 = &st.velocities[0];
        assert!((v0.x - 0.816909).abs() < 1e-4, "atom0 vx: {}", v0.x);
        assert!((v0.y - (-0.184407)).abs() < 1e-4, "atom0 vy: {}", v0.y);
        assert!((v0.z - 0.448161).abs() < 1e-4, "atom0 vz: {}", v0.z);

        // Forces are not present in this CPT file
        assert!(st.forces.is_empty(), "CPT should have no forces");
    }

    #[test]
    fn test_cpt_pick() {
        use crate::io::FileHandler;

        let mut h = match FileHandler::open("tests/state.cpt") {
            Ok(h) => h,
            Err(e) => {
                let is_not_found =
                    std::iter::successors(Some(&e as &dyn std::error::Error), |e| e.source())
                        .any(|e| e.to_string().contains("plugin not found"));
                if is_not_found {
                    eprintln!("Skipping test_cpt_pick: Gromacs plugin not available");
                    return;
                }
                panic!("unexpected error: {e}");
            }
        };

        // Read coords only — velocities should be empty even though the file has them
        let st = h.read_state_pick(true, false, false).unwrap();
        assert_eq!(st.len(), 96027);
        assert!(
            st.velocities.is_empty(),
            "velocities should be empty when not requested"
        );
        assert!(st.forces.is_empty());
    }
}
