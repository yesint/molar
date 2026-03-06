#[cfg(not(feature = "netcdf"))]
pub use disabled::*;
#[cfg(feature = "netcdf")]
pub use enabled::*;

//----------------------------------------------------------------------
// Stub used when the `netcdf` feature is not enabled
//----------------------------------------------------------------------
#[cfg(not(feature = "netcdf"))]
mod disabled {
    use std::path::Path;
    use thiserror::Error;

    use crate::io::FileFormatHandler;

    pub struct NetCdfFileHandler {}

    #[derive(Error, Debug)]
    pub enum NetCdfHandlerError {
        #[error("molar was compiled without NetCDF support (enable the `netcdf` feature)")]
        Disabled,
    }

    impl FileFormatHandler for NetCdfFileHandler {
        fn open(_fname: impl AsRef<Path>) -> Result<Self, super::super::FileFormatError>
        where
            Self: Sized,
        {
            Err(NetCdfHandlerError::Disabled.into())
        }
    }
}

//----------------------------------------------------------------------
// Full implementation when the `netcdf` feature is enabled
//----------------------------------------------------------------------
#[cfg(feature = "netcdf")]
mod enabled {
    use std::path::Path;
    use thiserror::Error;

    use crate::io::{FileFormatError, FileFormatHandler};
    use crate::periodic_box::{PeriodicBox, PeriodicBoxError};
    use crate::{Pos, State};

    pub struct NetCdfFileHandler {
        file: netcdf::File,
        n_atoms: usize,
        n_frames: usize,
        cur_frame: usize,
        has_time: bool,
        has_box: bool,
    }

    #[derive(Error, Debug)]
    pub enum NetCdfHandlerError {
        #[error("netcdf error: {0}")]
        Nc(#[from] netcdf::Error),

        #[error("missing required dimension '{0}'")]
        MissingDimension(&'static str),

        #[error("missing required variable 'coordinates'")]
        MissingCoordinates,

        #[error("invalid periodic box")]
        Pbc(#[from] PeriodicBoxError),
    }

    impl FileFormatHandler for NetCdfFileHandler {
        fn open(fname: impl AsRef<Path>) -> Result<Self, FileFormatError>
        where
            Self: Sized,
        {
            let file = netcdf::open(fname).map_err(NetCdfHandlerError::Nc)?;

            let n_atoms = file
                .dimension("atom")
                .ok_or(NetCdfHandlerError::MissingDimension("atom"))?
                .len();

            let n_frames = file
                .dimension("frame")
                .ok_or(NetCdfHandlerError::MissingDimension("frame"))?
                .len();

            // Validate that coordinates variable exists
            if file.variable("coordinates").is_none() {
                return Err(NetCdfHandlerError::MissingCoordinates.into());
            }

            let has_time = file.variable("time").is_some();
            let has_box = file.variable("cell_lengths").is_some()
                && file.variable("cell_angles").is_some();

            Ok(NetCdfFileHandler {
                file,
                n_atoms,
                n_frames,
                cur_frame: 0,
                has_time,
                has_box,
            })
        }

        fn read_state(&mut self) -> Result<State, FileFormatError> {
            if self.cur_frame >= self.n_frames {
                return Err(FileFormatError::Eof);
            }

            let fr = self.cur_frame;

            // Read coordinates: (frame, atom, spatial) in Ångström → nm
            let coords_var = self.file.variable("coordinates").unwrap();
            let raw: Vec<f32> = coords_var
                .get_values(([fr, 0, 0], [1, self.n_atoms, 3]))
                .map_err(NetCdfHandlerError::Nc)?;

            let coords: Vec<Pos> = raw
                .chunks_exact(3)
                .map(|c| Pos::new(c[0] * 0.1, c[1] * 0.1, c[2] * 0.1))
                .collect();

            // Read time (ps)
            let time = if self.has_time {
                let t_var = self.file.variable("time").unwrap();
                t_var
                    .get_value::<f32, _>([fr])
                    .map_err(NetCdfHandlerError::Nc)?
            } else {
                fr as f32
            };

            // Read periodic box
            let pbox = if self.has_box {
                let len_var = self.file.variable("cell_lengths").unwrap();
                let ang_var = self.file.variable("cell_angles").unwrap();

                // cell_lengths in Å → nm; cell_angles in degrees
                let lengths: Vec<f64> = len_var
                    .get_values(([fr, 0], [1, 3]))
                    .map_err(NetCdfHandlerError::Nc)?;
                let angles: Vec<f64> = ang_var
                    .get_values(([fr, 0], [1, 3]))
                    .map_err(NetCdfHandlerError::Nc)?;

                Some(
                    PeriodicBox::from_vectors_angles(
                        (lengths[0] * 0.1) as f32,
                        (lengths[1] * 0.1) as f32,
                        (lengths[2] * 0.1) as f32,
                        angles[0] as f32,
                        angles[1] as f32,
                        angles[2] as f32,
                    )
                    .map_err(NetCdfHandlerError::Pbc)?,
                )
            } else {
                None
            };

            self.cur_frame += 1;
            Ok(State { coords, time, pbox })
        }

        fn seek_frame(&mut self, fr: usize) -> Result<(), FileFormatError> {
            if fr >= self.n_frames {
                return Err(FileFormatError::SeekFrame(fr));
            }
            self.cur_frame = fr;
            Ok(())
        }
    }
}
