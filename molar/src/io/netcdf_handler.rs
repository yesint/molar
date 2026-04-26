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

        fn create(_fname: impl AsRef<Path>) -> Result<Self, super::super::FileFormatError>
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
    use std::path::{Path, PathBuf};
    use thiserror::Error;

    use crate::io::{FileFormatError, FileFormatHandler, SaveState};
    use crate::periodic_box::{PeriodicBox, PeriodicBoxError};
    use crate::{Float, Pos, State};

    //------------------------------------------------------------------
    // Sub-types
    //------------------------------------------------------------------

    pub(crate) struct NetCdfReader {
        file: netcdf::File,
        n_atoms: usize,
        n_frames: usize,
        cur_frame: usize,
        has_time: bool,
        has_box: bool,
    }

    pub(crate) struct NetCdfWriter {
        path: PathBuf,
        /// `None` until the first `write_state` call (lazy init, needs atom count).
        file: Option<netcdf::FileMut>,
        n_atoms: usize,
        cur_frame: usize,
    }

    pub enum NetCdfFileHandler {
        Reader(NetCdfReader),
        Writer(NetCdfWriter),
    }

    //------------------------------------------------------------------
    // Error type
    //------------------------------------------------------------------

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

    //------------------------------------------------------------------
    // Write helper: create a properly structured AMBER NCTRAJ file
    //------------------------------------------------------------------

    fn init_file(path: &Path, n_atoms: usize) -> Result<netcdf::FileMut, NetCdfHandlerError> {
        let mut file = netcdf::create(path).map_err(NetCdfHandlerError::Nc)?;

        // Global attributes (FileMut::add_attribute)
        file.add_attribute("Conventions", "AMBER")
            .map_err(NetCdfHandlerError::Nc)?;
        file.add_attribute("ConventionVersion", "1.0")
            .map_err(NetCdfHandlerError::Nc)?;
        file.add_attribute("program", "molar")
            .map_err(NetCdfHandlerError::Nc)?;
        file.add_attribute("programVersion", env!("CARGO_PKG_VERSION"))
            .map_err(NetCdfHandlerError::Nc)?;
        file.add_attribute("title", "molar generated trajectory")
            .map_err(NetCdfHandlerError::Nc)?;

        // Dimensions
        file.add_unlimited_dimension("frame")
            .map_err(NetCdfHandlerError::Nc)?;
        file.add_dimension("atom", n_atoms)
            .map_err(NetCdfHandlerError::Nc)?;
        file.add_dimension("spatial", 3)
            .map_err(NetCdfHandlerError::Nc)?;
        file.add_dimension("cell_spatial", 3)
            .map_err(NetCdfHandlerError::Nc)?;
        file.add_dimension("cell_angular", 3)
            .map_err(NetCdfHandlerError::Nc)?;

        // Variables
        file.add_variable::<f32>("coordinates", &["frame", "atom", "spatial"])
            .map_err(NetCdfHandlerError::Nc)?;
        file.variable_mut("coordinates")
            .unwrap()
            .put_attribute("units", "angstrom")
            .map_err(NetCdfHandlerError::Nc)?;

        file.add_variable::<f32>("time", &["frame"])
            .map_err(NetCdfHandlerError::Nc)?;
        file.variable_mut("time")
            .unwrap()
            .put_attribute("units", "picosecond")
            .map_err(NetCdfHandlerError::Nc)?;

        file.add_variable::<f64>("cell_lengths", &["frame", "cell_spatial"])
            .map_err(NetCdfHandlerError::Nc)?;
        file.variable_mut("cell_lengths")
            .unwrap()
            .put_attribute("units", "angstrom")
            .map_err(NetCdfHandlerError::Nc)?;

        file.add_variable::<f64>("cell_angles", &["frame", "cell_angular"])
            .map_err(NetCdfHandlerError::Nc)?;
        file.variable_mut("cell_angles")
            .unwrap()
            .put_attribute("units", "degree")
            .map_err(NetCdfHandlerError::Nc)?;

        Ok(file)
    }

    //------------------------------------------------------------------
    // FileFormatHandler implementation
    //------------------------------------------------------------------

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

            if file.variable("coordinates").is_none() {
                return Err(NetCdfHandlerError::MissingCoordinates.into());
            }

            let has_time = file.variable("time").is_some();
            let has_box = file.variable("cell_lengths").is_some()
                && file.variable("cell_angles").is_some();

            Ok(NetCdfFileHandler::Reader(NetCdfReader {
                file,
                n_atoms,
                n_frames,
                cur_frame: 0,
                has_time,
                has_box,
            }))
        }

        fn create(fname: impl AsRef<Path>) -> Result<Self, FileFormatError>
        where
            Self: Sized,
        {
            Ok(NetCdfFileHandler::Writer(NetCdfWriter {
                path: fname.as_ref().to_path_buf(),
                file: None,
                n_atoms: 0,
                cur_frame: 0,
            }))
        }

        fn read_state(&mut self) -> Result<State, FileFormatError> {
            match self {
                Self::Reader(r) => {
                    if r.cur_frame >= r.n_frames {
                        return Err(FileFormatError::Eof);
                    }

                    let fr = r.cur_frame;

                    // Read coordinates: (frame, atom, spatial) in Ångström → nm
                    let coords_var = r.file.variable("coordinates").unwrap();
                    let raw: Vec<f32> = coords_var
                        .get_values(([fr, 0, 0], [1, r.n_atoms, 3]))
                        .map_err(NetCdfHandlerError::Nc)?;

                    // NetCDF coords are f32 on disk; cast to Float at the boundary.
                    let coords: Vec<Pos> = raw
                        .chunks_exact(3)
                        .map(|c| Pos::new(c[0] as Float * 0.1, c[1] as Float * 0.1, c[2] as Float * 0.1))
                        .collect();

                    // Read time (ps)
                    let time = if r.has_time {
                        let t_var = r.file.variable("time").unwrap();
                        t_var
                            .get_value::<f32, _>([fr])
                            .map_err(NetCdfHandlerError::Nc)? as Float
                    } else {
                        fr as Float
                    };

                    // Read periodic box (f64 on disk)
                    let pbox = if r.has_box {
                        let len_var = r.file.variable("cell_lengths").unwrap();
                        let ang_var = r.file.variable("cell_angles").unwrap();

                        let lengths: Vec<f64> = len_var
                            .get_values(([fr, 0], [1, 3]))
                            .map_err(NetCdfHandlerError::Nc)?;
                        let angles: Vec<f64> = ang_var
                            .get_values(([fr, 0], [1, 3]))
                            .map_err(NetCdfHandlerError::Nc)?;

                        Some(
                            PeriodicBox::from_vectors_angles(
                                (lengths[0] * 0.1) as Float,
                                (lengths[1] * 0.1) as Float,
                                (lengths[2] * 0.1) as Float,
                                angles[0] as Float,
                                angles[1] as Float,
                                angles[2] as Float,
                            )
                            .map_err(NetCdfHandlerError::Pbc)?,
                        )
                    } else {
                        None
                    };

                    r.cur_frame += 1;
                    Ok(State { coords, time, pbox, ..Default::default() })
                }
                Self::Writer(_) => unreachable!(),
            }
        }

        fn write_state(&mut self, data: &dyn SaveState) -> Result<(), FileFormatError> {
            match self {
                Self::Writer(w) => {
                    // Lazy init: create the file on the first write when atom count is known
                    if w.file.is_none() {
                        w.n_atoms = data.len();
                        w.file = Some(init_file(&w.path, w.n_atoms).map_err(FileFormatError::from)?);
                    }

                    let fr = w.cur_frame;
                    let n_atoms = w.n_atoms;
                    let file = w.file.as_mut().unwrap();

                    // Write coordinates: nm → Å (NetCDF coord variable is f32 on disk)
                    let raw: Vec<f32> = data
                        .iter_pos_dyn()
                        .flat_map(|p| [(p.x * 10.0) as f32, (p.y * 10.0) as f32, (p.z * 10.0) as f32])
                        .collect();
                    file.variable_mut("coordinates")
                        .unwrap()
                        .put_values(&raw, ([fr, 0, 0], [1, n_atoms, 3]))
                        .map_err(NetCdfHandlerError::Nc)?;

                    // Write time (ps) — netcdf variable is f32 on disk.
                    file.variable_mut("time")
                        .unwrap()
                        .put_value(data.get_time() as f32, [fr])
                        .map_err(NetCdfHandlerError::Nc)?;

                    // Write periodic box if present
                    if let Some(b) = data.get_box() {
                        let (lengths, angles) = b.to_vectors_angles();
                        let cell_lengths: Vec<f64> =
                            lengths.iter().map(|&l| (l * 10.0) as f64).collect();
                        let cell_angles: Vec<f64> =
                            angles.iter().map(|&a| a as f64).collect();

                        file.variable_mut("cell_lengths")
                            .unwrap()
                            .put_values(&cell_lengths, ([fr, 0], [1, 3]))
                            .map_err(NetCdfHandlerError::Nc)?;
                        file.variable_mut("cell_angles")
                            .unwrap()
                            .put_values(&cell_angles, ([fr, 0], [1, 3]))
                            .map_err(NetCdfHandlerError::Nc)?;
                    }

                    w.cur_frame += 1;
                    Ok(())
                }
                Self::Reader(_) => unreachable!(),
            }
        }

        fn seek_frame(&mut self, fr: usize) -> Result<(), FileFormatError> {
            match self {
                Self::Reader(r) => {
                    if fr >= r.n_frames {
                        return Err(FileFormatError::SeekFrame(fr));
                    }
                    r.cur_frame = fr;
                    Ok(())
                }
                Self::Writer(_) => Err(FileFormatError::NotRandomAccessFormat),
            }
        }

        fn seek_time(&mut self, t: Float) -> Result<(), FileFormatError> {
            match self {
                Self::Reader(r) => {
                    if !r.has_time {
                        return Err(FileFormatError::NotRandomAccessFormat);
                    }
                    // Read all time stamps at once — cheap since it's a 1-D variable.
                    let times: Vec<f32> = r
                        .file
                        .variable("time")
                        .unwrap()
                        .get_values(..)
                        .map_err(NetCdfHandlerError::Nc)?;
                    // Find the last frame whose time is strictly less than `t`,
                    // matching the serial-read fallback semantics in FileHandler.
                    let t_disk = t as f32;
                    let fr = times.partition_point(|&ts| ts < t_disk);
                    if fr >= r.n_frames {
                        return Err(FileFormatError::SeekTime(t));
                    }
                    r.cur_frame = fr;
                    Ok(())
                }
                Self::Writer(_) => Err(FileFormatError::NotRandomAccessFormat),
            }
        }

        fn seek_last(&mut self) -> Result<(), FileFormatError> {
            match self {
                Self::Reader(r) => {
                    if r.n_frames == 0 {
                        return Err(FileFormatError::Eof);
                    }
                    r.cur_frame = r.n_frames - 1;
                    Ok(())
                }
                Self::Writer(_) => Err(FileFormatError::NotRandomAccessFormat),
            }
        }
    }
}
