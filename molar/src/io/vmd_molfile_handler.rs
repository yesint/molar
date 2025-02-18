use super::{StateProvider, TopologyProvider};
use crate::core::*;
use molar_molfile::molfile_bindings::*;
use std::default::Default;
use std::ffi::{c_void, CStr, CString, NulError};
use std::ptr::{self, null_mut};
use std::str::Utf8Error;
use thiserror::Error;

#[derive(PartialEq)]
enum OpenMode {
    Read,
    Write,
    None,
}

pub enum VmdMolFileType {
    Pdb,
    Xyz,
    Dcd,
}

pub struct VmdMolFileHandler {
    // Plugin pointer
    plugin: *mut molfile_plugin_t,
    // File handle for C
    file_handle: *mut c_void,
    /// C string for defered file opening
    fname: CString,
    // Number of atoms for defered file opening
    natoms: usize,
    // Mode it which file is open
    mode: OpenMode,
}

// Allow sending handler between threads
unsafe impl Send for VmdMolFileHandler {}

#[derive(Error, Debug)]
pub enum VmdHandlerError {
    #[error("unexpected null characted")]
    CStringNull(#[from] NulError),

    #[error("invalid utf8")]
    CStringUtf8(#[from] Utf8Error),

    #[error("can't get plugin pointer")]
    NullPluginPtr,

    #[error("plugin can't open file for reading")]
    OpenRead,

    #[error("plugin can't open file for writing")]
    OpenWrite,

    #[error("writing {0} atoms to file opened for {1}")]
    NatomsMismatch(usize, usize),

    #[error("plugin failed to read structure")]
    ReadStructure,

    #[error("plugin failed to write structure")]
    WriteStructure,

    #[error("plugin failed to read state")]
    ReadState,

    #[error("plugin failed to write state")]
    WriteState,

    #[error("invalid periodic box")]
    Pbc(#[from] PeriodicBoxError),

    #[error("fixed size field is {0} while needed {1}")]
    FixedSizeFieldOverflow(usize, usize),
}

// Helper convertion function from C fixed-size string to String
fn char_slice_to_str(buf: &[::std::os::raw::c_char]) -> Result<String, VmdHandlerError> {
    unsafe { Ok(CStr::from_ptr(buf.as_ptr()).to_str()?.to_owned()) }
}

/// Universal handler of different VMD molfile file formats
impl VmdMolFileHandler {
    fn new_for_type(ftype: VmdMolFileType) -> Result<Self, VmdHandlerError> {
        // Get plugin pointer
        // C funtion registers plugin on first call
        // and returns stored pointer on later invocations
        let plugin = unsafe {
            match ftype {
                VmdMolFileType::Pdb => pdb_get_plugin_ptr(),
                VmdMolFileType::Xyz => xyz_get_plugin_ptr(),
                VmdMolFileType::Dcd => dcd_get_plugin_ptr(),
            }
        };

        if plugin.is_null() {
            Err(VmdHandlerError::NullPluginPtr)
        } else {
            Ok(VmdMolFileHandler {
                plugin,
                file_handle: ptr::null_mut(),
                fname: CString::default(),
                natoms: 0,
                mode: OpenMode::None,
            })
        }
    }

    fn open_read(&mut self, fname: &str) -> Result<(), VmdHandlerError> {
        // Prepare c-strings for file opening
        let f_name = CString::new(fname)?;

        // Open file and get file pointer
        let mut n: i32 = 0;
        self.file_handle = unsafe {
            self.plugin.as_ref().unwrap().open_file_read.unwrap()(
                f_name.as_ptr(),
                c"".as_ptr(), // Pass empty file type
                &mut n,
            ) // Call function
        };
        self.natoms = n as usize;

        if self.file_handle == ptr::null_mut() {
            return Err(VmdHandlerError::OpenRead);
        }

        self.mode = OpenMode::Read;

        Ok(())
    }

    fn open_write_if_needed(&mut self, natoms: usize) -> Result<(), VmdHandlerError> {
        // If file hanlder is already set just return
        if self.mode == OpenMode::Write {
            return Ok(());
        }

        // Set number of atoms
        self.natoms = natoms;

        // Open file and get file handle
        self.file_handle = unsafe {
            self.plugin.as_ref().unwrap().open_file_write.unwrap()(
                self.fname.as_ptr(),
                c"".as_ptr(), // Pass empty file type
                self.natoms as i32,
            ) // Call function
        };

        if self.file_handle == ptr::null_mut() {
            return Err(VmdHandlerError::OpenWrite);
        }

        self.mode = OpenMode::Write;

        Ok(())
    }

    pub fn open(fname: &str, ftype: VmdMolFileType) -> Result<Self, VmdHandlerError> {
        let mut instance = Self::new_for_type(ftype)?;
        instance.open_read(fname)?;
        Ok(instance)
    }

    pub fn create(fname: &str, ftype: VmdMolFileType) -> Result<Self, VmdHandlerError> {
        let mut instance = Self::new_for_type(ftype)?;
        instance.fname = CString::new(fname)?;
        // We can't open for writing here because we don't know
        // the number of atoms to write yet. Defer it to
        // actual writing operation
        Ok(instance)
    }

    pub fn read_topology(&mut self) -> Result<Topology, VmdHandlerError> {
        let mut optflags: i32 = 0;
        // Prepare array of atoms
        let mut vmd_atoms = Vec::<molfile_atom_t>::with_capacity(self.natoms);

        let ret = unsafe {
            self.plugin.as_ref().unwrap().read_structure.unwrap()(
                self.file_handle,
                &mut optflags,
                vmd_atoms.as_mut_ptr(),
            )
        };

        if ret != MOLFILE_SUCCESS {
            return Err(VmdHandlerError::ReadStructure);
        }

        // C function populated the atoms, set the vector size for Rust
        unsafe { vmd_atoms.set_len(self.natoms) }

        // Convert to Structure
        let mut top: TopologyStorage = Default::default();
        top.atoms.reserve(self.natoms);

        for ref vmd_at in vmd_atoms {
            let mut at = Atom {
                name: char_slice_to_str(&vmd_at.name)?,
                resid: vmd_at.resid,
                resname: char_slice_to_str(&vmd_at.resname)?,
                chain: char_slice_to_str(&vmd_at.chain)?.chars().next().unwrap(),
                charge: vmd_at.charge,
                occupancy: vmd_at.occupancy,
                bfactor: vmd_at.bfactor,
                ..Default::default()
            };
            // See if the element number and mass are set
            // and guess if required
            if vmd_at.atomicnumber == 0 || vmd_at.mass == 0.0 {
                at.guess_element_and_mass_from_name();
            } else {
                at.atomic_number = vmd_at.atomicnumber as u8;
                at.mass = vmd_at.mass;
            }
            top.atoms.push(at);
        }

        // Assign resindexes
        let top: Topology = top.into();
        top.assign_resindex();

        Ok(top)
    }

    pub fn write_topology(&mut self, data: &impl TopologyProvider) -> Result<(), VmdHandlerError> {
        let n = data.num_atoms();
        // Open file if not yet opened
        self.open_write_if_needed(n)?;

        let mut vmd_atoms = Vec::<molfile_atom_t>::with_capacity(n);
        for at in data.iter_atoms() {
            let mut vmd_at = molfile_atom_t::default();
            copy_str_to_c_buffer(&at.name, &mut vmd_at.name)?;
            copy_str_to_c_buffer(&at.resname, &mut vmd_at.resname)?;
            vmd_at.resid = at.resid;
            vmd_at.chain[0] = at.chain as i8;
            vmd_at.chain[1] = '\0' as i8;
            vmd_at.occupancy = at.occupancy;
            vmd_at.bfactor = at.bfactor;
            vmd_at.mass = at.mass;
            vmd_at.charge = at.charge;
            vmd_at.atomicnumber = at.atomic_number as i32;
            vmd_atoms.push(vmd_at);
        }

        let flags: u32 = MOLFILE_OCCUPANCY
            | MOLFILE_BFACTOR
            | MOLFILE_ATOMICNUMBER
            | MOLFILE_CHARGE
            | MOLFILE_MASS;
        let ret = unsafe {
            self.plugin.as_ref().unwrap().write_structure.unwrap()(
                self.file_handle,
                flags as i32,
                vmd_atoms.as_ptr(),
            )
        };

        match ret {
            MOLFILE_SUCCESS => Ok(()),
            _ => Err(VmdHandlerError::WriteStructure),
        }
    }

    pub fn read_state(&mut self) -> Result<Option<State>, VmdHandlerError> {
        let mut state: StateStorage = Default::default();

        // Allocate storage for coordinates, but don't initialize them
        // This doesn't waste time for initialization, which will be overwritten anyway
        state.coords = Vec::with_capacity(self.natoms as usize);

        let mut ts = molfile_timestep_t {
            coords: state.coords.as_mut_ptr().cast::<f32>(), // Raw ptr to allocated storage
            velocities: ptr::null_mut(),                     // Don't read velocities
            A: 0.0,
            B: 0.0,
            C: 0.0,
            alpha: 0.0,
            beta: 0.0,
            gamma: 0.0,
            physical_time: 0.0,
        };

        // Read the time step
        let ret = unsafe {
            self.plugin.as_ref().unwrap().read_next_timestep.unwrap()(
                self.file_handle,
                self.natoms as i32,
                &mut ts,
            )
        };

        // In case of successfull read populate rust State
        if ret == MOLFILE_SUCCESS {
            // C function populated the coordinates, set the vector size for Rust
            unsafe { state.coords.set_len(self.natoms as usize) }
            // Convert the box
            state.pbox = Some(PeriodicBox::from_vectors_angles(
                ts.A * 0.1,
                ts.B * 0.1,
                ts.C * 0.1,
                ts.alpha,
                ts.beta,
                ts.gamma,
            )?);
            // time
            state.time = ts.physical_time as f32;
            // Convert to nm
            for c in state.coords.iter_mut() {
                c.coords *= 0.1;
            }
        }

        match ret {
            MOLFILE_SUCCESS => Ok(Some(state.into())),
            MOLFILE_EOF => Ok(None),
            _ => Err(VmdHandlerError::ReadState),
        }
    }

    pub fn write_state(&mut self, data: &impl StateProvider) -> Result<(), VmdHandlerError> {
        let n = data.num_coords();

        // Open file if not yet opened
        self.open_write_if_needed(n)?;

        // Buffer for coordinates allocated on heap
        let mut buf = Vec::from_iter(
            data.iter_pos()
                .map(|p| p.coords.iter().cloned())
                .flatten()
                .map(|el| el * 10.0),
        );

        // Periodic box
        let (box_vec, box_ang) = match data.get_box() {
            Some(b) => b.to_vectors_angles(),
            None => (Vector3f::zeros(), Vector3f::zeros()),
        };

        let ts = molfile_timestep_t {
            coords: buf.as_mut_ptr(), // MolFile API requires a mut ptr for some reason
            velocities: null_mut(),
            A: box_vec[0] * 10.0,
            B: box_vec[1] * 10.0,
            C: box_vec[2] * 10.0,
            alpha: box_ang[0],
            beta: box_ang[1],
            gamma: box_ang[2],
            physical_time: data.get_time() as f64,
        };

        let ret = unsafe {
            self.plugin.as_ref().unwrap().write_timestep.unwrap() // get function ptr
            (self.file_handle, &ts) // Call function
        };

        match ret {
            MOLFILE_SUCCESS => Ok(()),
            _ => Err(VmdHandlerError::WriteState),
        }
    }

    // pub fn get_file_name(&self) -> &str {
    //     &self.file_name
    // }
}

impl Drop for VmdMolFileHandler {
    fn drop(&mut self) {
        if self.file_handle != ptr::null_mut() {
            match self.mode {
                OpenMode::Read => unsafe {
                    self.plugin.as_ref().unwrap().close_file_read.unwrap()(self.file_handle)
                },
                OpenMode::Write => unsafe {
                    self.plugin.as_ref().unwrap().close_file_write.unwrap()(self.file_handle)
                },
                OpenMode::None => (),
            };
        }
    }
}

fn copy_str_to_c_buffer(st: &str, cbuf: &mut [i8]) -> Result<(), VmdHandlerError> {
    let n = st.len();
    if n + 1 >= cbuf.len() {
        return Err(VmdHandlerError::FixedSizeFieldOverflow(cbuf.len(), n + 1));
    }

    let bytes = st.as_bytes();
    for i in 0..n {
        cbuf[i] = bytes[i] as i8;
    }
    cbuf[n + 1] = '\0' as i8;
    Ok(())
}
