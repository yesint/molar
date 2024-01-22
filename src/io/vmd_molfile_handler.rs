use super::{IoReader, IoTrajectoryReader, IoTrajectoryWriter, IoTopologyReader, IoTopologyWriter, IoWriter, IoTopologyProvider, IoIndexProvider, IoStateProvider, IoOnceReader, IoOnceWriter};
use crate::core::*;
use crate::io::get_ext;
use anyhow::{bail, Result};
use ascii::{AsciiStr, AsciiString};
use molar_molfile::molfile_bindings::*;
use std::default::Default;
use std::ffi::{c_void, CStr, CString};
use std::ptr::{self, null_mut};

enum OpenMode {
    Read,
    Write,
    None,
}

pub struct VmdMolFileHandler<'a> {
    // File
    file_name: String,
    // Plugin pointer
    plugin: &'a molfile_plugin_t,
    // File handle for C
    file_handle: *mut c_void,
    mode: OpenMode,
    natoms: usize,
}

// Helper convertion function from C fixed-size string to AsciiString
fn char_slice_to_ascii_str(buf: &[::std::os::raw::c_char]) -> AsciiString {
    let cstr = unsafe { CStr::from_ptr(buf.as_ptr()).to_bytes() };
    let s = unsafe { AsciiString::from_ascii_unchecked(cstr) };
    s
}

/// Universal handler of different VMD molfile file formats
impl VmdMolFileHandler<'_> {
    fn new(fname: &str) -> Result<Self> {
        // Get plugin pointer
        // C funtion registers plugin on first call
        // and returns stored pointer on later invocations
        let plugin = unsafe {
            let ext = get_ext(fname)?;
            match ext {
                "pdb" => pdb_get_plugin_ptr(),
                "xyz" => xyz_get_plugin_ptr(),
                "dcd" => dcd_get_plugin_ptr(),
                &_ => bail!("Unrecognized extention {ext}!"),
            }
            .as_ref()
            .unwrap()
        };

        Ok(VmdMolFileHandler {
            file_name: fname.to_owned(),
            plugin,
            file_handle: ptr::null_mut(),
            mode: OpenMode::None,
            natoms: 0,
        })
    }

    fn open_read(&mut self) -> Result<()> {
        // Prepare c-strings for file opening
        let f_type = CString::new("")?; // Not used
        let f_name = CString::new(self.file_name.clone())?;

        // Open file and get file pointer
        let mut n: i32 = 0;
        self.file_handle = unsafe {
            self.plugin.open_file_read.unwrap() // get function ptr
            (f_name.as_ptr(), f_type.as_ptr(), &mut n) // Call function
        };
        self.natoms = n as usize;

        if self.file_handle == ptr::null_mut() {
            bail!("Can't open file {} for reading!", self.file_name);
        }

        self.mode = OpenMode::Read;
        Ok(())
    }

    fn open_write(&mut self, natoms: usize) -> Result<()> {
        // Prepare c-strings for file opening
        let f_type = CString::new("")?; // Not used
        let f_name = CString::new(self.file_name.clone())?;

        // Set number of atoms
        self.natoms = natoms;

        // Open file and get file handle
        self.file_handle = unsafe {
            self.plugin.open_file_write.unwrap()(
                f_name.as_ptr(),
                f_type.as_ptr(),
                self.natoms as i32,
            ) // Call function
        };

        if self.file_handle == ptr::null_mut() {
            bail!("Can't open file {} for writing!", self.file_name);
        }

        self.mode = OpenMode::Write;

        Ok(())
    }

    fn try_open_write(&mut self, natoms: usize) -> Result<()> {
        match self.mode {
            OpenMode::None => self.open_write(natoms),
            OpenMode::Write => {
                if natoms != self.natoms {
                    bail!(
                        "Number of atoms mismatch: given {}, file opened for writing with {}!",
                        natoms,
                        self.natoms
                    )
                } else {
                    Ok(())
                }
            }
            OpenMode::Read => unreachable!(),
        }
    }
}

impl IoReader for VmdMolFileHandler<'_> {
    fn open(fname: &str) -> Result<Self> {
        let mut instance = Self::new(fname)?;
        instance.open_read()?;
        Ok(instance)
    }
}

impl IoWriter for VmdMolFileHandler<'_> {
    fn create(fname: &str) -> Result<Self> {
        let instance = Self::new(fname)?;
        // We can't open for writing here because we don't know
        // the number of atoms to write yet. Defer it to
        // actual writing operation
        Ok(instance)
    }
}

impl IoTopologyReader for VmdMolFileHandler<'_> {
    fn read_topology(&mut self) -> Result<Topology> {
        let mut optflags: i32 = 0;
        // Prepare array of atoms
        let mut vmd_atoms = Vec::<molfile_atom_t>::with_capacity(self.natoms);

        let ret = unsafe {
            self.plugin.read_structure.unwrap()(
                self.file_handle,
                &mut optflags,
                vmd_atoms.as_mut_ptr(),
            )
        };

        if ret != MOLFILE_SUCCESS {
            bail!("Plugin error reading structure!");
        }

        // C function populated the atoms, set the vector size for Rust
        unsafe { vmd_atoms.set_len(self.natoms) }

        // Convert to Structure
        let mut topology: Topology = Default::default();
        topology.atoms.reserve(self.natoms);

        for ref vmd_at in vmd_atoms {
            let mut at = Atom {
                name: char_slice_to_ascii_str(&vmd_at.name),
                resid: vmd_at.resid,
                resname: char_slice_to_ascii_str(&vmd_at.resname),
                chain: char_slice_to_ascii_str(&vmd_at.chain).first().unwrap(),
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
            topology.atoms.push(at);
        }

        // Assign resindexes
        topology.assign_resindex();

        Ok(topology)
    }
}

fn copy_str_to_c_buffer(st: &AsciiStr, cbuf: &mut [i8]) {
    let n = st.len();
    if n + 1 >= cbuf.len() {
        panic!("VMD fixed size field is too short!");
    }
    for i in 0..n {
        cbuf[i] = st[i] as i8;
    }
    cbuf[n + 1] = '\0' as i8;
}

impl IoTopologyWriter for VmdMolFileHandler<'_> {
    fn write_topology(
        &mut self,
        data: &(impl IoIndexProvider+IoTopologyProvider),
    ) -> Result<()> {
        let index = data.get_index();
        let top = data.get_topology();
        let n = index.len();
        // Open file if not yet opened
        self.try_open_write(n)?;

        let mut vmd_atoms = Vec::<molfile_atom_t>::with_capacity(n);
        for ind in index {
            let at = &top.atoms[ind];
            let mut vmd_at = molfile_atom_t::default();
            copy_str_to_c_buffer(&at.name, &mut vmd_at.name);
            copy_str_to_c_buffer(&at.resname, &mut vmd_at.resname);
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
            self.plugin.write_structure.unwrap()(self.file_handle, flags as i32, vmd_atoms.as_ptr())
        };

        match ret {
            MOLFILE_SUCCESS => Ok(()),
            _ => bail!("Error writing structure"),
        }
    }
}

impl IoTrajectoryReader for VmdMolFileHandler<'_> {
    fn read_state(&mut self) -> Result<Option<State>> {
        let mut state: State = Default::default();

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
            self.plugin.read_next_timestep.unwrap()(self.file_handle, self.natoms as i32, &mut ts)
        };

        // In case of successfull read populate rust State
        if ret == MOLFILE_SUCCESS {
            // C function populated the coordinates, set the vector size for Rust
            unsafe { state.coords.set_len(self.natoms as usize) }
            // Convert the box
            state.box_ = PeriodicBox::from_vectors_angles(
                ts.A * 0.1,
                ts.B * 0.1,
                ts.C * 0.1,
                ts.alpha,
                ts.beta,
                ts.gamma,
            )
            .ok();
            // time
            state.time = ts.physical_time as f32;
            // Convert to nm
            for c in state.coords.iter_mut() {
                c.coords *= 0.1;
            }
        }

        match ret {
            MOLFILE_SUCCESS => Ok(Some(state)),
            MOLFILE_EOF => Ok(None),
            _ => bail!("Error reading timestep!"),
        }
    }
}

impl IoTrajectoryWriter for VmdMolFileHandler<'_> {
    fn write_state(
        &mut self,
        data: &(impl IoIndexProvider+IoStateProvider),
    ) -> Result<()> {
        let index = data.get_index();
        let st = data.get_state();
        let n = index.len();

        // Open file if not yet opened
        self.try_open_write(n)?;

        // Buffer for coordinates allocated on heap
        let mut buf = Vec::<f32>::with_capacity(3 * n);
        // Fill the buffer and convert to angstroms
        for ind in index {
            for dim in 0..3 {
                buf.push(st.coords[ind][dim] * 10.0);
            }
        }

        // Periodic box
        let (box_vec, box_ang) = match st.box_.as_ref() {
            Some(b) => b.to_vectors_angles(),
            None => (Vector3f::zeros(), Vector3f::zeros()),
        };

        let ts = molfile_timestep_t {
            coords: buf.as_mut_ptr(),
            velocities: null_mut(),
            A: box_vec[0] * 10.0,
            B: box_vec[1] * 10.0,
            C: box_vec[2] * 10.0,
            alpha: box_ang[0],
            beta: box_ang[1],
            gamma: box_ang[2],
            physical_time: st.time as f64,
        };

        let ret = unsafe {
            self.plugin.write_timestep.unwrap() // get function ptr
            (self.file_handle, &ts)
        };

        match ret {
            MOLFILE_SUCCESS => Ok(()),
            _ => bail!("Error writing timestep!"),
        }
    }
}

impl Drop for VmdMolFileHandler<'_> {
    fn drop(&mut self) {
        if self.file_handle != ptr::null_mut() {
            match self.mode {
                OpenMode::Read => unsafe { self.plugin.close_file_read.unwrap()(self.file_handle) },
                OpenMode::Write => unsafe {
                    self.plugin.close_file_write.unwrap()(self.file_handle)
                },
                OpenMode::None => (),
            };
        }
    }
}

impl IoOnceReader for VmdMolFileHandler<'_> {
    fn read(&mut self) -> Result<(Topology,State)> {
        // Read topology and first frame at once
        let top = self.read_topology()?;
        let st = self.read_state()?.unwrap();
        Ok((top,st))
    }
}

impl IoOnceWriter for VmdMolFileHandler<'_> {
    fn write(&mut self, data: &(impl IoIndexProvider + IoTopologyProvider + IoStateProvider)) -> Result<()> {
        // Write topology and first frame at once
        self.write_topology(data)?;
        self.write_state(data)?;
        Ok(())
    }
}