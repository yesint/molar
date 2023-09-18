use super::{IoReader, IoWriter, IoStructureReader, IoStateReader, IoStructureWriter, IoStateWriter};
use crate::core::*;
use ascii::{AsciiString,AsciiStr};

use crate::io::get_ext;

use molar_molfile::molfile_bindings::{
    dcd_get_plugin_ptr, molfile_atom_t, molfile_plugin_t, molfile_timestep_t, 
    pdb_get_plugin_ptr, xyz_get_plugin_ptr, MOLFILE_EOF, MOLFILE_SUCCESS,
};

use std::ffi::{c_void, CStr, CString};
use std::ptr;
use std::iter::zip;

use std::default::Default;

use anyhow::{bail, Result};
use nalgebra::Matrix3;
//use num_traits::pow::Pow;

fn box_from_vmd(a: f32, b: f32, c: f32, alpha: f32, beta: f32, gamma: f32) -> Matrix3<f32> {
    let mut box_ = Matrix3::<f32>::zeros();

    box_[(0,0)] = a;

    if alpha != 90.0 || beta != 90.0 || gamma != 90.0 {
        let cosa = if alpha != 90.0 {
            alpha.to_radians().cos()
        } else {
            0.0
        };
        let cosb = if beta != 90.0 {
            beta.to_radians().cos()
        } else {
            0.0
        };
        let (sing, cosg) = if gamma != 90.0 {
            gamma.to_radians().sin_cos()
        } else {
            (1.0, 0.0)
        };
        box_[(0,1)] = b * cosg;
        box_[(1,1)] = b * sing;
        box_[(0,2)] = c * cosb;
        box_[(1,2)] = c * (cosa - cosb * cosg) / sing;
        box_[(2,2)] = (c * c - box_[(0,2)].powf(2.0) - box_[(1,2)].powf(2.0)).sqrt();
    } else {
        box_[(1,1)] = b;
        box_[(2,2)] = c;
    }
    box_ / 10.0 // Convert to nm
}


enum OpenMode {
    Read,
    Write,
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

// Helper convertion function from C-wrapped fixed-size string to AsciiString
fn char_slice_to_ascii_str(buf: &[::std::os::raw::c_char]) -> AsciiString {
    let cstr = unsafe { CStr::from_ptr(buf.as_ptr()).to_bytes() };
    let s = unsafe { AsciiString::from_ascii_unchecked(cstr) };
    s
}

fn char_slice_to_ascii_slice(buf: &[::std::os::raw::c_char]) -> &AsciiStr {
    let cstr = unsafe { CStr::from_ptr(buf.as_ptr()).to_bytes() };
    unsafe { AsciiStr::from_ascii_unchecked(cstr) }
}


#[doc = "Universal handler of different VMD molfile file formats"]
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

        // We can't open file here because for writing we need to know natoms,
        // which is only visible in actuall call to write.
        // For reading we can open, but for consistency we'll defer it as well.

        Ok(VmdMolFileHandler {
            file_name: fname.to_owned(),
            plugin,
            file_handle: ptr::null_mut(),
            mode: OpenMode::Read, // Defaults to read
            natoms: 0,
        })
    }

    fn open_read(&mut self) -> Result<()> {
        // Open file and get file pointer
        // Prepare c-strings for file opening
        let f_type = CString::new("")?; // Not used
        let f_name = CString::new(self.file_name.clone())?;

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

    fn open_write(&mut self) -> Result<()> {
        // Open file and get file pointer
        // Prepare c-strings for file opening
        let f_type = CString::new("")?; // Not used
        let f_name = CString::new(self.file_name.clone())?;

        self.file_handle = unsafe {
            self.plugin.open_file_write.unwrap() // get function ptr
                (f_name.as_ptr(), f_type.as_ptr(), self.natoms as i32) // Call function
        };

        if self.file_handle == ptr::null_mut() {
            bail!("Can't open file {} for writing!", self.file_name);
        }

        self.mode = OpenMode::Write;

        Ok(())
    }
}

impl IoReader for VmdMolFileHandler<'_> {
    fn new_reader(fname: &str) -> Result<Self> {
        let mut instance = Self::new(fname)?;
        instance.open_read()?;
        Ok(instance)
    }
}

impl IoWriter for VmdMolFileHandler<'_> {
    fn new_writer(fname: &str) -> Result<Self> {
        let mut instance = Self::new(fname)?;
        instance.open_write()?;
        Ok(instance)
    }
}

impl IoStructureReader for VmdMolFileHandler<'_> {
    fn read_structure(&mut self) -> Result<Structure> {
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
        let mut structure: Structure = Default::default();
        structure.atoms.reserve(self.natoms);

        for ref at in vmd_atoms {
            let mut new_atom = Atom {
                name: char_slice_to_ascii_str(&at.name),
                resid: at.resid,
                resname: char_slice_to_ascii_str(&at.resname),
                chain: char_slice_to_ascii_str(&at.chain).first().unwrap(),
                charge: at.charge,
                occupancy: at.occupancy,
                bfactor: at.bfactor,
                ..Default::default()
            };
            // See if the element number and mass are set
            // and guess if required
            if at.atomicnumber == 0 || at.mass==0.0 {
                new_atom.guess_element_from_name();
            } else {
                new_atom.atomic_number = at.atomicnumber as usize;
                new_atom.mass = at.mass;
            }
            structure.atoms.push(new_atom);
        }

        Ok(structure)
    }
}

impl IoStructureWriter for VmdMolFileHandler<'_> {
    fn write_structure(&mut self, data: &Structure) -> Result<()> {
        let vmd_atoms = Vec::<molfile_atom_t>::with_capacity(data.atoms.len());
        for ref vmd_at in zip(vmd_atoms) {
            let a = char_slice_to_ascii_slice(&vmd_at.name);
        }

        for(int i=0; i<sel.size(); ++i){
            strcpy( atoms[i].name, sel.name(i).c_str() );
            strcpy( atoms[i].resname, sel.resname(i).c_str() );
            atoms[i].resid = sel.resid(i);
            //stringstream ss;
            //ss << sel.chain(i);
            //strcpy( atoms[i].chain, ss.str().c_str() );
            atoms[i].chain[0] = sel.chain(i);
            atoms[i].chain[1] = '\0';
            atoms[i].occupancy = sel.occupancy(i);
            atoms[i].bfactor = sel.beta(i);
            atoms[i].mass = sel.mass(i);
            atoms[i].charge = sel.charge(i);
            atoms[i].atomicnumber = sel.atomic_number(i);

            // For MOL2 we also need to set atom type as a string
            strcpy( atoms[i].type, sel.element_name(i).c_str() );
        }
        int flags = MOLFILE_OCCUPANCY | MOLFILE_BFACTOR | MOLFILE_ATOMICNUMBER
                    | MOLFILE_CHARGE | MOLFILE_MASS;
        plugin->write_structure(handle,flags,&atoms.front());
        Ok(())
    }
}

impl IoStateReader for VmdMolFileHandler<'_> {
    fn read_next_state(&mut self) -> Result<Option<State>> {
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
            state.box_ = box_from_vmd(ts.A, ts.B, ts.C, ts.alpha, ts.beta, ts.gamma);
            // time
            state.time = ts.physical_time as f32;
        }

        match ret {
            MOLFILE_SUCCESS => Ok(Some(state)),
            MOLFILE_EOF => Ok(None),
            _ => bail!("Error reading timestep!"),
        }
    }
}


impl IoStateWriter for VmdMolFileHandler<'_> {
    fn write_next_state(&mut self, data: &State) -> Result<()> {
        
        Ok(())
    }
}


impl Drop for VmdMolFileHandler<'_> {
    fn drop(&mut self) {
        if self.file_handle != ptr::null_mut() {
            match self.mode {
                OpenMode::Read => unsafe {
                    self.plugin.close_file_read.unwrap()(self.file_handle)
                },
                OpenMode::Write => unsafe {
                    self.plugin.close_file_write.unwrap()(self.file_handle)
                },
            };
        }
    }
}
