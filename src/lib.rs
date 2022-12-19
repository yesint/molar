use std::ffi::{CString,CStr};
use std::os::raw::{c_char,c_void,c_int};
use std::{ptr, error};
use molfile_bindings::{molfile_plugin_t, molfile_atom_t,
    pdb_get_plugin_ptr, xyz_get_plugin_ptr, dcd_get_plugin_ptr,
    MOLFILE_SUCCESS, MOLFILE_EOF, molfile_timestep_t};
use std::io::{Error,ErrorKind};
use std::path::Path;
use std::default::Default;
use ascii::{AsciiString,AsciiChar};

#[allow(dead_code)]
#[allow(non_upper_case_globals)]
#[allow(non_camel_case_types)]
#[allow(non_snake_case)]
mod molfile_bindings;


#[derive(Debug,Default,Clone)]
pub struct Atom {
    name: AsciiString,
    resname: AsciiString,
    resid: isize,
    //resindex: usize,
    mass: f32,
    charge: f32,
    chain: AsciiChar,
}

impl Atom {
    pub fn new() -> Self {
        Default::default()
    }
}


#[derive(Debug,Default,Clone)]
pub struct Structure {
    pub atoms: Vec<Atom>,
    bonds: Vec<[usize;2]>,
    angles: Vec<[usize;3]>,
    dihedrals: Vec<[usize;4]>,
}

impl Structure {
    pub fn new() -> Self {
        Default::default()
    }
}


#[derive(Debug,Default)]
pub struct State {
    pub coords: Vec<[f32;3]>,
}

impl State {
    pub fn new() -> Self {
        Default::default()
    }
}


// Traits for different file types
trait MolfileStructure {
    fn read_structure(&self) -> Result<Structure,std::io::Error>;
    fn write_structure(&self, data: &Structure) -> Result<(),std::io::Error>;
}

trait MolfileSingleFrame {
    fn read_state(&self) -> Result<State,std::io::Error>;
    fn write_state(&self, data: &State) -> Result<(),std::io::Error>;
}


enum OpenMode {Read, Write, None}

pub struct MolFileHandler<'a> {
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
fn c_buf_to_ascii_str(buf: &[::std::os::raw::c_char]) -> Result<AsciiString,Box<dyn error::Error>> {
    let cstr = unsafe { CStr::from_ptr(buf.as_ptr()).to_bytes() };
    let s = unsafe { AsciiString::from_ascii_unchecked(cstr) };
    Ok(s)
}

#[doc = "Universal handler of different file formats"]
impl MolFileHandler<'_> {
    
    pub fn new(fname: &str) -> Self {
        // Get extention
        let ext = Path::new(fname)
            .extension().expect("File with extension expected!")
            .to_str().unwrap();
        
        // Get plugin pointer
        // C funtion registers plugin on first call 
        // and returns stored pointer on later invocations
        let plugin = unsafe {
            match ext {
                "pdb" => pdb_get_plugin_ptr(),
                "xyz" => xyz_get_plugin_ptr(),
                "dcd" => dcd_get_plugin_ptr(),
                &_ => panic!("Unrecognized extention {ext}!"),
            }.as_ref().unwrap()
        };

        // We can't open file here because for writing we need to know natoms,
        // which is only visible in actuall call to write.
        // For reading we can open, but for consistency we'll defer it as well.
        
        // Return file handler
        MolFileHandler {
            file_name: fname.to_owned(),
            plugin,
            file_handle: ptr::null_mut(),
            mode: OpenMode::None, // Not opened yet
            natoms: 0
        }
    }
    
    // Returns number of atoms
    fn open_read(&mut self) -> Result<(),Box<dyn error::Error>> {
        // Open file and get file pointer
        // Prepare c-strings for file opening
        let f_type = CString::new("").unwrap(); // Not used
        let f_name = CString::new(self.file_name.clone()).unwrap();

        let mut n: i32 = 0;
        self.file_handle = unsafe {
            self.plugin.open_file_read.unwrap() // get function ptr
                (f_name.as_ptr(), f_type.as_ptr(), &mut n) // Call function
        };
        self.natoms = n as usize;

        if self.file_handle == ptr::null_mut() {
            Err(format!("Can't open file {} for reading!",self.file_name))?;
        }

        self.mode = OpenMode::Read;
        
        Ok(())
    }

    fn open_write(&mut self, natoms: i32) {
        // Open file and get file pointer
        // Prepare c-strings for file opening
        let f_type = CString::new("").unwrap(); // Not used
        let f_name = CString::new(self.file_name.clone()).unwrap();

        self.file_handle = unsafe {
            self.plugin.open_file_write.unwrap() // get function ptr
                (f_name.as_ptr(), f_type.as_ptr(), natoms) // Call function
        };

        if self.file_handle == ptr::null_mut() {
            panic!("Can't open file {} for writing!",self.file_name);
        }

        self.mode = OpenMode::Write;
    }
    
    pub fn read_structure(&mut self) -> Result<Structure,Box<dyn error::Error>> {
        // Open file for reading
        match self.mode {
            OpenMode::None => self.open_read(), // Open file
            OpenMode::Read => Ok(()),
            OpenMode::Write => Err("File is already opened for writing")?,
        };

        let mut optflags: i32 = 0;
        let el: molfile_atom_t = Default::default();
        let mut atoms = vec![el; self.natoms];

        let ret = unsafe{
            self.plugin.read_structure.unwrap()
            (self.file_handle, &mut optflags, atoms.as_mut_ptr()) 
        };

        if ret != MOLFILE_SUCCESS {
            Err("Plugin error reading structure!")?;
        }

        // Convert to Structure
        let mut structure: Structure = Default::default();
        structure.atoms.reserve(self.natoms);
        for ref at in atoms {
            let atom = Atom {
                name: c_buf_to_ascii_str(&at.name)?,
                resid: at.resid as isize,
                resname: c_buf_to_ascii_str(&at.resname)?,
                chain: c_buf_to_ascii_str(&at.chain)?.first().unwrap(),
                mass: at.mass,
                charge: at.charge,
            };
            structure.atoms.push(atom);
        }

        Ok(structure)
    }
    
    pub fn read_state(&mut self) -> Result<State,Box<dyn error::Error>> {
        // Open file for reading
        match self.mode {
            OpenMode::None => Err("Can't read state before structure")?,
            OpenMode::Read => (),
            OpenMode::Write => Err("File is already opened for writing")?,
        }

        let mut st: State = Default::default();
        st.coords = vec![[0.0,0.0,0.0]; self.natoms as usize]; 

        //let mut coord = vec![0.0; 3*self.natoms as usize];
        let mut ts = molfile_timestep_t{
            coords: st.coords.as_mut_ptr().cast::<f32>(),
            velocities: ptr::null_mut(),
            A: 0.0, B: 0.0, C: 0.0, alpha: 0.0, beta: 0.0, gamma: 0.0,
            physical_time: 0.0,
        };

        let ret = unsafe{
            self.plugin.read_next_timestep.unwrap()
            (self.file_handle,self.natoms.try_into().unwrap(),&mut ts) 
        };

        match ret {
            MOLFILE_SUCCESS => (), // pass trough
            MOLFILE_EOF => (), // return without reading
            _ => Err("Error reading PDB timestep!")?,
        }

        // Convert to pteros data structures
         

        Ok(st)
    }

}

impl Drop for MolFileHandler<'_> {
    fn drop(&mut self) {
        if self.file_handle != ptr::null_mut() {
            match self.mode {
                OpenMode::Read => unsafe{ self.plugin.close_file_read.unwrap()(self.file_handle) },
                OpenMode::Write => unsafe{ self.plugin.close_file_write.unwrap()(self.file_handle) },
                OpenMode::None => panic!("This should never happen!"),
            };
        }
    }
}
