use std::ffi::CString;
use std::os::raw::{c_char,c_void,c_int};
use std::ptr;
use molfile_bindings::{molfile_plugin_t, pdb_get_plugin_ptr, molfile_atom_t, MOLFILE_SUCCESS, MOLFILE_EOF, molfile_timestep_t};
use std::mem::MaybeUninit;

#[allow(dead_code)]
#[allow(non_upper_case_globals)]
#[allow(non_camel_case_types)]
#[allow(non_snake_case)]
mod molfile_bindings;

enum OpenMode {Read, Write}

pub struct PdbFile<'a> {
    // Plugin pointer
    plugin: &'a molfile_plugin_t,
    // File handle
    //file_handle: *mut c_void,
    file_handle: *mut c_void,
    // Number of atoms
    natoms: i32,
    mode: OpenMode,
}

impl PdbFile<'_> {
    pub fn new() -> Self {
        // Get plugin pointer
        // C funtion registers plugin on first call 
        // and returns stored pointer on later invocations
        let plugin = unsafe { pdb_get_plugin_ptr().as_ref().unwrap() };
        println!("Rust: {:?}",unsafe{pdb_get_plugin_ptr()});
        PdbFile {plugin, file_handle: ptr::null_mut(), natoms: 0, mode: OpenMode::Read}
    }

    
    pub fn open_read(&mut self, fname: &str) -> &Self {
        // Prepare c-strings for file opening
        let f_type = CString::new("pdb").unwrap();
        let f_name = CString::new(fname).unwrap();
        
        self.file_handle = unsafe {
            self.plugin.open_file_read.unwrap() // get function ptr
            (f_name.as_ptr(), f_type.as_ptr(), &mut self.natoms) // call function
        };

        if self.file_handle == ptr::null_mut() {
            panic!("Can't open file {fname} for reading!");
        }
        
        self.mode = OpenMode::Read;

        self
    }

    pub fn read(&self) -> &Self {
        let mut optflags: i32 = 0;
        let el: molfile_atom_t = Default::default();
        let mut atoms = vec![el; self.natoms as usize];

        let ret = unsafe{
            self.plugin.read_structure.unwrap()
            (self.file_handle, &mut optflags, atoms.as_mut_ptr()) 
        };

        if ret != MOLFILE_SUCCESS {
            panic!("Error reading PDB structure!");
        }

        let mut coord = vec![0.0; 3*self.natoms as usize];
        let mut ts = molfile_timestep_t{
            coords: coord.as_mut_ptr(),
            velocities: ptr::null_mut(),
            A: 0.0, B: 0.0, C: 0.0, alpha: 0.0, beta: 0.0, gamma: 0.0,
            physical_time: 0.0,
        };

        let ret = unsafe{
            self.plugin.read_next_timestep.unwrap()
            (self.file_handle,self.natoms,&mut ts) 
        };

        match ret {
            MOLFILE_SUCCESS => (), // pass trough
            MOLFILE_EOF => return self, // return without reading
            _ => panic!("Error reading PDB timestep!"),
        }

        // Convert to pteros data structures

        self
    }
    

}

impl Drop for PdbFile<'_> {
    fn drop(&mut self) {
        if self.file_handle != ptr::null_mut() {
            match self.mode {
                OpenMode::Read => unsafe{ self.plugin.close_file_read.unwrap()(self.file_handle) },
                OpenMode::Write => unsafe{ self.plugin.close_file_write.unwrap()(self.file_handle) },
            };
        }
    }
}
