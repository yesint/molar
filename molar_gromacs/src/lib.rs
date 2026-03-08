use std::{ffi::CStr, os::raw::c_char};

/// C-ABI atom data, layout matches `TprAtom` in `wrapper.hpp`.
#[repr(C)]
#[derive(Clone, Copy)]
pub struct TprAtom {
    pub name:          [u8; 8],
    pub resname:       [u8; 8],
    pub type_name:     [u8; 8],
    pub chain:         u8,
    // 3 bytes implicit padding (repr(C) follows C alignment rules)
    pub resind:        u32,
    pub type_id:       u32,
    pub atomic_number: u32,
    pub charge:        f32,
    pub mass:          f32,
    pub occupancy:     f32,
    pub bfactor:       f32,
}

/// C-ABI bond pair.
#[repr(C)]
#[derive(Clone, Copy)]
pub struct TprBond {
    pub atom1: u32,
    pub atom2: u32,
}

/// C-ABI molecule range (both ends inclusive).
#[repr(C)]
#[derive(Clone, Copy)]
pub struct TprMolecule {
    pub start: u32,
    pub end:   u32,
}

/// Opaque handle returned by `tpr_open`.
pub enum TprHandle {}

/// Function pointer table resolved from the plugin shared library.
pub struct TprPluginFns {
    pub tpr_open:          unsafe extern "C" fn(*const c_char) -> *mut TprHandle,
    pub tpr_close:         unsafe extern "C" fn(*mut TprHandle),
    pub tpr_last_error:    unsafe extern "C" fn() -> *const c_char,
    pub tpr_natoms:        unsafe extern "C" fn(*mut TprHandle) -> usize,
    pub tpr_nbonds:        unsafe extern "C" fn(*mut TprHandle) -> usize,
    pub tpr_nmolecules:    unsafe extern "C" fn(*mut TprHandle) -> usize,
    pub tpr_fill_atoms:    unsafe extern "C" fn(*mut TprHandle, *mut TprAtom),
    pub tpr_fill_bonds:    unsafe extern "C" fn(*mut TprHandle, *mut TprBond),
    pub tpr_fill_molecules:unsafe extern "C" fn(*mut TprHandle, *mut TprMolecule),
    pub tpr_fill_coords:   unsafe extern "C" fn(*mut TprHandle, *mut f32),
    pub tpr_fill_box:      unsafe extern "C" fn(*mut TprHandle, *mut f32),
}

/// Holds the dynamically loaded plugin library and the resolved function pointers.
///
/// The `Library` must outlive all use of the function pointers.  Wrapping both
/// in a single struct guarantees that invariant.
pub struct TprPlugin {
    _lib: libloading::Library,
    pub fns: TprPluginFns,
}

// TprHandle is heap-allocated C++ data; TprPlugin is just a library + fn ptrs.
unsafe impl Send for TprPlugin {}
unsafe impl Sync for TprPlugin {}

impl TprPlugin {
    /// Try to load the plugin shared library.
    ///
    /// Search order:
    /// 1. `MOLAR_GROMACS_PLUGIN` env var (runtime override, full path).
    /// 2. Path baked in at compile time via `MOLAR_GROMACS_PLUGIN_DEFAULT`
    ///    (set by `molar_gromacs/build.rs` when Gromacs env vars were present).
    /// 3. System library search path (`libmolar_gromacs_plugin.so`).
    pub fn load() -> Result<Self, libloading::Error> {
        let lib = Self::open_library()?;
        let fns = unsafe { Self::resolve(&lib)? };
        Ok(TprPlugin { _lib: lib, fns })
    }

    fn open_library() -> Result<libloading::Library, libloading::Error> {
        // 1. Runtime env var
        if let Ok(path) = std::env::var("MOLAR_GROMACS_PLUGIN") {
            return unsafe { libloading::Library::new(path) };
        }
        // 2. Compile-time default (present only when molar_gromacs was built
        //    with Gromacs env vars set)
        if let Some(path) = option_env!("MOLAR_GROMACS_PLUGIN_DEFAULT") {
            if let Ok(lib) = unsafe { libloading::Library::new(path) } {
                return Ok(lib);
            }
        }
        // 3. System search
        let name = libloading::library_filename("molar_gromacs_plugin");
        unsafe { libloading::Library::new(name) }
    }

    unsafe fn resolve(lib: &libloading::Library) -> Result<TprPluginFns, libloading::Error> {
        macro_rules! sym {
            ($name:literal, $ty:ty) => {
                *lib.get::<$ty>($name)?
            };
        }
        Ok(TprPluginFns {
            tpr_open:           sym!(b"tpr_open\0",           unsafe extern "C" fn(*const c_char) -> *mut TprHandle),
            tpr_close:          sym!(b"tpr_close\0",          unsafe extern "C" fn(*mut TprHandle)),
            tpr_last_error:     sym!(b"tpr_last_error\0",     unsafe extern "C" fn() -> *const c_char),
            tpr_natoms:         sym!(b"tpr_natoms\0",         unsafe extern "C" fn(*mut TprHandle) -> usize),
            tpr_nbonds:         sym!(b"tpr_nbonds\0",         unsafe extern "C" fn(*mut TprHandle) -> usize),
            tpr_nmolecules:     sym!(b"tpr_nmolecules\0",     unsafe extern "C" fn(*mut TprHandle) -> usize),
            tpr_fill_atoms:     sym!(b"tpr_fill_atoms\0",     unsafe extern "C" fn(*mut TprHandle, *mut TprAtom)),
            tpr_fill_bonds:     sym!(b"tpr_fill_bonds\0",     unsafe extern "C" fn(*mut TprHandle, *mut TprBond)),
            tpr_fill_molecules: sym!(b"tpr_fill_molecules\0", unsafe extern "C" fn(*mut TprHandle, *mut TprMolecule)),
            tpr_fill_coords:    sym!(b"tpr_fill_coords\0",    unsafe extern "C" fn(*mut TprHandle, *mut f32)),
            tpr_fill_box:       sym!(b"tpr_fill_box\0",       unsafe extern "C" fn(*mut TprHandle, *mut f32)),
        })
    }

    /// Return the error message from the last failed `tpr_open` call.
    pub fn last_error(&self) -> String {
        unsafe {
            let ptr = (self.fns.tpr_last_error)();
            if ptr.is_null() {
                "unknown error".to_owned()
            } else {
                CStr::from_ptr(ptr).to_string_lossy().into_owned()
            }
        }
    }
}
