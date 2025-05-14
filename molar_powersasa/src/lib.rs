
#[allow(dead_code)]
#[allow(non_upper_case_globals)]
#[allow(non_camel_case_types)]
#[allow(non_snake_case)]
pub mod powersasa_bindings;

// Helper for wrapping closure into C callback with context.
// From here: 
// https://users.rust-lang.org/t/passing-a-closure-to-an-external-c-ffi-library/100271/2

use std::cell::OnceCell;
use std::os::raw::c_void;
use std::marker::PhantomData;

struct CCallback<'closure, Arg, Ret> {
    pub function: unsafe extern "C" fn(*mut c_void, Arg) -> Ret,
    pub user_data: *mut c_void,

    // without this it's too easy to accidentally drop the closure too soon
    _lifetime: PhantomData<&'closure mut c_void>,
}

impl<'closure, Arg, Ret> CCallback<'closure, Arg, Ret> {
    pub fn new<F>(closure: &'closure mut F) -> Self 
    where 
        F: FnMut(Arg) -> Ret 
    {
        let function: unsafe extern "C" fn(*mut c_void, Arg) -> Ret = Self::call_closure::<F>;

        debug_assert_eq!(std::mem::size_of::<&'closure mut F>(), std::mem::size_of::<*const c_void>());
        debug_assert_eq!(std::mem::size_of_val(&function), std::mem::size_of::<*const c_void>());

        Self {
            function,
            user_data: closure as *mut F as *mut c_void,
            _lifetime: PhantomData,
        }
    }

    unsafe extern "C" fn call_closure<F>(user_data: *mut c_void, arg: Arg) -> Ret 
    where 
        F: FnMut(Arg) -> Ret
    {
        let cb: &mut F = user_data.cast::<F>().as_mut().unwrap();
        (*cb)(arg)
    }
}

/// Result of SASA calculation containing areas and volumes (total and per atom)
#[derive(Debug)]
pub struct SasaResults {
    areas: Vec<f32>,
    volumes: Vec<f32>,
    total_area: OnceCell<f32>,
    total_volume: OnceCell<f32>,
}

impl SasaResults {
    fn new(natoms: usize) -> Self {
        Self {
            areas: Vec::with_capacity(natoms),
            volumes: Vec::with_capacity(natoms),
            total_area: Default::default(),
            total_volume: Default::default(),
        }
    }

    pub fn areas(&self) -> &[f32] {
        &self.areas
    }

    pub fn volumes(&self) -> &[f32] {
        &self.volumes
    }

    pub fn total_area(&self) -> f32 {
        *self.total_area.get_or_init(|| self.areas.iter().sum())
    }

    pub fn total_volume(&self) -> f32 {
        *self.total_volume.get_or_init(|| self.volumes.iter().sum())
    }
}

pub fn compute_sasa(
    natoms: usize,
    probe_r: f32,
    mut pos_callback: impl FnMut(usize)-> *mut f32, 
    mut vdw_callback: impl FnMut(usize)-> f32,
) -> SasaResults {
    // results
    let mut res = SasaResults::new(natoms);
    
    let crd_cb = CCallback::new(&mut pos_callback);

    let vdw_closure = &mut |i: usize| { probe_r + vdw_callback(i) };
    let vdw_cb = CCallback::new(vdw_closure);

    unsafe {
        powersasa_bindings::run_powersasa(
            res.areas.as_mut_ptr(),
            res.volumes.as_mut_ptr(),
            Some(crd_cb.function),
            Some(vdw_cb.function),
            crd_cb.user_data,
            vdw_cb.user_data,
            natoms
        );
        // Resize vectors accordingly
        res.areas.set_len(natoms);
        res.volumes.set_len(natoms);
    }

    res
}