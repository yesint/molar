
#[allow(dead_code)]
#[allow(non_upper_case_globals)]
#[allow(non_camel_case_types)]
#[allow(non_snake_case)]
pub mod powersasa_bindings;

// Helper for wrapping closure into C callback with context.
// From here: 
// https://users.rust-lang.org/t/passing-a-closure-to-an-external-c-ffi-library/100271/2

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

pub fn sasa(
    natoms: usize,
    probe_r: f32,
    mut pos_callback: impl FnMut(usize)-> *mut f32, 
    mut vdw_callback: impl FnMut(usize)-> f32,
) -> (Vec<f32>,Vec<f32>) {
    // Arrays for areas and volumes
    let mut areas = Vec::<f32>::with_capacity(natoms);
    let mut volumes = Vec::<f32>::with_capacity(natoms);       
    
    let crd_cb = CCallback::new(&mut pos_callback);

    let vdw_closure = &mut |i: usize| { probe_r + vdw_callback(i) };
    let vdw_cb = CCallback::new(vdw_closure);

    unsafe {
        powersasa_bindings::run_powersasa(
            areas.as_mut_ptr(),
            volumes.as_mut_ptr(),
            Some(crd_cb.function),
            Some(vdw_cb.function),
            crd_cb.user_data,
            vdw_cb.user_data,
            natoms
        );
        // Resize vectors accordingly
        areas.set_len(natoms);
        volumes.set_len(natoms);
    }

    (areas,volumes)
}