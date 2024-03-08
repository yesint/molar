
#[allow(dead_code)]
#[allow(non_upper_case_globals)]
#[allow(non_camel_case_types)]
#[allow(non_snake_case)]
pub mod powersasa_bindings {
    include!(concat!(env!("OUT_DIR"), "/powersasa_bindings.rs"));
}

// Helper for wrapping closure into C callback with context.
// From here: 
// https://users.rust-lang.org/t/passing-a-closure-to-an-external-c-ffi-library/100271/2

use std::os::raw::c_void;
use std::marker::PhantomData;

pub struct CCallback<'closure, Arg, Ret> {
    pub function: unsafe extern "C" fn(*mut c_void, Arg) -> Ret,
    pub user_data: *mut c_void,

    // without this it's too easy to accidentally drop the closure too soon
    _lifetime: PhantomData<&'closure mut c_void>,
}

impl<'closure, Arg, Ret> CCallback<'closure, Arg, Ret> {
    pub fn new<F>(closure: &'closure mut F) -> Self where F: FnMut(Arg) -> Ret {
        let function: unsafe extern "C" fn(*mut c_void, Arg) -> Ret = Self::call_closure::<F>;

        debug_assert_eq!(std::mem::size_of::<&'closure mut F>(), std::mem::size_of::<*const c_void>());
        debug_assert_eq!(std::mem::size_of_val(&function), std::mem::size_of::<*const c_void>());

        Self {
            function,
            user_data: closure as *mut F as *mut c_void,
            _lifetime: PhantomData,
        }
    }

    unsafe extern "C" fn call_closure<F>(user_data: *mut c_void, arg: Arg) -> Ret where F: FnMut(Arg) -> Ret {
        let cb: &mut F = user_data.cast::<F>().as_mut().unwrap();
        (*cb)(arg)
    }
}

/*
fn main() {
    let mut v = Vec::new();

    // must assign to a variable to ensure it lives until the end of scope
    let closure = &mut |x: i32| { v.push(x) };
    let c = CCallback::new(closure);

    unsafe { (c.function)(c.user_data, 123) };
    
    assert_eq!(v, [123]);
}
*/