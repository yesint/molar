
#[allow(dead_code)]
#[allow(non_upper_case_globals)]
#[allow(non_camel_case_types)]
#[allow(non_snake_case)]
pub mod gromacs_bindings {
    #[cfg(feature = "gen_bindings")]
    include!(concat!(env!("OUT_DIR"), "/gromacs_bindings.rs"));
    #[cfg(not(feature = "gen_bindings"))]
    include!("gromacs_bindings.rs");
}



