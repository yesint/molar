fn main() {
    // Check if Gromacs env variables are provided. If not disable gromacs support
    let src_env = option_env!("GROMACS_SOURCE_DIR");
    let bld_env = option_env!("GROMACS_BUILD_DIR");
    let lib_env = option_env!("GROMACS_LIB_DIR");
    // Tell cargo to check for 'gromacs' config option
    println!("cargo::rustc-check-cfg=cfg(gromacs)");
    // We don't check if the pathes are correct. If not molar_gromacs will just fail to compile.
    if let (Some(_), Some(_), Some(_)) = (src_env, bld_env, lib_env) {
        println!("cargo:rustc-cfg=gromacs");
        println!("cargo::warning=\"Gromacs support enabled\"");
    } else {
        println!("cargo::warning=\"Gromacs support disabled\"");
    }
}