use std::env;
use std::path::PathBuf;
use std::fs::read_to_string;

fn main() { 
    // Configure Gromacs
    let mut cfg = cmake::Config::new("gromacs");

    // Transfer env vars to cmake if they are provided to cargo
    let src_env = option_env!("GROMACS_SOURCE_DIR");
    let bin_env = option_env!("GROMACS_BINARY_DIR");
    let lib_env = option_env!("GROMACS_LIB_DIR");

    let external_gmx = src_env.is_some() && bin_env.is_some() && lib_env.is_some();

    #[cfg(not(feature="build_gromacs"))]
    if !external_gmx {
        let msg = "Either set Gromacs ENV variables in .cargo/config.toml or enable build_gromacs feature!";
        println!("cargo:warning={msg}");
        panic!("{msg}")
    }

    let mut gmx_source_dir = String::new();
    let mut gmx_binary_dir = String::new();

    if external_gmx {
        gmx_source_dir = src_env.unwrap().to_owned();
        gmx_binary_dir = bin_env.unwrap().to_owned();
        let gmx_lib_dir = lib_env.unwrap().to_owned();
        
        cfg.configure_arg(format!("-DGROMACS_SOURCE_DIR={gmx_source_dir}"));
        cfg.configure_arg(format!("-DGROMACS_BINARY_DIR={gmx_binary_dir}"));
        cfg.configure_arg(format!("-DGROMACS_LIB_DIR={gmx_lib_dir}"));
        
        println!("cargo:rustc-link-search=native={gmx_lib_dir}");
    }

    // Do CMAKE build (could be very slow if Gromacs is built in place)
    let dst = cfg.profile("Release").build();

    // Link pathes and libs
    println!("cargo:rustc-link-search=native={}", dst.display());
    println!("cargo:rustc-link-lib=stdc++");
    println!("cargo:rustc-link-lib=gromacs");
    println!("cargo:rustc-link-lib=muparser");
    println!("cargo:rustc-link-lib=static=gromacs_wrapper");

    // Generate the bindings
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());

    if !external_gmx {
        // Read text file wih pathes written by CMake
        let content = read_to_string(out_path.join("cmake_to_cargo.txt")).unwrap();
        let mut it = content.lines();
        gmx_source_dir = it.next().unwrap().to_owned();
        gmx_binary_dir = it.next().unwrap().to_owned();
    }
    
    let bindings = bindgen::Builder::default()
        .header("gromacs/wrapper.hpp")
        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .allowlist_type("t_topology")
        .allowlist_type("TprHelper")
        .allowlist_var("F_.*")
        .opaque_type("std::.*")
        .opaque_type("t_state")
        .opaque_type("t_inputrec")
        .opaque_type("gmx_mtop_t")
        .clang_args(["-x","c++"])
        .clang_arg("-std=c++20")
        .clang_arg(format!("-I{gmx_source_dir}/src"))
        .clang_arg(format!("-I{gmx_source_dir}/src/gromacs/utility/include"))
        .clang_arg(format!("-I{gmx_source_dir}/src/gromacs/math/include"))
        .clang_arg(format!("-I{gmx_source_dir}/src/gromacs/topology/include"))
        .clang_arg(format!("-I{gmx_source_dir}/api/legacy/include"))
        .clang_arg(format!("-I{gmx_binary_dir}/api/legacy/include"))
        .clang_arg(format!("-I{gmx_source_dir}/src/external"))
        .layout_tests(false)
        // Finish the builder and generate the bindings.
        .generate()
        .expect("able to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    bindings
        .write_to_file(out_path.join("gromacs_bindings.rs"))
        .expect("able to write bindings!");
}
