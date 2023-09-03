use cmake::Config;
use std::env;
use std::path::PathBuf;
use std::fs::read_to_string;

fn main() { 
    // Compilation of Gromacs. SLOW!
    let dst = Config::new("gromacs").profile("Release").build();

    println!("cargo:rustc-link-search=native={}", dst.display());
    println!("cargo:rustc-link-lib=stdc++");
    println!("cargo:rustc-link-lib=static=gromacs");
    println!("cargo:rustc-link-lib=static=muparser");
    println!("cargo:rustc-link-lib=static=gromacs_wrapper");

    // Bindgen
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());

    // Read text file wih pathes written by CMake
    let content = read_to_string(out_path.join("cmake_to_cargo.txt")).unwrap();
    let mut it = content.lines();
    let gmx_source_dir = it.next().unwrap();
    let gmx_binary_dir = it.next().unwrap();

    let bindings = bindgen::Builder::default()
        .header("gromacs/wrapper.hpp")
        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
        .allowlist_type("t_topology")
        .allowlist_type("TprHelper")
        //.allowlist_function("drop_top_mtop")
        .opaque_type("std::.*")
        .opaque_type("t_state")
        .opaque_type("t_inputrec")
        .opaque_type("gmx_mtop_t")
        .clang_args(["-x","c++"])
        .clang_arg("-std=c++17")
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
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    bindings
        .write_to_file(out_path.join("gromacs_bindings.rs"))
        .expect("Couldn't write bindings!");
}