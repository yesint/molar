use std::env;

fn main() {
    let src_env = env!("GROMACS_SOURCE_DIR");
    let bin_env = env!("GROMACS_BINARY_DIR");
    let lib_env = env!("GROMACS_LIB_DIR");

    cc::Build::new()
        .cpp(true)
        .file("gromacs/wrapper.cpp")
        .include("gromacs")
        .include(format!("{src_env}/src"))
        .include(format!("{src_env}/src/gromacs/utility/include"))
        .include(format!("{src_env}/src/gromacs/math/include"))
        .include(format!("{src_env}/src/gromacs/topology/include"))
        .include(format!("{src_env}/api/legacy/include"))
        .include(format!("{src_env}/src/external"))
        .include(format!("{bin_env}/api/legacy/include")) // For generated headers)
        .pic(true)
        .warnings(false)
        .compile("gromacs_wrapper");

    // Link pathes and libs
    println!("cargo:rustc-link-search=native={lib_env}");
    //println!("cargo:rustc-link-lib=stdc++");
    println!("cargo:rustc-link-lib=gromacs");
    println!("cargo:rustc-link-lib=muparser");

    #[cfg(feature = "gen_bindings")]
    {
        // Generate the bindings
        let out_path = std::path::PathBuf::from(env::var("OUT_DIR").unwrap());
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
            .clang_args(["-x", "c++"])
            .clang_arg("-std=c++20")
            .clang_arg(format!("-I{src_env}/src"))
            .clang_arg(format!("-I{src_env}/src/gromacs/utility/include"))
            .clang_arg(format!("-I{src_env}/src/gromacs/math/include"))
            .clang_arg(format!("-I{src_env}/src/gromacs/topology/include"))
            .clang_arg(format!("-I{src_env}/api/legacy/include"))
            .clang_arg(format!("-I{src_env}/src/external"))
            .clang_arg(format!("-I{bin_env}/api/legacy/include"))
            .layout_tests(false)
            // Finish the builder and generate the bindings.
            .generate()
            .expect("able to generate bindings");

        // Write the bindings to the $OUT_DIR/bindings.rs file.
        bindings
            .write_to_file(out_path.join("gromacs_bindings.rs"))
            .expect("should be able to write bindings!");
    }
}
