use std::{env, path::PathBuf};

use cmake::Config;

fn main() {
    //println!("cargo:rerun-if-changed=build.rs");
    //println!("cargo:rerun-if-changed=powersasa/wrapper.cpp");
    //println!("cargo:rerun-if-changed=powersasa/wrapper.hpp");

    let dst = Config::new("powersasa").profile("Release").build();
    println!("cargo:rustc-link-search=native={}", dst.display());
    println!("cargo:rustc-link-lib=static=powersasa_wrapper");
    
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    
    // We are generating C bindings, so no C++ options needed
    let bindings = bindgen::Builder::default()
        .header("powersasa/wrapper.h")
        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))        
        .allowlist_function("run_powersasa")
        .layout_tests(false)
        // Finish the builder and generate the bindings.
        .generate()
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    bindings
        .write_to_file(out_path.join("powersasa_bindings.rs"))
        .expect("Couldn't write bindings!"); 
} 