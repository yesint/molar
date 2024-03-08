use std::{env, path::PathBuf};

use cmake::Config;

fn main() {
    let dst = Config::new("powersasa").profile("Release").build();
    println!("cargo:rustc-link-search=native={}", dst.display());
    //println!("cargo:rustc-link-lib=stdc++");
    println!("cargo:rustc-link-lib=static=powersasa_wrapper");
    
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    
    let bindings = bindgen::Builder::default()
        .header("powersasa/wrapper.hpp")
        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))        
        .allowlist_function("powersasa")
        .clang_args(["-x","c++"])
        .clang_arg("-std=c++20")        
        .layout_tests(false)
        // Finish the builder and generate the bindings.
        .generate()
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    bindings
        .write_to_file(out_path.join("powersasa_bindings.rs"))
        .expect("Couldn't write bindings!"); 
}