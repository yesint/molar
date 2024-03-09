use std::{env, path::PathBuf};

fn main() {
    cxx_build::bridge("src/main.rs") // returns a cc::Build
        .file("powersasa/wrapper.cpp")
        .std("c++17")
        .compile("cxxbridge");

    println!("cargo:rerun-if-changed=src/main.rs");
    println!("cargo:rerun-if-changed=powersasa/wrapper.cpp");
    println!("cargo:rerun-if-changed=powersasa/wrapper.h");
}
