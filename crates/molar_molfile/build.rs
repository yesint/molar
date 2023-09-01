use cmake::Config;

fn main() {    
    let dst = Config::new(format!("molfile")).profile("Release").build();
    println!("cargo:rustc-link-search=native={}", dst.display());
    println!("cargo:rustc-link-lib=static=molfile");
    println!("cargo:rerun-if-changed=build.rs");
}