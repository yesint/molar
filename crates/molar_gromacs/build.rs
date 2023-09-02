use cmake::Config;

fn main() { 
    // Compilation of Gromacs. SLOW!
    let dst = Config::new("gromacs").profile("Release").build();

    println!("cargo:rustc-link-search=native={}", dst.display());
    println!("cargo:rustc-link-lib=stdc++");
    println!("cargo:rustc-link-lib=static=gromacs");
    println!("cargo:rustc-link-lib=static=muparser");
    println!("cargo:rustc-link-lib=static=gromacs_wrapper");
}