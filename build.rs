fn main() {
    println!("cargo:rustc-link-lib=static=gromacs");
    println!("cargo:rustc-link-lib=static=muparser");
    println!("cargo:rustc-link-lib=static=gromacs_wrapper");
}