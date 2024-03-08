use cmake::Config;

fn main() {
    let dst = Config::new("powersasa").profile("Release").build();
    println!("cargo:rustc-link-search=native={}", dst.display());
    println!("cargo:rustc-link-lib=static=powersasa_wrapper");
}