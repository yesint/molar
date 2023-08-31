use cmake::Config;

fn main() {
    // Run CMake for dependencies
    for dep_name in ["molfile","xdrfile"] {
        let dst = Config::new(format!("src/third_party/{dep_name}"))
            //.define("CMAKE_BUILD_TYPE","Release")
            .profile("Release")
            .build();
        println!("cargo:rustc-link-search=native={}", dst.display());
        println!("cargo:rustc-link-lib=static={}", dep_name);
    }

    let dst = Config::new("src/third_party/gromacs").build();
    println!("cargo:rustc-link-search=native={}", dst.display());
    println!("cargo:rustc-link-lib=static={}", "gromacs");
    println!("cargo:rustc-link-lib=static={}", "muparser");
}