use cmake::Config;

fn main() {
    // Run CMake for dependencies
    for dep_name in ["molfile","xdrfile"] {
        let dst = Config::new(format!("src/{}",dep_name))
            .define("CMAKE_BUILD_TYPE","Release")
            .build();
        println!("cargo:rustc-link-search=native={}", dst.display());
        println!("cargo:rustc-link-lib=static={}", dep_name);
    }
}