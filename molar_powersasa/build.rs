use flate2::read::GzDecoder;
use std::fs::File;
use std::path::PathBuf;
use tar::Archive;

fn main() {
    let eigen_path = if let Ok(local_eigen) = std::env::var("EIGEN_DIR") {
        PathBuf::from(local_eigen)
    } else {
        // Unpack bundled eigen if not already unpacked
        let dst = PathBuf::from(std::env::var("OUT_DIR").unwrap()).join("eigen/");
        let include_path = dst.join("eigen-3.4.0");
        if !include_path.exists() {
            println!("cargo::warning=unpacking bundled Eigen headers...");
            let tar_gz = File::open("eigen-3.4.0.tar.gz").unwrap();
            let tar = GzDecoder::new(tar_gz);
            let mut archive = Archive::new(tar);
            archive.unpack(dst).unwrap();
        }
        include_path
    };
    println!("cargo::warning=Eigen headers in {:?}", eigen_path);

    cc::Build::new()
        .cpp(true)
        .file("powersasa/wrapper.cpp")
        .include(eigen_path)
        .pic(true)
        .flag_if_supported("-std=c++17")
        .warnings(false)
        .compile("powersasa");
}
