use cached_path::{Cache, Options};
use std::path::PathBuf;

fn main() {
    // Download Eigen
    let cache = Cache::builder()
        .dir(PathBuf::from(std::env::var("OUT_DIR").unwrap()).join("eigen/"))
        .connect_timeout(std::time::Duration::from_secs(5))
        .build()
        .unwrap();

    let eigen_path = if let Ok(local_eigen) = std::env::var("EIGEN_DIR") {
        PathBuf::from(local_eigen)
    } else {
        cache.cached_path_with_options(
                "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz",
                &Options::default().extract(),
            )
            .unwrap()
            .join("eigen-3.4.0/")
    };
    println!("cargo::warning=Eigen headers in {:?}",eigen_path);

    cc::Build::new()
        .cpp(true)
        .file("powersasa/wrapper.cpp")
        .include(eigen_path)
        .pic(true)
        .flag_if_supported("-std=c++17")
        .warnings(false)
        .compile("powersasa");
}
