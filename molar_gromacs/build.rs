use std::{path::PathBuf, process::Command};

fn main() {
    let src_env = option_env!("GROMACS_SOURCE_DIR");
    let bld_env = option_env!("GROMACS_BUILD_DIR");
    let lib_env = option_env!("GROMACS_LIB_DIR");

    if let (Some(src_env), Some(bld_env), Some(lib_env)) = (src_env, bld_env, lib_env) {
        let out_dir = PathBuf::from(std::env::var("OUT_DIR").unwrap());

        // Determine shared library filename for the current target OS.
        let target_os = std::env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
        let so_name = match target_os.as_str() {
            "macos" => "libmolar_gromacs_plugin.dylib",
            _       => "libmolar_gromacs_plugin.so",
        };
        let so_path = out_dir.join(so_name);

        // Use the cc crate only to locate the C++ compiler.
        let compiler = cc::Build::new().cpp(true).get_compiler();

        let status = Command::new(compiler.path())
            .arg("-shared")
            .arg("-fPIC")
            .arg("-std=c++20")
            .arg("-w")                              // suppress warnings
            .arg("gromacs/wrapper.cpp")
            .arg("-Igromacs")
            .arg(format!("-I{src_env}/src"))
            .arg(format!("-I{src_env}/src/gromacs/utility/include"))
            .arg(format!("-I{src_env}/src/gromacs/math/include"))
            .arg(format!("-I{src_env}/src/gromacs/topology/include"))
            .arg(format!("-I{src_env}/api/legacy/include"))
            .arg(format!("-I{src_env}/src/external"))
            .arg(format!("-I{bld_env}/api/legacy/include"))
            .arg(format!("-L{lib_env}"))
            .arg("-lgromacs")
            .arg("-lmuparser")
            .arg("-o")
            .arg(&so_path)
            .status()
            .expect("failed to invoke C++ compiler");

        assert!(status.success(), "failed to build {so_name}");

        println!("cargo::warning=Gromacs plugin built: {}", so_path.display());
        // Bake the path into the binary so the loader can find it by default.
        println!("cargo::rustc-env=MOLAR_GROMACS_PLUGIN_DEFAULT={}", so_path.display());
    }
}
