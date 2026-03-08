use std::path::PathBuf;

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

        // Step 1: compile wrapper.cpp to an object file via the cc crate.
        // This lets cc handle compiler detection, CXX env var, PIC, etc.
        let obj_files = cc::Build::new()
            .cpp(true)
            .file("gromacs/wrapper.cpp")
            .include("gromacs")
            .include(format!("{src_env}/src"))
            .include(format!("{src_env}/src/gromacs/utility/include"))
            .include(format!("{src_env}/src/gromacs/math/include"))
            .include(format!("{src_env}/src/gromacs/topology/include"))
            .include(format!("{src_env}/api/legacy/include"))
            .include(format!("{src_env}/src/external"))
            .include(format!("{bld_env}/api/legacy/include"))
            .pic(true)
            .warnings(false)
            .compile_intermediates();

        // Step 2: link object file(s) into a shared library using the same compiler.
        let compiler = cc::Build::new().cpp(true).get_compiler();
        let status = compiler
            .to_command()
            .arg("-shared")
            .args(&obj_files)
            .arg(format!("-L{lib_env}"))
            .arg("-lgromacs")
            .arg("-lmuparser")
            .arg("-o")
            .arg(&so_path)
            .status()
            .expect("failed to link shared library");

        assert!(status.success(), "failed to build {so_name}");

        // Bake the plugin path as the default value of MOLAR_GROMACS_PLUGIN.
        // Users can override it at runtime by setting the same env var.
        println!("cargo::rustc-env=MOLAR_GROMACS_PLUGIN={}", so_path.display());
        println!("cargo::warning=Gromacs plugin built: {}", so_path.display());
    }
}
