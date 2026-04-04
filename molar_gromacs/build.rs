use std::path::PathBuf;

fn main() {
    let src_env = option_env!("GROMACS_SOURCE_DIR");
    let bld_env = option_env!("GROMACS_BUILD_DIR");
    let lib_env = option_env!("GROMACS_LIB_DIR");

    if let (Some(src_env), Some(bld_env), Some(lib_env)) = (src_env, bld_env, lib_env) {
        let out_dir = PathBuf::from(std::env::var("OUT_DIR").unwrap());

        // Determine shared library filename for the current target OS.
        let so_name = libloading::library_filename("molar_gromacs_plugin");
        let so_path = out_dir.join(&so_name);

        // Step 1: compile wrapper.cpp to an object file via the cc crate.
        // This lets cc handle compiler detection, CXX env var, PIC, etc.
        let mut build = cc::Build::new();
        build
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
            .include(format!("{bld_env}/src/include"))
            .pic(true)
            .warnings(false);

        // thread_mpi headers are only needed when Gromacs is built with
        // GMX_THREAD_MPI=1 (indicated in config.h).  Add the path only if
        // the directory is present so builds without thread-MPI are unaffected.
        let tmpi_inc = format!("{src_env}/src/external/thread_mpi/include");
        if std::path::Path::new(&tmpi_inc).is_dir() {
            build.include(&tmpi_inc);
        }

        let obj_files = build.compile_intermediates();

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

        assert!(status.success(), "failed to build {}", so_path.display());

        // Bake the plugin path as the default value of MOLAR_GROMACS_PLUGIN.
        // Users can override it at runtime by setting the same env var.
        println!("cargo::rustc-env=MOLAR_GROMACS_PLUGIN={}", so_path.display());
        println!("cargo::warning=Gromacs plugin built: {}", so_path.display());
    } else {
        println!("cargo::warning=Gromacs plugin in NOT built");
    }
}
