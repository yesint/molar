#[cfg(feature = "netcdf")]
use molar::prelude::*;

#[cfg(feature = "netcdf")]
const NC_FILE: &str = "tests/benzene.nc";
#[cfg(feature = "netcdf")]
const XTC_FILE: &str = "tests/benzene.xtc";

/// Read all frames from a trajectory file and return them as a Vec<State>.
#[cfg(feature = "netcdf")]
fn read_all_frames(path: &str) -> Vec<State> {
    FileHandler::open(path)
        .expect("failed to open file")
        .into_iter()
        .collect()
}

/// Verify that coordinates from the NetCDF trajectory match those from the
/// reference XTC trajectory within the expected numeric tolerance.
///
/// XTC uses lossy float32 compression (default precision 1000 pos/nm),
/// so a tolerance of 1e-3 nm is appropriate.
#[test]
#[cfg(feature = "netcdf")]
fn netcdf_matches_xtc() {
    let nc_frames = read_all_frames(NC_FILE);
    let xtc_frames = read_all_frames(XTC_FILE);

    assert_eq!(
        nc_frames.len(),
        xtc_frames.len(),
        "frame count mismatch: nc={} xtc={}",
        nc_frames.len(),
        xtc_frames.len()
    );

    for (i, (nc, xtc)) in nc_frames.iter().zip(xtc_frames.iter()).enumerate() {
        assert_eq!(
            nc.coords.len(),
            xtc.coords.len(),
            "frame {i}: atom count mismatch"
        );

        assert!(
            (nc.time - xtc.time).abs() < 0.01,
            "frame {i}: time mismatch: nc={} xtc={}",
            nc.time,
            xtc.time
        );

        for (j, (p_nc, p_xtc)) in nc.coords.iter().zip(xtc.coords.iter()).enumerate() {
            let dist = (p_nc - p_xtc).norm();
            assert!(
                dist < 1e-3,
                "frame {i}, atom {j}: coordinate mismatch {dist:.6} nm > 1e-3 nm\n  nc={p_nc}\n  xtc={p_xtc}"
            );
        }
    }
}
