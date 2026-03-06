#[cfg(feature = "netcdf")]
use molar::prelude::*;

#[cfg(feature = "netcdf")]
const NC_FILE: &str = "tests/benzene.nc";
#[cfg(feature = "netcdf")]
const XTC_FILE: &str = "tests/benzene.xtc";
#[cfg(feature = "netcdf")]
const NC_OUT_FILE: &str = "/tmp/benzene_out.nc";

/// Read all frames from a trajectory file and return them as a Vec<State>.
#[cfg(feature = "netcdf")]
fn read_all_frames(path: &str) -> Vec<State> {
    FileHandler::open(path)
        .expect("failed to open file")
        .into_iter()
        .collect()
}

/// Write a slice of states to a NetCDF file and return them.
#[cfg(feature = "netcdf")]
fn write_frames(path: &str, frames: &[State]) {
    let mut fh = FileHandler::create(path).expect("failed to create file");
    for state in frames {
        fh.write_state(state).expect("failed to write frame");
    }
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

    for (i, (nc, xtc)) in nc_frames.iter().zip(&xtc_frames).enumerate() {
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

/// Write the reference NetCDF trajectory to a new file, then re-read it and
/// verify that coordinates, time, and box round-trip exactly (within f32 precision).
#[test]
#[cfg(feature = "netcdf")]
fn netcdf_write_roundtrip() {
    let original = read_all_frames(NC_FILE);

    write_frames(NC_OUT_FILE, &original);

    let roundtrip = read_all_frames(NC_OUT_FILE);

    assert_eq!(original.len(), roundtrip.len(), "frame count changed after write");

    for (i, (orig, rt)) in original.iter().zip(&roundtrip).enumerate() {
        assert_eq!(orig.coords.len(), rt.coords.len(), "frame {i}: atom count changed");

        assert!(
            (orig.time - rt.time).abs() < 1e-4,
            "frame {i}: time mismatch after round-trip: orig={} rt={}",
            orig.time,
            rt.time
        );

        for (j, (p_orig, p_rt)) in orig.coords.iter().zip(rt.coords.iter()).enumerate() {
            let dist = (p_orig - p_rt).norm();
            assert!(
                dist < 1e-4,
                "frame {i}, atom {j}: coordinate mismatch after round-trip {dist:.6} nm\n  orig={p_orig}\n  rt={p_rt}"
            );
        }

        // Verify box round-trip
        match (&orig.pbox, &rt.pbox) {
            (Some(a), Some(b)) => {
                let diff = (a.get_matrix() - b.get_matrix()).abs().max();
                assert!(diff < 1e-4, "frame {i}: box mismatch after round-trip: diff={diff}");
            }
            (None, None) => {}
            _ => panic!("frame {i}: box presence changed after round-trip"),
        }
    }
}
