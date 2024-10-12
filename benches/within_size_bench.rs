#[cfg(test)]
mod tests {
    use molar::core::Source;
    use std::io::prelude::*;
    use std::{fs::File, io::BufWriter};

    fn run(mut code: impl FnMut()) -> f32 {
        let t = std::time::Instant::now();
        code();
        t.elapsed().as_secs_f32()
    }

    #[test]
    fn within_size_benchmark() -> anyhow::Result<()> {
        let src = Source::serial_from_file("tests/albumin.pdb")?;

        for pbc in vec!["",] {
        //for pbc in vec!["","pbc"] {
            let pbc_str = if pbc != "" {
                "pbc"
            } else {
                ""
            };

            for n_res in vec![1, 20, 40, 60] {
                let mut out =
                    BufWriter::new(File::create(format!("target/molar2_ref_{}{}.dat", pbc, n_res))?);
                for i in 0..40 {
                    let d = 0.3 + 0.1 * i as f32;
                    let t = run(|| {
                        for start_res in (0..100).step_by(10) {
                            src.select_str(format!(
                                "within {} {} of resid {}:{}",
                                d,
                                pbc_str,
                                start_res,
                                start_res + n_res
                            ))
                            .unwrap();
                        }
                    });
                    writeln!(out, "{d} {t}")?;
                    println!("{d} {t}");
                }
            }
        }

        Ok(())
    }
}
