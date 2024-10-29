use std::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use molar::prelude::*;

fn molar_benchmark(c: &mut Criterion) {
    c.bench_function("align", |b| b.iter(
        black_box(
        || {
            let src = Source::serial_from_file("tests/protein.pdb").unwrap();
            let ref_sel = src.select_all().unwrap();
            let mut cur_sel = src.select_all().unwrap();

            let mut rmsd = vec![];

            let trj = FileHandler::open("tests/protein.xtc").unwrap().into_iter();
            for st in trj {
                cur_sel.set_state(st).unwrap();
                let tr = MeasureMasses::fit_transform(&cur_sel, &ref_sel).unwrap();
                cur_sel.apply_transform(&tr);
                rmsd.push( MeasurePos::rmsd(&cur_sel, &ref_sel).unwrap() );
            }
            //println!("{:?}",&rmsd[..10]);
        }))
    );

    c.bench_function("within", |b| b.iter(
        black_box(|| {
            let src = Source::serial_from_file("tests/protein.pdb").unwrap();
            let mut sel = src.select_str("within 1.0 of resid 560").unwrap();
            let mut cm = vec![];
            let trj = FileHandler::open("tests/protein.xtc").unwrap().into_iter();
            for st in trj {
                sel.set_state(st).unwrap();
                cm.push( sel.center_of_mass().unwrap() );
            }
            //println!("{:?}",&cm[..10]);
        }))
    );

    c.bench_function("trjconv", |b| b.iter(
        black_box(|| {
            let src = Source::serial_from_file("tests/protein.pdb").unwrap();
            let mut sel = src.select_str("resid 560").unwrap();

            let in_trj = FileHandler::open("tests/protein.xtc").unwrap().into_iter();
            let mut out_trj = FileHandler::create("target/.extracted.dcd").unwrap();
            for st in in_trj {
                sel.set_state(st).unwrap();
                out_trj.write_state(&sel).unwrap();
            }
        }))
    );    
}

criterion_group!{
    name = comparison_small;
    config = Criterion::default().sample_size(20).warm_up_time(Duration::from_secs(5));
    targets = molar_benchmark
}
criterion_main!(comparison_small);