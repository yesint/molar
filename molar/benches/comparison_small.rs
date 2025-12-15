use std::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use molar::prelude::*;

fn molar_benchmark(c: &mut Criterion) {
    c.bench_function("align", |b| {
        b.iter(black_box(|| {
            let mut sys = System::from_file("tests/protein.pdb").unwrap();
            let ref_ind = sys.select_all_as_index();
            let cur_ind = sys.select_all_as_index();

            let mut rmsds = vec![];

            let trj = FileHandler::open("tests/protein.xtc").unwrap().into_iter();
            for st in trj {
                sys.set_state(st).unwrap();
                let tr = fit_transform(&sys.bind(&cur_ind).unwrap(), &sys.bind(&ref_ind).unwrap())
                    .unwrap();
                sys.bind_mut(&cur_ind).unwrap().apply_transform(&tr);
                rmsds.push(
                    rmsd(&sys.bind(&cur_ind).unwrap(), &sys.bind(&ref_ind).unwrap()).unwrap(),
                );
            }
            //println!("{:?}",&rmsd[..10]);
        }))
    });

    c.bench_function("within", |b| {
        b.iter(black_box(|| {
            let mut sys = System::from_file("tests/protein.pdb").unwrap();
            let ind = sys.select_as_index("within 1.0 of resid 560").unwrap();
            let mut cm = vec![];
            let trj = FileHandler::open("tests/protein.xtc").unwrap().into_iter();
            for st in trj {
                sys.set_state(st).unwrap();
                cm.push(sys.bind(&ind).unwrap().center_of_mass().unwrap());
            }
            //println!("{:?}",&cm[..10]);
        }))
    });

    c.bench_function("trjconv", |b| {
        b.iter(black_box(|| {
            let mut sys = System::from_file("tests/protein.pdb").unwrap();
            let ind = sys.select_as_index("resid 560").unwrap();

            let in_trj = FileHandler::open("tests/protein.xtc").unwrap().into_iter();
            let mut out_trj = FileHandler::create("target/.extracted.dcd").unwrap();
            for st in in_trj {
                sys.set_state(st).unwrap();
                out_trj.write_state(&sys.bind(&ind).unwrap()).unwrap();
            }
        }))
    });
}

criterion_group! {
    name = comparison_small;
    config = Criterion::default().sample_size(20).warm_up_time(Duration::from_secs(5));
    targets = molar_benchmark
}
criterion_main!(comparison_small);
