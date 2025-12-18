use std::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use molar::{prelude::*, bind};

const PDB: &str = "inp.pdb";
const XTC: &str = "traj_comp.xtc";

fn molar_benchmark_large(c: &mut Criterion) {
    c.bench_function("align_large", |b| {
        b.iter(black_box(|| {
            let mut sys = System::from_file(PDB).unwrap();
            let ref_ind = sys.select_as_index("protein").unwrap();
            let cur_ind = sys.select_as_index("protein").unwrap();

            let mut rmsds = vec![];

            let trj = FileHandler::open(XTC).unwrap().into_iter();
            for st in trj {
                sys.set_state(st).unwrap();
                let tr = bind!(sys, ref_ind, cur_ind, {
                    fit_transform(&cur_ind, &ref_ind).unwrap()
                });

                sys.bind_mut(&cur_ind).unwrap().apply_transform(&tr);

                bind!(sys, ref_ind, cur_ind, {
                    rmsds.push(rmsd(&cur_ind, &ref_ind).unwrap())
                });
            }
            //println!("{:?}",&rmsd[..10]);
        }))
    });

    c.bench_function("within_large", |b| {
        b.iter(black_box(|| {
            let mut sys = System::from_file(PDB).unwrap();
            let ind = sys.select_as_index("within 1.0 of protein").unwrap();
            let mut cm = vec![];
            let trj = FileHandler::open(XTC).unwrap().into_iter();
            for st in trj {
                sys.set_state(st).unwrap();
                cm.push(sys.bind(&ind).unwrap().center_of_mass().unwrap());
            }
            //println!("{:?}",&cm[..10]);
        }))
    });

    c.bench_function("trjconv_large", |b| {
        b.iter(black_box(|| {
            let mut sys = System::from_file(PDB).unwrap();
            let ind = sys.select_as_index("protein").unwrap();

            let in_trj = FileHandler::open(XTC).unwrap().into_iter();
            let mut out_trj = FileHandler::create("target/.extracted.dcd").unwrap();
            for st in in_trj {
                sys.set_state(st).unwrap();
                out_trj.write_state(&sys.bind(&ind).unwrap()).unwrap();
            }
        }))
    });
}

criterion_group! {
    name = comparison_large;
    config = Criterion::default().sample_size(10).warm_up_time(Duration::from_secs(5));
    targets = molar_benchmark_large
}
criterion_main!(comparison_large);
