use std::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use molar::prelude::*;

const PDB: &str = "inp.pdb";
const XTC: &str = "traj_comp.xtc";

fn molar_benchmark_large(c: &mut Criterion) {
    c.bench_function("align_large", |b| b.iter(
        black_box(
        || {
            let src = Source::serial_from_file(PDB).unwrap();
            let ref_sel = src.select_str("protein").unwrap();
            let mut cur_sel = src.select_str("protein").unwrap();

            let mut rmsd = vec![];

            let trj = FileHandler::open(XTC).unwrap().into_iter();
            for st in trj {
                cur_sel.set_state(st).unwrap();
                let tr = MeasureMasses::fit_transform(&cur_sel, &ref_sel).unwrap();
                cur_sel.apply_transform(&tr);
                rmsd.push( MeasurePos::rmsd(&cur_sel, &ref_sel).unwrap() );
            }
            //println!("{:?}",&rmsd[..10]);
        }))
    );

    c.bench_function("within_large", |b| b.iter(
        black_box(|| {
            let src = Source::serial_from_file(PDB).unwrap();
            let mut sel = src.select_str("within 1.0 of protein").unwrap();
            let mut cm = vec![];
            let trj = FileHandler::open(XTC).unwrap().into_iter();
            for st in trj {
                sel.set_state(st).unwrap();
                cm.push( sel.center_of_mass().unwrap() );
            }
            //println!("{:?}",&cm[..10]);
        }))
    );

    c.bench_function("trjconv_large", |b| b.iter(
        black_box(|| {
            let src = Source::serial_from_file(PDB).unwrap();
            let mut sel = src.select_str("protein").unwrap();

            let in_trj = FileHandler::open(XTC).unwrap().into_iter();
            let mut out_trj = FileHandler::create("target/.extracted.dcd").unwrap();
            for st in in_trj {
                sel.set_state(st).unwrap();
                out_trj.write_state(&sel).unwrap();
            }
        }))
    );    
}

criterion_group!{
    name = comparison_large;
    config = Criterion::default().sample_size(10).warm_up_time(Duration::from_secs(5));
    targets = molar_benchmark_large
}
criterion_main!(comparison_large);