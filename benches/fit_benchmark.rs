use criterion::{black_box, criterion_group, criterion_main, Criterion};
use molar::prelude::*;
use nalgebra::Unit;

fn read_test_pdb() -> (triomphe::UniqueArc<Topology>, triomphe::UniqueArc<State>) {
    let mut h = FileHandler::open("tests/colored.pdb").unwrap();
    let top = h.read_topology().unwrap();
    let state = h.read_state().unwrap().unwrap();
    (top, state)
}


fn make_sel_prot() -> anyhow::Result<Sel<MutableSerial>> {
    let (top,st) = read_test_pdb();
    let mut b = Source::new_serial(top, st).unwrap();
    //let sel = b.select_str("not resname TIP3 POT CLA").unwrap();
    let sel = b.select_all().unwrap();
    Ok(sel)
}

#[allow(dead_code)]
fn test_fit(c: &mut Criterion) {
    let sel1 = make_sel_prot().unwrap();
    let sel2 = make_sel_prot().unwrap();   
    sel2.rotate(&Unit::new_normalize(Vector3f::x()), 80.0_f32.to_radians());   
       
    //c.bench_function("fit gmx", |b| b.iter(
    //    || fit_transform_matrix(black_box(sel1.query().iter_particles()), sel2.query().iter_particles()).unwrap())
    //);
    
    //c.bench_function("fit quad", |b| b.iter(
    //    || fit_transform_gmx(black_box(sel1.query().iter_particles()), sel2.query().iter_particles()).unwrap())
    //);

    c.bench_function("fit kabsch ref", |b| b.iter(
        || Sel::fit_transform(black_box(&sel1), &sel2).unwrap())
    );

    c.bench_function("fit kabsch at origin", |b| b.iter(
        || Sel::fit_transform_at_origin(black_box(&sel1), &sel2).unwrap())
    );
}

fn search_par(c: &mut Criterion) {
    let sel = make_sel_prot().unwrap();

    let mut searcher = DistanceSearcherSingle::new_periodic(
        0.3, 
        sel.iter_particle().map(|p| (p.id,*p.pos)), 
        sel.get_box().unwrap(), 
        &PBC_FULL
    );

    c.bench_function("serial search_single", |b| b.iter(
        || {
            searcher.set_serial_limit(1e10 as usize);
            let _: Vec<usize> = black_box(searcher.search());
        })
    );

    c.bench_function("parallel search_single", |b| b.iter(
        || {
            searcher.set_serial_limit(0);
            let _: Vec<usize> = black_box(searcher.search());
        })
    );
}

fn paper_benchmark_align(c: &mut Criterion) {
    c.bench_function("align", |b| b.iter(
        black_box(
        || {
            let mut src = Source::serial_from_file("tests/protein.pdb").unwrap();
            let ref_sel = src.select_all().unwrap();
            let mut cur_sel = src.select_all().unwrap();

            let mut rmsd = vec![];

            let trj = FileHandler::open("tests/protein.xtc").unwrap().into_iter();
            for st in trj {
                cur_sel.set_shared_state(st.shareable()).unwrap();
                let tr = MeasureMasses::fit_transform(&cur_sel, &ref_sel).unwrap();
                cur_sel.apply_transform(&tr);
                rmsd.push( MeasurePos::rmsd(&cur_sel, &ref_sel).unwrap() );
            }
            //println!("{:?}",&rmsd[..10]);
        }))
    );

    c.bench_function("within", |b| b.iter(
        black_box(|| {
            let mut src = Source::serial_from_file("tests/protein.pdb").unwrap();
            let mut sel = src.select_str("within 1.0 of resid 560").unwrap();
            let mut cm = vec![];
            let trj = FileHandler::open("tests/protein.xtc").unwrap().into_iter();
            for st in trj {
                sel.set_shared_state(st.shareable()).unwrap();
                cm.push( sel.center_of_mass().unwrap() );
            }
            //println!("{:?}",&cm[..10]);
        }))
    );

    c.bench_function("trjconv", |b| b.iter(
        black_box(|| {
            let mut src = Source::serial_from_file("tests/protein.pdb").unwrap();
            let mut sel = src.select_str("resid 560").unwrap();

            let in_trj = FileHandler::open("tests/protein.xtc").unwrap().into_iter();
            let mut out_trj = FileHandler::create("tests/.extracted.xtc").unwrap();
            for st in in_trj {
                sel.set_shared_state(st.shareable()).unwrap();
                out_trj.write_state(&sel).unwrap();
            }
        }))
    );
}

criterion_group!{
    name = benches;
    // This can be any expression that returns a `Criterion` object.
    config = Criterion::default().sample_size(10);
    targets = paper_benchmark_align
}
criterion_main!(benches);