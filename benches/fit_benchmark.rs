use criterion::{black_box, criterion_group, criterion_main, Criterion};
use molar::{core::{providers::BoxProvider, MeasureMasses, ModifyPos, OverlappingMut, Sel, Source, State, Topology, Vector3f, PBC_FULL}, distance_search::DistanceSearcherSingle, io::FileHandler};
use nalgebra::Unit;

fn read_test_pdb() -> (triomphe::UniqueArc<Topology>, triomphe::UniqueArc<State>) {
    let mut h = FileHandler::open("tests/colored.pdb").unwrap();
    let top = h.read_topology_raw().unwrap().to_rc();
    let state = h.read_state_raw().unwrap().unwrap().to_rc();
    (top, state)
}


fn make_sel_prot() -> anyhow::Result<Sel<OverlappingMut>> {
    let (top,st) = read_test_pdb();
    let mut b = Source::new_overlapping_mut(top, st)?;
    //let sel = b.select_str("not resname TIP3 POT CLA")?;
    let sel = b.select_all()?;
    Ok(sel)
}

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
        sel.iter().map(|(i,_,p)| (i,p)), 
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

criterion_group!(benches_par, search_par);
criterion_main!(benches_par);