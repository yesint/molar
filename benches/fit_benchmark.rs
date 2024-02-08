use criterion::{black_box, criterion_group, criterion_main, Criterion};
use lazy_static::lazy_static;
use molar::{core::{Topology, State, SelectionRc, Select, Vector3f, fit_transform, fit_transform_at_origin, ModifyPos}, io::FileHandler};
use nalgebra::Unit;

fn read_test_pdb() -> (Topology, State) {
    let mut h = FileHandler::open("tests/no_ATP.pdb").unwrap();
    let top = h.read_topology().unwrap();
    let state = h.read_state().unwrap().unwrap();
    (top, state)
}

// Read the test PDB file once and provide the content for tests
lazy_static! {
    static ref SS: (Topology, State) = read_test_pdb();
}

fn make_sel_prot() -> anyhow::Result<SelectionRc> {
    let t = SS.0.clone().to_rc();
    let s = SS.1.clone().to_rc();
    let sel = "not resname TIP3 POT CLA".select(&t, &s)?;
    Ok(sel)
}

fn test_fit(c: &mut Criterion) {
    let sel1 = make_sel_prot().unwrap();
    let sel2 = make_sel_prot().unwrap();   
    sel2.modify().rotate(&Unit::new_normalize(Vector3f::x()), 80.0_f32.to_radians());   
       
    //c.bench_function("fit gmx", |b| b.iter(
    //    || fit_transform_matrix(black_box(sel1.query().iter_particles()), sel2.query().iter_particles()).unwrap())
    //);
    
    //c.bench_function("fit quad", |b| b.iter(
    //    || fit_transform_gmx(black_box(sel1.query().iter_particles()), sel2.query().iter_particles()).unwrap())
    //);

    c.bench_function("fit kabsch ref", |b| b.iter(
        || fit_transform(black_box(sel1.query()), sel2.query()).unwrap())
    );

    c.bench_function("fit kabsch at origin", |b| b.iter(
        || fit_transform_at_origin(black_box(sel1.query()), sel2.query()).unwrap())
    );
}

criterion_group!(benches, test_fit);
criterion_main!(benches);