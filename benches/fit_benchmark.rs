use criterion::{black_box, criterion_group, criterion_main, Criterion};
use molar::{core::{MeasureMasses, ModifyPos, Sel, SelBuilder, Serial, State, Topology, Vector3f}, io::FileHandler};
use nalgebra::Unit;

fn read_test_pdb() -> (triomphe::UniqueArc<Topology>, triomphe::UniqueArc<State>) {
    let mut h = FileHandler::open("tests/no_ATP.pdb").unwrap();
    let top = h.read_topology_raw().unwrap().to_rc();
    let state = h.read_state_raw().unwrap().unwrap().to_rc();
    (top, state)
}


fn make_sel_prot() -> anyhow::Result<Sel<Serial>> {
    let (top,st) = read_test_pdb();
    let b = SelBuilder::new(top, st)?;
    let sel = b.select_str("not resname TIP3 POT CLA")?;
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

criterion_group!(benches, test_fit);
criterion_main!(benches);