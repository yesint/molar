use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use molar::prelude::*;
use std::{hint::black_box, time::Duration};
#[path = "support/legacy_pbc.rs"]
mod legacy_pbc;

fn periodicity(c: &mut Criterion) {
    let boxes = [
        (
            "diagonal",
            Matrix3f::from_diagonal(&Vector3f::new(10., 10., 10.)),
        ),
        (
            "skew",
            Matrix3f::new(10., 4., -4., 0., 10., 0., 0., 0., 10.),
        ),
        ("thin", Matrix3f::new(10., 0., -2., 0., 10., 0., 0., 0., 1.)),
    ];
    for (name, m) in boxes {
        let b = PeriodicBox::from_matrix(m).unwrap();
        let old = legacy_pbc::LegacyBox::new(m);
        for (distribution, pbc) in [
            ("cell", PBC_FULL),
            ("near", PBC_FULL),
            ("xy", PbcDims::new(true, true, false)),
        ] {
            let mut seed = 2026_u64;
            let data: Vec<_> = (0..4096)
                .map(|_| {
                    let v = Vector3f::from_fn(|_, _| {
                        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                        (seed >> 32) as u32 as Float / u32::MAX as Float - 0.5
                    });
                    if distribution == "near" { v } else { m * v }
                })
                .collect();
            let mut group = c.benchmark_group(format!("periodicity/{name}/{distribution}"));
            group.throughput(Throughput::Elements(data.len() as u64));
            group.bench_function("legacy", |bencher| {
                bencher.iter(|| {
                    for v in &data {
                        black_box(black_box(&old).shortest(black_box(v), black_box(pbc)));
                    }
                })
            });
            if pbc == PBC_FULL {
                group.bench_function("paired_legacy", |bencher| {
                    bencher.iter(|| {
                        for v in &data {
                            black_box(black_box(&old).paired(black_box(v)));
                        }
                    })
                });
            }
            group.bench_function("exact", |bencher| {
                bencher.iter(|| {
                    for v in &data {
                        black_box(black_box(&b).shortest_vector_dims(black_box(v), black_box(pbc)));
                    }
                })
            });
            group.finish();
        }
    }
}

fn distance_grid(c: &mut Criterion) {
    let m = Matrix3f::new(10., 4., -4., 0., 10., 0., 0., 0., 10.);
    let b = PeriodicBox::from_matrix(m).unwrap();
    let mut seed = 417_u64;
    let points: Vec<_> = (0..1024)
        .map(|_| {
            Pos::from(
                m * Vector3f::from_fn(|_, _| {
                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                    (seed >> 32) as u32 as Float / u32::MAX as Float
                }),
            )
        })
        .collect();
    for (name, pbc) in [("full", PBC_FULL), ("xy", PbcDims::new(true, true, false))] {
        let mut group = c.benchmark_group(format!("distance_grid/{name}"));
        group.bench_function("grid", |bencher| {
            bencher.iter(|| {
                let pairs: Vec<(usize, usize)> = distance_search_single_pbc(
                    black_box(1.0),
                    points.iter(),
                    0..points.len(),
                    black_box(&b),
                    black_box(pbc),
                );
                black_box(pairs);
            })
        });
        group.bench_function("direct", |bencher| {
            bencher.iter(|| {
                let mut pairs = Vec::new();
                for i in 0..points.len() {
                    for j in i + 1..points.len() {
                        if b.distance_squared(
                            black_box(&points[i]),
                            black_box(&points[j]),
                            black_box(pbc),
                        ) <= 1.0
                        {
                            pairs.push((i, j));
                        }
                    }
                }
                black_box(pairs);
            })
        });
        group.finish();
    }
}

fn box_setup(c: &mut Criterion) {
    for (name, m) in [
        ("diagonal", Matrix3f::from_diagonal(&Vector3f::repeat(10.))),
        (
            "skew",
            Matrix3f::new(10., 4., -4., 0., 10., 0., 0., 0., 10.),
        ),
        ("thin", Matrix3f::new(10., 0., -2., 0., 10., 0., 0., 0., 1.)),
    ] {
        let mut group = c.benchmark_group(format!("box_setup/{name}"));
        group.bench_function("legacy", |b| {
            b.iter(|| black_box(legacy_pbc::LegacyBox::original(black_box(m))))
        });
        group.bench_function("exact", |b| {
            b.iter(|| black_box(PeriodicBox::from_matrix(black_box(m)).unwrap()))
        });
        group.finish();
    }
}

criterion_group! {name=benches;config=Criterion::default().sample_size(20).warm_up_time(Duration::from_millis(300)).measurement_time(Duration::from_secs(1));targets=periodicity,distance_grid,box_setup}
criterion_main!(benches);
