use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use molar::prelude::*;
use std::hint::black_box;

fn distance_search(c: &mut Criterion) {
    let mut group = c.benchmark_group("distance_search");
    for n in [32, 10_000, 100_000] {
        let side = (n as Float / 30.).cbrt().max(1.);
        let matrix = Matrix3f::from_diagonal(&Vector3f::repeat(side));
        let pbox = PeriodicBox::from_matrix(matrix).unwrap();
        let mut seed = 1_u64;
        let data = State {
            coords: (0..n)
                .map(|_| {
                    Pos::from(Vector3f::from_fn(|_, _| {
                        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                        ((seed >> 32) as u32 as Float / u32::MAX as Float) * side
                    }))
                })
                .collect(),
            ..State::default()
        };
        let image_shift = matrix * Vector3f::repeat(1.5);
        let wrapped_data = State {
            coords: data.coords.iter().map(|p| p + image_shift).collect(),
            ..State::default()
        };
        group.throughput(Throughput::Elements(n as u64));
        group.bench_with_input(BenchmarkId::new("single", n), &n, |b, &n| {
            b.iter(|| {
                distance_search_single::<(usize, usize), Vec<_>>(
                    black_box(0.5),
                    black_box(&data),
                    0..n,
                )
                .unwrap()
            })
        });
        group.bench_with_input(BenchmarkId::new("single_pbc", n), &n, |b, &n| {
            b.iter(|| {
                distance_search_single_pbc::<(usize, usize), Vec<_>>(
                    black_box(0.5),
                    data.iter_pos(),
                    0..n,
                    black_box(&pbox),
                    PBC_FULL,
                )
                .unwrap()
            })
        });
        group.bench_with_input(BenchmarkId::new("single_partial_pbc", n), &n, |b, &n| {
            b.iter(|| {
                distance_search_single_pbc::<(usize, usize), Vec<_>>(
                    black_box(0.5),
                    data.iter_pos(),
                    0..n,
                    black_box(&pbox),
                    PbcDims::new(true, true, false),
                )
                .unwrap()
            })
        });
        group.bench_with_input(BenchmarkId::new("single_wrapped_pbc", n), &n, |b, &n| {
            b.iter(|| {
                distance_search_single_pbc::<(usize, usize), Vec<_>>(
                    black_box(0.5),
                    wrapped_data.iter_pos(),
                    0..n,
                    black_box(&pbox),
                    PBC_FULL,
                )
                .unwrap()
            })
        });
    }
    group.finish();
}

criterion_group!(benches, distance_search);
criterion_main!(benches);
