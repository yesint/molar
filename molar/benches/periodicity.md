Periodicity benchmark and implementation decision

The implementation uses exact minimum-image distances for full and partial periodicity. Bits X, Y and Z select box vectors a, b and c. Nonperiodic grid coordinates use the data bounds; they are not restricted to the primary unit cell.

The box stores an LLL-reduced basis and a thin QR factorization for each two- or three-vector periodic lattice. Reduction preserves the lattice. Nearest-plane reduction supplies the initial image. A complete correction list tests opposite shifts together. Its coefficient bounds come from the reduced inverse basis and the diameter of the nearest-plane brick. These bounds replace the old fixed -1 to +1 range. If the list is too large, bounded sphere enumeration supplies the exact result. One-vector periodicity uses projection. Diagonal boxes have a separate fast path.

Run the benchmark from the workspace root:

```sh
RAYON_NUM_THREADS=1 cargo bench -p molar --bench periodicity -- --noplot
```

Recorded on 2026-09-09 with rustc 1.97.1, release optimization, f32, and an AMD Ryzen 7 PRO 6850U. No concurrent builds or tests were run during the final measurements. Criterion used 20 samples, 300 ms warm-up, and a 1 s measurement interval per case. Values below are medians. CPU frequency changes and system load can affect repeat runs.

Each distance sample processes 4,096 deterministic displacements. “Cell” uses uniform fractional coordinates in [-0.5, 0.5]. “Near” uses Cartesian components in [-0.5, 0.5] nm. “XY” uses the cell distribution with only a and b periodic.

Box vectors are columns:

- Diagonal: a=(10,0,0), b=(0,10,0), c=(0,0,10).
- Skew: a=(10,0,0), b=(4,10,0), c=(-4,0,10).
- Thin: a=(10,0,0), b=(0,10,0), c=(-2,0,1).

| Box / distribution | Previous method, ns/vector | Paired previous method, ns/vector | Adopted exact method, ns/vector | Previous / exact |
|---|---:|---:|---:|---:|
| diagonal / cell | 16.94 | 13.35 | 17.51 | 0.97× |
| diagonal / near | 16.94 | 13.54 | 17.51 | 0.97× |
| diagonal / xy | 15.43 | — | 15.56 | 0.99× |
| skew / cell | 78.33 | 43.31 | 51.75 | 1.51× |
| skew / near | 72.07 | 35.77 | 17.38 | 4.15× |
| skew / xy | 15.43 | — | 20.63 | 0.75× |
| thin / cell | 94.78 | 49.09 | 41.49 | 2.28× |
| thin / near | 75.92 | 38.69 | 10.48 | 7.25× |
| thin / xy | 15.50 | — | 17.70 | 0.88× |

The previous and paired previous methods are timing baselines, not correctness references. Both retain the incomplete fixed search range. The previous partial-PBC method omits image corrections. Its lower partial-PBC times do not represent equivalent work. The paired previous method was not adopted as-is.

The adopted method improves the tested full-periodic skew and thin cases. Diagonal times remain close to the previous method. Exact partial periodicity adds some cost.

| Box setup | Previous method, ns | Adopted method, ns |
|---|---:|---:|
| diagonal | 26.6 | 49.6 |
| skew | 345.5 | 1619.2 |
| thin | 357.2 | 784.2 |

The setup baseline omits paired-list construction and thus represents the original setup work. The adopted skew-box setup adds about 1.3 µs; the measured cell-displacement savings recover this cost after approximately 50 distance calls. Setup is repeated when box vectors change.

The grid benchmark uses 1,024 deterministic points in the skew box, a 1 nm cutoff, and one Rayon thread. Both methods allocate their output pair list. The grid measurement includes construction and population. The reference is direct testing of all pairs with the adopted exact distance routine.

| Periodicity | Grid, ms | Direct pair testing, ms | Direct / grid |
|---|---:|---:|---:|
| full | 1.507 | 38.855 | 25.8× |
| xy | 1.223 | 24.356 | 19.9× |

Validation completed:

- `cargo test -p molar`: 216 library tests and 7 documentation tests passed; 23 documentation tests remain ignored.
- `cargo check --workspace`: passed, with an existing unused-function warning in the Python crate.
- `cargo test -p molar --lib --features f64 periodic_box::tests`: 21 passed.
- `cargo test -p molar --lib --features f64 distance_search::tests`: 5 passed.

The new tests cover the two-vector-shift counterexample, all eight periodicity masks, permitted integer shifts, rotated boxes, large nonperiodic offsets, cache rebuilding after scaling, failed scaling, negative wrapping, face spacings, cached and fallback searches near image boundaries, the missed 0.8 nm grid pair, single/double/van der Waals searches, outside points, and duplicate prevention for small grids and large cutoffs.
