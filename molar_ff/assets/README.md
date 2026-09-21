# Espaloma charge-model assets

`espaloma_charge.onnx` is the source model. It is retained for provenance and
for an explicit model update. It is not linked into `molar_ff` and is not parsed
at run time.

`espaloma_charge.weights` is the generated production asset. Its first eight
bytes are the ASCII format marker `ESPCHG01`. The rest is little-endian `f32`
data. The tensors occur in this fixed order:

1. input matrix `[116, 128]`;
2. input bias `[128]`;
3. four pairs of self and neighbor matrices `[128, 128]`;
4. output matrix `[128, 2]`;
5. output bias `[2]`.

The generated file is 585,232 bytes. The reviewed asset checksums are:

- source ONNX: `0ea0dc48742931980854836e6a2546a107890911ccba27750576ab04a0da0d83`;
- generated weights: `c5b93ad94dea3e28f61e4e98c68e57aadace98aea9df9bcbab6f0786111e416c`.

The fixed layout is intentional. The model has only `MatMul`, `Add`, `Tanh`,
`Relu`, and output-column selection operations. A general ONNX parser and graph
executor added a large dependency and compile-time cost without adding useful
run-time flexibility.

## Regeneration

Run this command from the workspace root:

```text
python molar_ff/tools/extract_espaloma_weights.py
```

The script uses only the Python standard library. It pins the full source file
checksum and checks all tensor names, shapes, types, and byte counts. A model
change must cause a corresponding review of `molar_ff/src/charge.rs` and an
explicit update of the checksum in the extractor.

After regeneration, run:

```text
cargo test -p molar_ff
cargo test -p molar_ff large_sparse_inference_smoke --release -- --ignored --nocapture
```

The first command checks the Python reference fixture and the full charge
corpus. The second command is an opt-in large sparse-system check. It reports
elapsed time and verifies finite output without making a timing assertion.
