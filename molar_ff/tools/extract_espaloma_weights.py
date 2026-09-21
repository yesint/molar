#!/usr/bin/env python3
"""Extract the fixed espaloma-charge tensors from the source ONNX model.

The production Rust code does not parse ONNX. This maintenance tool reads only
the small protobuf subset that is necessary to copy named, little-endian f32
initializers. It has no third-party Python dependencies.

Run from the workspace root after a deliberate model update:

    python molar_ff/tools/extract_espaloma_weights.py

The script rejects a model whose tensor names, shapes, types, or byte counts do
not match the network implemented in ``molar_ff/src/charge.rs``. After it runs,
run the charge fixture and corpus tests. Do not edit the generated file by hand.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
import struct


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "assets" / "espaloma_charge.onnx"
OUTPUT = ROOT / "assets" / "espaloma_charge.weights"
MAGIC = b"ESPCHG01"
EXPECTED_SOURCE_SHA256 = "0ea0dc48742931980854836e6a2546a107890911ccba27750576ab04a0da0d83"

# Output order is the execution order. ONNX stores these tensors in a different
# order, so the extractor also removes all graph-format and Gather-index data.
TENSORS = (
    ("permute", (116, 128)),
    ("fin_b", (128,)),
    ("permute_1", (128, 128)),
    ("permute_2", (128, 128)),
    ("permute_3", (128, 128)),
    ("permute_4", (128, 128)),
    ("permute_5", (128, 128)),
    ("permute_6", (128, 128)),
    ("permute_7", (128, 128)),
    ("permute_8", (128, 128)),
    ("permute_9", (128, 2)),
    ("fr_b", (2,)),
)


def read_varint(data: bytes, offset: int) -> tuple[int, int]:
    """Read one unsigned protobuf varint."""
    value = 0
    shift = 0
    while True:
        if offset >= len(data) or shift >= 70:
            raise ValueError("invalid protobuf varint")
        byte = data[offset]
        offset += 1
        value |= (byte & 0x7F) << shift
        if byte < 0x80:
            return value, offset
        shift += 7


def fields(data: bytes):
    """Yield ``(field_number, wire_type, value)`` for one protobuf message."""
    offset = 0
    while offset < len(data):
        key, offset = read_varint(data, offset)
        number, wire_type = key >> 3, key & 7
        if wire_type == 0:
            value, offset = read_varint(data, offset)
        elif wire_type == 1:
            value = data[offset : offset + 8]
            offset += 8
        elif wire_type == 2:
            size, offset = read_varint(data, offset)
            value = data[offset : offset + size]
            offset += size
        elif wire_type == 5:
            value = data[offset : offset + 4]
            offset += 4
        else:
            raise ValueError(f"unsupported protobuf wire type {wire_type}")
        yield number, wire_type, value


def one_bytes(message: bytes, tag: int) -> bytes:
    values = [value for number, wire, value in fields(message) if number == tag and wire == 2]
    if len(values) != 1:
        raise ValueError(f"expected one length-delimited field {tag}, got {len(values)}")
    return values[0]


def unpack_packed_varints(data: bytes) -> tuple[int, ...]:
    values = []
    offset = 0
    while offset < len(data):
        value, offset = read_varint(data, offset)
        values.append(value)
    return tuple(values)


def main() -> None:
    model = SOURCE.read_bytes()
    source_sha256 = hashlib.sha256(model).hexdigest()
    if source_sha256 != EXPECTED_SOURCE_SHA256:
        raise ValueError(
            "source model checksum changed; review the graph and fixed Rust kernel, "
            "then update EXPECTED_SOURCE_SHA256 deliberately"
        )
    graph = one_bytes(model, 7)  # ModelProto.graph
    initializers: dict[str, tuple[tuple[int, ...], int, bytes]] = {}
    for number, wire, tensor in fields(graph):
        if number != 5 or wire != 2:  # GraphProto.initializer
            continue
        name = one_bytes(tensor, 8).decode("utf-8")
        dims_values = []
        for field, kind, value in fields(tensor):
            if field == 1 and kind == 0:
                dims_values.append(value)
            elif field == 1 and kind == 2:
                dims_values.extend(unpack_packed_varints(value))
        dims = tuple(dims_values)
        data_type = next(
            value for field, kind, value in fields(tensor) if field == 2 and kind == 0
        )
        raw = one_bytes(tensor, 9)
        initializers[name] = (dims, data_type, raw)

    output = bytearray(MAGIC)
    for name, expected_dims in TENSORS:
        try:
            dims, data_type, raw = initializers[name]
        except KeyError as error:
            raise ValueError(f"source model has no initializer {name!r}") from error
        if dims != expected_dims:
            raise ValueError(f"{name}: expected shape {expected_dims}, got {dims}")
        if data_type != 1:  # TensorProto.FLOAT
            raise ValueError(f"{name}: expected f32 data type 1, got {data_type}")
        expected_bytes = 4
        for dim in dims:
            expected_bytes *= dim
        if len(raw) != expected_bytes:
            raise ValueError(f"{name}: expected {expected_bytes} bytes, got {len(raw)}")
        # Decode and encode once to reject partial or non-little-endian data.
        count = len(raw) // 4
        output.extend(struct.pack(f"<{count}f", *struct.unpack(f"<{count}f", raw)))

    OUTPUT.write_bytes(output)
    output_sha256 = hashlib.sha256(output).hexdigest()
    print(
        f"wrote {OUTPUT} ({len(output)} bytes, {len(output[8:]) // 4} f32 values, "
        f"sha256={output_sha256})"
    )


if __name__ == "__main__":
    main()
