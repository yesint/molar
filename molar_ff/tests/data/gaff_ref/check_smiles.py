#!/usr/bin/env python3
"""Cross-check molar's SMILES writer against RDKit over the 597-molecule GAFF corpus.

This is the real correctness gate for the writer. The in-tree Rust checks
(`molar_ff/tests/smiles_parity.rs`) verify that the emitted string encodes molar's own
molecular graph; only an independent toolkit can confirm the string *means* the same
molecule to everyone else. Not part of `cargo test` — RDKit is a development dependency,
in the same spirit as `molar/tests/gen_tests_pteros.py` and `gen_tests_vmd.tcl`.

Usage:
    cargo test -p molar_ff --test smiles_parity -- --ignored dump
    python tests/data/gaff_ref/check_smiles.py

Comparison is against **stereo-stripped** reference SMILES: `references.json` carries
RDKit-canonical strings with `[C@H]` / `/C=C/` markers, and molar's v1 writer emits no
stereochemistry at all. A mismatch here therefore means a real constitutional error —
wrong element, bond order, connectivity, charge, or hydrogen count — not a missing
stereocentre.

Two normalizations are needed to make that comparison mean what it should:

* **Explicit hydrogens are removed from both sides.** Some reference strings anchor an E/Z
  marker on an explicit hydrogen (`[H]/N=c1...`). Stripping stereochemistry leaves that
  hydrogen behind as a real graph atom, so the canonical forms would differ purely in how
  the same hydrogen is written.
* **Hydrogen-deficient inputs are reported separately.** A few corpus SDFs do not carry a
  complete hydrogen set. RDKit fills the gap from its valence model when it reads the file,
  so the reference string describes the *completed* molecule; molar has no implicit-hydrogen
  concept — hydrogens are ordinary atoms — so it faithfully writes the graph the file
  actually contains, using bracket atoms such as `[CH]` to state the real count. Neither is
  wrong, and lumping these in with genuine errors would hide the errors.
"""

import json
import os
import sys

try:
    from rdkit import Chem, RDLogger
except ImportError:
    sys.exit("RDKit is required: conda install -c conda-forge rdkit")

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))


def canonical(smiles):
    """Canonical SMILES with stereochemistry and explicit hydrogens normalized away, or None if
    RDKit cannot parse the string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    Chem.RemoveStereochemistry(mol)
    mol = Chem.RemoveHs(mol)
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def hydrogens_added_by_rdkit(name):
    """How many hydrogens RDKit's valence model invents when reading this SDF.

    Nonzero means the file is hydrogen-deficient, so the reference SMILES describes a molecule
    with hydrogens the file never had.
    """
    path = os.path.join(HERE, "sdf", "%s.sdf" % name)
    raw = Chem.MolFromMolFile(path, removeHs=False, sanitize=False)
    done = Chem.MolFromMolFile(path, removeHs=False)
    if raw is None or done is None:
        return 0
    in_file = sum(1 for a in raw.GetAtoms() if a.GetAtomicNum() == 1)
    completed = sum(
        a.GetTotalNumHs() + (1 if a.GetAtomicNum() == 1 else 0) for a in done.GetAtoms()
    )
    return completed - in_file


def main():
    with open(os.path.join(HERE, "references.json")) as f:
        reference = {m["name"]: m["smiles"] for m in json.load(f)["molecules"]}

    dump_path = os.path.join(HERE, "smiles_molar.json")
    if not os.path.exists(dump_path):
        sys.exit(
            "missing %s -- run:\n"
            "  cargo test -p molar_ff --test smiles_parity -- --ignored dump" % dump_path
        )
    with open(dump_path) as f:
        ours = {m["name"]: m["smiles"] for m in json.load(f)}

    agree, write_err, unparsable, mismatch, deficient = 0, [], [], [], []

    for name, smiles in sorted(ours.items()):
        if smiles.startswith("!ERROR"):
            write_err.append((name, smiles))
            continue

        got = canonical(smiles)
        if got is None:
            unparsable.append((name, smiles))
            continue

        want = canonical(reference[name])
        if want is None:
            # A bad reference string, not a writer failure -- report separately.
            unparsable.append((name, "REFERENCE: " + reference[name]))
            continue

        if got == want:
            agree += 1
        elif hydrogens_added_by_rdkit(name) > 0:
            deficient.append((name, hydrogens_added_by_rdkit(name)))
        else:
            mismatch.append((name, smiles, got, want))

    total = len(ours)
    comparable = total - len(deficient)
    print("\n=== molar SMILES vs RDKit over %d corpus molecules ===" % total)
    print("  agree            : %d/%d (%.2f%%)"
          % (agree, comparable, 100.0 * agree / max(comparable, 1)))
    print("  writer errors    : %d" % len(write_err))
    print("  RDKit unparsable : %d" % len(unparsable))
    print("  constitutional   : %d   <-- real failures" % len(mismatch))
    print("  H-deficient input: %d   (excluded: file lacks hydrogens the reference has)"
          % len(deficient))

    for name, err in write_err[:15]:
        print("  [write] %s: %s" % (name, err))
    for name, smiles in unparsable[:15]:
        print("  [parse] %s: %s" % (name, smiles))
    for name, added in deficient:
        print("  [defic] %s: RDKit invents %d hydrogen(s) the SDF does not contain"
              % (name, added))
    for name, raw, got, want in mismatch[:15]:
        print("  [diff]  %s\n            ours  %s\n            canon %s\n            want  %s"
              % (name, raw, got, want))
    dropped = len(write_err[15:]) + len(unparsable[15:]) + len(mismatch[15:])
    if dropped:
        print("  ... and %d more not shown" % dropped)
    print()

    return 0 if agree == comparable else 1


if __name__ == "__main__":
    sys.exit(main())
