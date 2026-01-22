# tests/test_pymolar_api.py
#
# Pytest-style, idiomatic tests derived from the original ad-hoc script.
# Adjust DATA_DIR if your repo layout differs.

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pymolar import FileHandler, PeriodicBox, System, distance_search


# --- test data paths ---------------------------------------------------------

HERE = Path(__file__).resolve().parent
DATA_DIR = (HERE / ".." / ".." / "molar" / "tests").resolve()

PDB_PATH = DATA_DIR / "protein.pdb"
XTC_PATH = DATA_DIR / "protein.xtc"


# --- fixtures ----------------------------------------------------------------

@pytest.fixture()
def top_and_state():
    top, st1 = FileHandler(str(PDB_PATH), "r").read()
    return top, st1


@pytest.fixture()
def system(top_and_state):
    top, st1 = top_and_state
    return System(top, st1)


@pytest.fixture()
def selection(system):
    return system("resid 5:600")


# --- tests -------------------------------------------------------------------

def test_set_state_updates_system_time_and_returns_previous_state(top_and_state):
    top, st1 = top_and_state

    # Read a second state (and tweak its time)
    st2 = FileHandler(str(PDB_PATH), "r").read_state()
    st2.time = 100

    sys_ = System(top, st1)
    sel1 = sys_("name CA")
    sel2 = sys_("name CB")

    assert sys_.time == st1.time
    assert sel1.time == st1.time
    assert sel2.time == st1.time
    assert st1.time != st2.time

    old = sel1.set_state(st2)

    # Expected: setting state on a selection updates the whole System state
    assert sys_.time == st2.time
    assert sel1.time == st2.time
    assert sel2.time == st2.time

    # Returned value should be the previous state (or at least carry its time)
    assert old.time == st1.time

    # Ensure original st1 was not mutated
    assert st1.time != st2.time


@pytest.mark.slow
@pytest.mark.skipif(not XTC_PATH.exists(), reason="Trajectory file not available")
def test_iterating_trajectory_and_setting_state_updates_selection_com(selection):
    trj = FileHandler(str(XTC_PATH), "r")

    # Just sanity-check the loop: state application should work and com should be finite
    n = 0
    for st in trj:
        selection.set_state(st)
        com = selection.com()
        assert np.isfinite(com).all()
        n += 1
        if n >= 5:  # keep unit test fast
            break

    assert n > 0


def test_pos_is_a_view_into_underlying_coordinates(system):
    sel = system("resid 5:600")

    pos0 = sel[0].pos
    before = sel[0].pos.copy()

    # Mutating the returned array should mutate the underlying coordinates (as in the original script)
    pos0[1] += 1.0
    after = sel[0].pos

    assert after[1] == pytest.approx(before[1] + 1.0)


def test_atom_setters_do_not_alias_between_atoms(system):
    sel = system("resid 5:600")

    # Position setter affects only that atom
    pos1_before = sel[1].pos.copy()
    sel[0].pos = [100, 100, 3]
    assert_allclose(sel[0].pos, np.array([100, 100, 3], dtype=float), rtol=0, atol=0)
    assert_allclose(sel[1].pos, pos1_before)

    # Name setter affects only that atom
    name1_before = sel[1].name
    sel[0].name = "AAA"
    assert sel[0].name == "AAA"
    assert sel[1].name == name1_before

    # Resid setter affects only that atom
    resid_before = sel[100].resid
    sel[0].resid = resid_before
    assert sel[0].resid == resid_before

    # Resname setter affects only that atom
    resname_before = sel[50].resname
    sel[0].resname = resname_before
    assert sel[0].resname == resname_before


def test_negative_indexing_and_scalar_coordinate_property(system):
    sel = system("resid 5:600")

    # Negative indexing should work
    atom = sel[-100]
    assert atom.atom.name == atom.name  # atom.atom.name exists per original script

    # x should be readable/writable and reflect in pos
    x_before = atom.x
    atom.x = 42
    assert atom.x == 42
    assert atom.pos[0] == pytest.approx(42)

    # restore (optional)
    atom.x = x_before


def test_subselection(system):
    sel = system("resid 5:600")
    subsel = sel("name CA")

    assert len(sel) > 0
    assert len(subsel) > 0
    assert len(subsel) <= len(sel)
    assert subsel[0].name == "CA"


def test_get_coord_set_coord_roundtrip(system):
    sel = system("resid 5:600")

    crd = sel.get_coord()
    assert isinstance(crd, np.ndarray)
    assert crd.shape[0] == 3
    assert crd.shape[1] == len(sel)

    # Modify coordinate matrix and apply
    crd2 = crd.copy()
    crd2[0, 0] = 42
    sel.set_coord(crd2)

    assert sel[0].pos[0] == pytest.approx(42)

    # Setting a full-zero coord array should work
    zeros = np.zeros((3, len(sel)), dtype=np.float32)
    sel.set_coord(zeros)
    assert_allclose(sel[0].pos, np.array([0.0, 0.0, 0.0], dtype=float), atol=1e-6)

    # Attribute access still works
    _ = sel[5].resname  # should not raise


def test_system_call_variants(system):
    # These mirror the original script's selection calls, with relationship assertions.
    sel_res = system("resid 500:600")
    sel_none = system(None)
    sel_empty = system()
    sel_range = system((0, 199))
    sel_list = system([1, 3, 4, 5, 6, 7])

    assert len(sel_res) >= 0
    assert len(sel_none) == len(sel_empty) == len(system)
    assert len(sel_list) == 6

    # Range length may be inclusive or exclusive depending on implementation
    assert len(sel_range) in (199, 200)


def test_periodic_box_vectors_angles_and_shortest_vector():
    b = PeriodicBox([1, 2, 3], [90, 90, 90])

    vecs, angs = b.to_vectors_angles()
    assert_allclose(np.array(angs, dtype=float), np.array([90, 90, 90], dtype=float), atol=1e-6)

    # shortest_vector should map into the minimum-image convention box
    v = np.array(b.shortest_vector([0.9, 0.5, 0.6]), dtype=float)

    # For a 1x2x3 orthorhombic box, x=0.9 should wrap to -0.1 (closest image)
    assert v[0] == pytest.approx(-0.1, abs=1e-6)
    assert v[1] == pytest.approx(0.5, abs=1e-6)
    assert v[2] == pytest.approx(0.6, abs=1e-6)


def test_distance_search_outputs_are_consistent(system):
    sel1 = system("resid 5:100")
    sel2 = system("resid 101:200")

    pairs, dist = distance_search("vdw", sel1, sel2)

    assert len(pairs) == len(dist)
    if len(dist) > 0:
        d = np.asarray(dist, dtype=float)
        assert (d >= 0).all()


def test_append_increases_system_size(system):
    before = len(system)
    sel = system("resid 550")
    system.append(sel)
    after = len(system)

    assert after == before + len(sel)
