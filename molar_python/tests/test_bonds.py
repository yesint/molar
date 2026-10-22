"""Bond table access at the Python interface: read, replace, validate, save."""

import numpy as np
import pytest

import pymolar
from pymolar import FileHandler, System

# Ethanol heavy atoms with two bonds, one of them with an explicit order.
CIF = """data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C C1 ETH A 1 0.000 0.000 0.000
ATOM 2 C C2 ETH A 1 1.520 0.000 0.000
ATOM 3 O O1 ETH A 1 2.000 1.300 0.000
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.pdbx_value_order
b1 covale ETH A 1 C1 ETH A 1 C2 sing
b2 covale ETH A 1 C2 ETH A 1 O1 ?
"""


@pytest.fixture
def ethanol(tmp_path):
    path = tmp_path / "ethanol.cif"
    path.write_text(CIF)
    return System(str(path))


def test_version_is_exported():
    assert isinstance(pymolar.__version__, str) and pymolar.__version__.count(".") == 2


def test_bonds_are_read_as_index_pairs(ethanol):
    bonds = ethanol.bonds
    assert bonds.shape == (2, 2)
    assert np.issubdtype(bonds.dtype, np.integer)
    assert sorted(map(sorted, bonds.tolist())) == [[0, 1], [1, 2]]
    np.testing.assert_array_equal(ethanol.topology.bonds, bonds)
    assert ethanol.bond_orders == ["single", "unspecified"]


def test_bonds_array_is_a_copy(ethanol):
    bonds = ethanol.bonds
    bonds[0] = [2, 2]
    assert sorted(map(sorted, ethanol.bonds.tolist())) == [[0, 1], [1, 2]]


def test_set_bonds_replaces_the_table(ethanol):
    ethanol.set_bonds([[0, 2]])
    assert ethanol.bonds.tolist() == [[0, 2]]
    assert ethanol.bond_orders is None  # no concrete order given
    ethanol.set_bonds(np.array([[0, 1], [1, 2]]), orders=["single", "double"])
    assert ethanol.bond_orders == ["single", "double"]
    ethanol.topology.set_bonds(np.zeros((0, 2), dtype=int))
    assert ethanol.bonds.shape == (0, 2)


def test_selections_survive_set_bonds(ethanol):
    sel = ethanol("name O1")
    ethanol.set_bonds([[0, 1]])
    assert [a.name for a in sel.iter_atoms()] == ["O1"]


@pytest.mark.parametrize(
    "pairs, orders, message",
    [
        ([[0, 3]], None, "bond 0-3"),  # index out of range
        ([[1, 1]], None, "bond 1-1"),  # self-bond
        ([[0, -1]], None, "negative"),
        ([[0, 1, 2]], None, "shape"),
        ([[0, 1]], ["single", "double"], "2 bond orders for 1 bonds"),
        ([[0, 1]], ["quadruple"], "unknown bond order"),
    ],
)
def test_set_bonds_rejects_invalid_tables(ethanol, pairs, orders, message):
    before = ethanol.bonds
    with pytest.raises(ValueError, match=message):
        ethanol.set_bonds(pairs, orders=orders)
    np.testing.assert_array_equal(ethanol.bonds, before)


@pytest.mark.parametrize("target", ["system", "selection"])
def test_saved_cif_keeps_bonds_and_orders(ethanol, tmp_path, target):
    ethanol.set_bonds([[0, 1], [1, 2]], orders=["single", "double"])
    out = tmp_path / "out.cif"
    writer = FileHandler(str(out), "w")
    writer.write(ethanol if target == "system" else ethanol("all"))
    del writer
    back = System(str(out))
    assert back.bonds.tolist() == [[0, 1], [1, 2]]
    assert back.bond_orders == ["single", "double"]
