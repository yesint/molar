"""LAMMPS data support through the Python interface."""
from pathlib import Path

import numpy as np
import pytest
import pymolar as mol

DATA = Path(__file__).resolve().parents[2] / "molar" / "tests" / "lammps"


def _bond_count(path):
    for line in path.read_text().splitlines():
        if line.endswith(" bonds"):
            return int(line.split()[0])
    raise AssertionError("missing bond count")


def test_lammps_read_scales_and_tuple_write(tmp_path):
    source = DATA / "wrapped_images.data"
    reader = mol.FileHandler(str(source), "r", lammps_length_scale=0.1,
                             lammps_mass_scale=2.0)
    top, state = reader.read()
    system = mol.System(top, state)
    original = mol.System(str(source))
    np.testing.assert_allclose(system().coords, original().coords * 0.1, atol=1e-5)
    output = tmp_path / "scaled.data"
    writer = mol.FileHandler(str(output), "w", lammps_length_scale=0.1,
                             lammps_mass_scale=2.0)
    writer.write((top, state))
    assert _bond_count(output) == 400
    np.testing.assert_allclose(mol.System(str(output))().coords,
                               original().coords, atol=1e-4)
    with pytest.raises(Exception):
        writer.write((top, state))


def test_lammps_system_and_selection_save_keep_bonds(tmp_path):
    system = mol.System(str(DATA / "unwrapped.data"))
    output = tmp_path / "system.data"
    system.save(str(output))
    assert _bond_count(output) == 400
    np.testing.assert_allclose(mol.System(str(output))().coords,
                               system().coords, atol=1e-4)
    output = tmp_path / "selection.data"
    selection = system("resid 2")
    mol.FileHandler(str(output), "w").write(selection)
    assert _bond_count(output) == 20
    reread = mol.System(str(output))
    assert len(reread) == 21
    np.testing.assert_allclose(reread().coords, selection.coords, atol=1e-4)


def test_lammps_scale_errors(tmp_path):
    output = tmp_path / "invalid.data"
    with pytest.raises(Exception):
        mol.FileHandler(str(output), "w", lammps_length_scale=0.0)
    assert not output.exists()
    with pytest.raises(ValueError):
        mol.FileHandler(str(tmp_path / "output.pdb"), "w", lammps_mass_scale=1.0)


@pytest.mark.parametrize("extension", ["sdf", "cif"])
def test_selection_bonds_in_other_formats(tmp_path, extension):
    # A selected middle fragment must keep its two bonds and their orders.
    source = tmp_path / "chain.sdf"
    source.write_text("""chain
  test

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  3  0  0  0  0
  4  5  1  0  0  0  0
M  END
$$$$
""")
    system = mol.System(str(source))
    selection = system([1, 2, 3])
    # CIF bond records use atom names to identify their endpoints.
    for i, atom in enumerate(system().iter_atoms()):
        atom.name = f"C{i}"
    output = tmp_path / f"fragment.{extension}"
    selection.save(str(output))
    reread = mol.System(str(output))
    assert len(reread) == 3
    np.testing.assert_allclose(reread().coords, selection.coords, atol=1e-5)
    # Convert back to SDF to check the retained bond types and endpoints.
    check = tmp_path / "check.sdf"
    reread.save(str(check))
    lines = check.read_text().splitlines()
    assert int(lines[3][3:6]) == 2
    assert [(int(line[0:3]), int(line[3:6]), int(line[6:9]))
            for line in lines[7:9]] == [(1, 2, 2), (2, 3, 3)]
