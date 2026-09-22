"""Bond/hydrogen perception at the Python interface."""

import math

import pytest
from pymolar import System


def _benzene_heavy_xyz(tmp_path):
    """Six bare carbons on a hexagon (Angstrom; the .xyz reader scales to nm)."""
    r = 1.39
    lines = ["6", "benzene heavy atoms"]
    for k in range(6):
        t = math.pi / 3.0 * k
        lines.append(f"C {r * math.cos(t):.4f} {r * math.sin(t):.4f} 0.0000")
    path = tmp_path / "benzene.xyz"
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_perception_pipeline_step_by_step(tmp_path):
    sys = System(_benzene_heavy_xyz(tmp_path))
    assert len(sys) == 6

    assert sys.perceive_connectivity() == 6  # a six-membered ring
    warnings = sys.perceive_bond_orders(infer_hydrogens=True)
    assert isinstance(warnings, list)
    assert sys.add_hydrogens() == 6  # C6H6
    assert len(sys) == 12


def test_prepare_for_ff_enables_typing(tmp_path):
    sys = System(_benzene_heavy_xyz(tmp_path))
    sys.prepare_for_ff(infer_hydrogens=True, add_hydrogens=True)
    assert len(sys) == 12
    # A raw structure prepared this way can now be GAFF-typed without error.
    sys.apply_ff("gaff")


def test_apply_ff_rejects_order_less_bonds(tmp_path):
    sys = System(_benzene_heavy_xyz(tmp_path))
    sys.perceive_connectivity()  # bonds now exist but have unspecified order
    # Strict typing refuses order-less bonds; preparation is required first.
    with pytest.raises(ValueError):
        sys.apply_ff("gaff")
