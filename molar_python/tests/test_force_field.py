"""Force-field API parity and option/error behavior for both precision builds.

Install both wheels before running this file. PYMOLAR_TEST_MODULE can select
one installed package for a single-wheel job.
"""

import importlib
import math
import os

import pytest

PACKAGES = ([os.environ["PYMOLAR_TEST_MODULE"]] if "PYMOLAR_TEST_MODULE" in os.environ
            else ["pymolar", "pymolar_f64"])


@pytest.fixture(params=PACKAGES)
def mol(request):
    return importlib.import_module(request.param)


@pytest.fixture
def benzene(tmp_path):
    path = tmp_path / "benzene.xyz"
    lines = ["6", "benzene heavy atoms"]
    for k in range(6):
        angle = math.pi * k / 3
        lines.append(f"C {1.39 * math.cos(angle):.4f} {1.39 * math.sin(angle):.4f} 0")
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def prepare(mol, path):
    sys = mol.System(path)
    opts = mol.PrepareOptions(
        bond_orders=mol.BondOrderOptions(
            input_orders=mol.InputOrders.PreserveKnown,
            hydrogens=mol.HydrogenPolicy.InferFromGeometry,
            total_charge=0,
        ),
        add_hydrogens=mol.HydrogenOptions(),
    )
    sys.prepare_for_ff(options=opts)
    return sys


def test_nested_defaults_and_shared_options(mol):
    first, second = mol.PrepareOptions(), mol.PrepareOptions()
    assert first.bond_orders.input_orders == mol.InputOrders.PreserveKnown
    assert mol.BondOrderOptions().input_orders == mol.InputOrders.ReassignAll
    assert first.bond_orders.hydrogens == mol.HydrogenPolicy.AllExplicit
    assert first.connectivity.pbc == [False, False, False]
    assert first.add_hydrogens is None
    first.connectivity.tolerance = 0.08
    first.bond_orders.limits.max_branches = 123
    assert second.connectivity.tolerance == pytest.approx(0.045)
    assert second.bond_orders.limits.max_branches == 5_000_000
    supplied = mol.SearchLimits(max_branches=456)
    first.bond_orders.limits = supplied
    supplied.max_branches = 789
    assert first.bond_orders.limits.max_branches == 789
    first.add_hydrogens = mol.HydrogenOptions(zero_fill_dynamics=False)
    assert first.add_hydrogens.zero_fill_dynamics is False


@pytest.mark.parametrize("target", ["system", "topology", "selection"])
@pytest.mark.parametrize("ff_name", ["gaff", "gaff2"])
def test_enum_and_string_typing_and_charging(mol, benzene, target, ff_name):
    sys = prepare(mol, benzene)
    assert len(sys) == 12
    obj = {"system": sys, "topology": sys.topology, "selection": sys()}[target]
    obj.apply_ff(ff_name.upper())
    expected = [a.type_name for a in sys.iter_atoms()]
    assert all(expected)
    enum = mol.FFType.Gaff if ff_name == "gaff" else mol.FFType.Gaff2
    obj.apply_ff(enum)
    assert [a.type_name for a in sys.iter_atoms()] == expected
    obj.apply_charges("ESPALOMA")
    expected = [a.charge for a in sys.iter_atoms()]
    obj.apply_charges(mol.ChargeModel.Espaloma)
    assert [a.charge for a in sys.iter_atoms()] == pytest.approx(expected)
    assert sum(expected) == pytest.approx(0, abs=2e-5)
    assert sys[0].atom.type_name is not None  # topology edits reach shared atom views


@pytest.mark.parametrize("target", ["system", "topology", "selection"])
def test_typed_errors_preserve_value_error_and_fields(mol, benzene, target):
    sys = mol.System(benzene)
    sys.perceive_connectivity()
    obj = {"system": sys, "topology": sys.topology, "selection": sys()}[target]
    for method, error_type in [(obj.apply_ff, mol.FFError), (obj.apply_charges, mol.ChargeError)]:
        with pytest.raises(error_type) as caught:
            method()
        assert isinstance(caught.value, ValueError)
        assert caught.value.kind == "MissingBondOrders"
        assert len(caught.value.details["atoms"]) == 2
    assert len(sys) == 6


@pytest.mark.parametrize("method,error_name", [("apply_ff", "FFError"), ("apply_charges", "ChargeError")])
def test_open_selection_reports_global_boundary(mol, benzene, method, error_name):
    sys = prepare(mol, benzene)
    with pytest.raises(getattr(mol, error_name)) as caught:
        getattr(sys([2]), method)()
    assert caught.value.kind == "OpenSelection"
    assert caught.value.details["atom"] == 2
    assert caught.value.details["neighbor"] != 2


def test_missing_box_option_is_typed_and_data_is_restored(mol, benzene):
    sys = mol.System(benzene)
    before = sys().coords.copy()
    opts = mol.PrepareOptions(connectivity=mol.ConnectivityOptions(pbc=[True, True, False]))
    with pytest.raises(mol.BondPerceptionError) as caught:
        sys.prepare_for_ff(options=opts)
    assert caught.value.kind == "MissingPeriodicBox"
    assert caught.value.details == {}
    assert len(sys) == 6
    assert (sys().coords == before).all()
    sys.prepare_for_ff(infer_hydrogens=True, add_hydrogens=True)
    assert len(sys) == 12


def test_connectivity_options_are_used(mol, benzene):
    sys = mol.System(benzene)
    opts = mol.PrepareOptions()
    opts.connectivity.tolerance = -0.01
    with pytest.raises(mol.BondPerceptionError) as caught:
        sys.prepare_for_ff(options=opts)
    assert caught.value.kind == "InvalidTolerance"
    assert caught.value.details["value"] == pytest.approx(-0.01)
    assert len(sys.topology) == 6


def test_options_cannot_silently_override_legacy_arguments(mol, benzene):
    sys = mol.System(benzene)
    with pytest.raises(ValueError, match="cannot be combined"):
        sys.prepare_for_ff(infer_hydrogens=True, options=mol.PrepareOptions())
    assert len(sys) == 6


def test_legacy_positional_preparation_is_unchanged(mol, benzene):
    sys = mol.System(benzene)
    sys.prepare_for_ff(True, True, 0)
    assert len(sys) == 12
    sys.apply_ff()


@pytest.mark.parametrize("method", ["apply_ff", "apply_charges"])
def test_invalid_choices_keep_value_error(mol, benzene, method):
    sys = mol.System(benzene)
    for obj in (sys, sys.topology, sys()):
        with pytest.raises(ValueError, match="unknown"):
            getattr(obj, method)("invalid")
        with pytest.raises(TypeError):
            getattr(obj, method)(object())


@pytest.mark.parametrize("formal_charge", [-1, 1])
@pytest.mark.parametrize("target", ["system", "topology", "selection"])
def test_charged_molecule_total_is_preserved(mol, tmp_path, formal_charge, target):
    # A single supported atom isolates charge equilibration from geometry.
    element = "Cl" if formal_charge < 0 else "H"
    path = tmp_path / "ion.sdf"
    path.write_text(
        "ion\n  pymolar\n\n  1  0  0  0  0  0            999 V2000\n"
        f"    0.0000    0.0000    0.0000 {element:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n"
        f"M  CHG  1   1 {formal_charge:3}\nM  END\n$$$$\n"
    )
    sys = mol.System(str(path))
    obj = {"system": sys, "topology": sys.topology, "selection": sys()}[target]
    obj.apply_charges(mol.ChargeModel.Espaloma)
    assert sum(a.charge for a in sys.iter_atoms()) == pytest.approx(formal_charge, abs=1e-6)
    # A second call also reads preserved formal charges, not old partial charges.
    obj.apply_charges()
    assert sys[0].charge == pytest.approx(formal_charge, abs=1e-6)


def test_unsupported_element_has_charge_error_details(mol, tmp_path):
    path = tmp_path / "sodium.xyz"
    path.write_text("1\nsodium\nNa 0 0 0\n")
    sys = mol.System(str(path))
    sys[0].atomic_number = 11  # isolate charge-model validation from XYZ name inference
    with pytest.raises(mol.ChargeError) as caught:
        sys.topology.apply_charges()
    assert caught.value.kind == "UnsupportedElement"
    assert caught.value.details == {"atomic_number": 11, "model": "espaloma"}



def test_nested_search_limit_is_used(mol, benzene):
    sys = mol.System(benzene)
    opts = mol.PrepareOptions()
    opts.bond_orders.hydrogens = mol.HydrogenPolicy.InferFromGeometry
    opts.bond_orders.use_functional_groups = False
    opts.bond_orders.use_residue_templates = False
    opts.bond_orders.limits.max_branches = 0
    with pytest.raises(mol.BondPerceptionError) as caught:
        sys.prepare_for_ff(options=opts)
    assert caught.value.kind == "SearchLimitExceeded"
    assert "atom" in caught.value.details


@pytest.mark.parametrize("zero_fill", [False, True])
def test_hydrogen_dynamics_option_is_used(mol, tmp_path, zero_fill):
    path = tmp_path / "carbon.gro"
    path.write_text(
        "carbon with velocity\n1\n"
        f"{1:5d}{'LIG':<5}{'C':>5}{1:5d}{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}"
        f"{0.1:8.4f}{0.2:8.4f}{0.3:8.4f}\n2.0 2.0 2.0\n"
    )
    sys = mol.System(str(path))
    opts = mol.PrepareOptions(
        bond_orders=mol.BondOrderOptions(hydrogens=mol.HydrogenPolicy.InferFromGeometry),
        add_hydrogens=mol.HydrogenOptions(zero_fill_dynamics=zero_fill),
    )
    if zero_fill:
        sys.prepare_for_ff(options=opts)
        assert len(sys) == 5
    else:
        with pytest.raises(mol.BondPerceptionError) as caught:
            sys.prepare_for_ff(options=opts)
        assert caught.value.kind == "DynamicsPresent"
        assert len(sys) == 1
