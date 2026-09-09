"""Distance search validation at the Python interface."""

import pytest
from pymolar import System, distance_search


@pytest.fixture
def selection(tmp_path):
    # XYZ has no periodic box.
    path = tmp_path / "points.xyz"
    path.write_text("2\nsearch test\nC 0 0 0\nC 2.5 0 0\n")
    system = System(str(path))
    return system("all")


@pytest.mark.parametrize("cutoff", [0.0, -0.5, float("nan"), float("inf")])
@pytest.mark.parametrize("double", [False, True])
def test_invalid_cutoff_is_value_error(selection, cutoff, double):
    with pytest.raises(ValueError, match="cutoff"):
        distance_search(cutoff, selection, selection if double else None)


@pytest.mark.parametrize("cutoff,double", [(0.5, False), (0.5, True), ("vdw", True)])
def test_missing_box_is_value_error(selection, cutoff, double):
    with pytest.raises(ValueError):
        distance_search(cutoff, selection, selection if double else None,
                        dims=[True, True, True])


def test_valid_search(selection):
    pairs, distances = distance_search(0.5, selection)
    assert pairs.shape == (1, 2)
    assert sorted(pairs[0].tolist()) == [0, 1]
    assert distances[0] == pytest.approx(0.25)
