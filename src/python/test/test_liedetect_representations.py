"""Regression tests for LieDetect lattice sampling and representation validation."""

import numpy as np
import pytest

from gudhi.datasets.linear_orbits import sample_from_lie_group, sample_from_lie_group_rep
from gudhi.liedetect import OrbitFitter
from gudhi.liedetect.algebra import get_random_lattice


def test_random_torus_allows_repeated_weights():
    def sample():
        return sample_from_lie_group("torus", 8, 10, frequency_max=1, group_dim=1, seed=1)

    points, weights = sample()
    repeated_points, repeated_weights = sample()
    assert points.shape == (10, 8)
    assert len(set(weights[0])) < len(weights[0])
    assert weights == repeated_weights
    np.testing.assert_array_equal(points, repeated_points)


def test_random_torus_validates_dimensions():
    with pytest.raises(ValueError, match="distinct nonzero frequency"):
        get_random_lattice(1, 3, 2, span_ambient_space=True, seed=1)
    with pytest.raises(ValueError, match="positive even integer"):
        sample_from_lie_group("torus", 7, 10, frequency_max=1, group_dim=1)
    with pytest.raises(ValueError, match="lattice_rank.*positive integer"):
        get_random_lattice(0, 4, 1)
    (weights,) = get_random_lattice(1, 12, 12, span_ambient_space=True, seed=1)
    assert sorted(map(abs, weights)) == list(range(1, 13))


@pytest.mark.parametrize(
    ("group", "rep_type", "message"),
    [
        ("SO(3)", (4,), "Invalid real irrep dimensions"),
        ("SU(2)", (2,), "Invalid real irrep dimensions"),
        ("SU(2)", (), "Invalid real irrep dimensions"),
        ("torus", ((1, 2), (1,)), "rectangular array of integer weights"),
        ("torus", ((1.5, 2),), "rectangular array of integer weights"),
    ],
)
def test_invalid_explicit_representation(group, rep_type, message):
    with pytest.raises(ValueError, match=message):
        sample_from_lie_group_rep(group, rep_type, 8)


def test_valid_explicit_representation():
    assert sample_from_lie_group_rep("SU(2)", (4, 1), 8).shape == (8, 5)
    assert sample_from_lie_group_rep("torus", ((2, 4),), 8).shape == (8, 4)
    with pytest.raises(ValueError, match="finite real vector of length 4"):
        sample_from_lie_group_rep("torus", ((1, 2),), 8, x=np.ones(3))


def test_explicit_candidates_are_validated_before_fitting():
    points = sample_from_lie_group_rep("torus", ((1, 2),), 20)
    fitter = OrbitFitter(points)
    fitter.lie_pca(nb_neighbors=4, orbit_dim=1)
    with pytest.raises(ValueError, match="does not match ambient dimension"):
        fitter.closest_algebra("torus", reps_to_test=[((1,),)])
    with pytest.raises(ValueError, match="must have row rank 2"):
        fitter.closest_algebra("torus", reps_to_test=[((1, 2), (2, 4))])
    with pytest.raises(ValueError, match="only supported for a torus"):
        fitter.closest_algebra("SO(3)", method="abelian")
