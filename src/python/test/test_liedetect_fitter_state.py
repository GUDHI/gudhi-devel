"""Regression tests for OrbitFitter recomputation and state invalidation."""

import numpy as np
import pytest

from gudhi.datasets.linear_orbits import sample_from_lie_group_rep
from gudhi.liedetect import OrbitFitter
from gudhi.liedetect import orbit_fitter


def fit_torus(fitter):
    return fitter.closest_algebra(
        "torus", group_dim=1, frequency_max=2, method="abelian"
    )


@pytest.fixture
def fitted_orbit():
    points = sample_from_lie_group_rep("torus", ((1, 2),), 40)
    fitter = OrbitFitter(points)
    fitter.lie_pca(nb_neighbors=5, orbit_dim=1)
    fit_torus(fitter)
    fitter.sample_orbit(24)
    return fitter


def test_fitter_copies_input_points():
    points = np.eye(3)
    fitter = OrbitFitter(points)
    points[:] = 0
    np.testing.assert_array_equal(fitter.points, np.eye(3))


def test_successful_recomputations_clear_downstream_state(fitted_orbit):
    fitter = fitted_orbit
    operator = fitter.lie_pca_operator_
    fit_torus(fitter)
    assert fitter.lie_pca_operator_ is operator
    assert fitter.is_closest_algebra
    assert not fitter.is_sample_orbit
    assert fitter.orbit_ is fitter.hausdorff_distances_ is None

    fitter.sample_orbit(24)
    fitter.lie_pca(nb_neighbors=6, orbit_dim=1)
    assert fitter.is_lie_pca
    assert not fitter.is_closest_algebra
    assert fitter.group is fitter.representation_type_ is fitter.algebra_ is None
    assert not fitter.is_sample_orbit
    assert fitter.orbit_ is fitter.hausdorff_distances_ is None
    with pytest.raises(RuntimeError, match="before sampling"):
        fitter.sample_orbit(24)


def test_failed_recomputations_retain_last_success(fitted_orbit, monkeypatch):
    fitter = fitted_orbit
    operator = fitter.lie_pca_operator_
    representation = fitter.representation_type_
    algebra = fitter.algebra_
    orbit = fitter.orbit_
    distances = fitter.hausdorff_distances_

    with pytest.raises(ValueError, match="nb_neighbors"):
        fitter.lie_pca(nb_neighbors=0, orbit_dim=2)
    with pytest.raises(ValueError, match="ambient dimension"):
        fitter.closest_algebra("SO(3)", reps_to_test=[(3,)])

    def fail_distance(*args, **kwargs):
        raise ValueError("distance failed")

    monkeypatch.setattr(orbit_fitter, "print_hausdorff_distance", fail_distance)
    with pytest.raises(ValueError, match="distance failed"):
        fitter.sample_orbit(16)

    assert fitter.is_lie_pca and fitter.is_closest_algebra and fitter.is_sample_orbit
    assert fitter.orbit_dim == 1
    assert fitter.lie_pca_operator_ is operator
    assert fitter.group == "torus"
    assert fitter.representation_type_ is representation
    assert fitter.algebra_ is algebra
    assert fitter.orbit_ is orbit
    assert fitter.hausdorff_distances_ is distances
