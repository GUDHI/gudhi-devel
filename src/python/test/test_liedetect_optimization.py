"""Regression tests for stationary and reproducible LieDetect optimization."""

import numpy as np
import pytest

from gudhi.datasets.linear_orbits import sample_from_lie_group_rep
from gudhi.liedetect import OrbitFitter
from gudhi.liedetect import optimization


def test_default_fit_on_planar_circle():
    points = sample_from_lie_group_rep("torus", ((1,),), 100)
    fitter = OrbitFitter(points)
    fitter.lie_pca(10, 1)
    representation, algebra = fitter.closest_algebra("torus", group_dim=1, frequency_max=2)
    assert representation == ((1,),)
    np.testing.assert_allclose(algebra[0] + algebra[0].T, 0, atol=1e-12)
    assert fitter.fit_diagnostics_["converged"]
    assert fitter.fit_diagnostics_["objective_value"] == pytest.approx(0, abs=1e-12)
    assert len(fitter.fit_diagnostics_["runs"]) == 2
    fitter.sample_orbit(100)
    assert fitter.hausdorff_distances_[0] < 1e-10


@pytest.mark.parametrize("method", ["bottom_lie_pca", "full_lie_pca", "abelian"])
def test_seeded_fit_and_sampling(method):
    points = sample_from_lie_group_rep("torus", ((1, 0, 1), (0, 1, 1)), 90)
    fitter = OrbitFitter(points)
    fitter.lie_pca(12, 2)
    options = {"max_iterations": 2}
    candidates = [((1, 0, 1), (0, 1, 1))]
    first_rep, first_algebra = fitter.closest_algebra(
        "torus", reps_to_test=candidates, method=method, seed=17, optimizer_options=options
    )
    first_cost = fitter.fit_diagnostics_["objective_value"]
    assert np.isfinite(first_cost)
    assert 1 <= fitter.fit_diagnostics_["iterations"] <= 2

    second_rep, second_algebra = fitter.closest_algebra(
        "torus", reps_to_test=candidates, method=method, seed=17, optimizer_options=options
    )
    assert second_rep == first_rep
    np.testing.assert_allclose(second_algebra, first_algebra, atol=1e-12)
    assert fitter.fit_diagnostics_["objective_value"] == pytest.approx(first_cost)
    first = fitter.sample_orbit(20, method="random_uniform", seed=1)
    np.testing.assert_array_equal(first, fitter.sample_orbit(20, method="random_uniform", seed=1))
    assert not np.allclose(first, fitter.sample_orbit(20, method="random_uniform", seed=2))


def test_initial_point_and_nonfinite_runs(monkeypatch):
    fitter = OrbitFitter(sample_from_lie_group_rep("torus", ((1,),), 24))
    fitter.lie_pca(5, 1)
    options = {"group_dim": 1, "frequency_max": 2}
    with pytest.raises(ValueError, match="initial_point"):
        fitter.closest_algebra("torus", initial_point=np.diag([-1, 1]), **options)

    original = optimization.pymanopt.optimizers.SteepestDescent.run
    reject_all = False
    starts = []

    def run(solver, problem, *, initial_point):
        starts.append(initial_point.copy())
        result = original(solver, problem, initial_point=initial_point)
        if reject_all or len(starts) == 1:
            result.cost = np.nan
        return result

    monkeypatch.setattr(optimization.pymanopt.optimizers.SteepestDescent, "run", run)
    fitter.closest_algebra("torus", initial_point=np.eye(2), **options)
    np.testing.assert_array_equal(starts, [np.eye(2), np.eye(2)])
    assert fitter.fit_diagnostics_["runs"][(((1,),), "+1")]["error"]
    diagnostics = fitter.fit_diagnostics_
    reject_all = True
    with pytest.raises(RuntimeError, match="all optimizer runs were nonfinite"):
        fitter.closest_algebra("torus", initial_point=np.eye(2), **options)
    assert fitter.fit_diagnostics_ is diagnostics


def test_su2_fitting_returns_generators_with_cyclic_brackets():
    fitter = OrbitFitter(sample_from_lie_group_rep("SU(2)", (4,), 100))
    fitter.lie_pca(14, 3)
    representation, generators = fitter.closest_algebra(
        "SU(2)", method="full_lie_pca", seed=9
    )
    assert representation == (4,)
    a, b, c = generators
    for first, second, bracket in ((a, b, c), (b, c, a), (c, a, b)):
        np.testing.assert_allclose(first + first.T, 0, atol=1e-12)
        np.testing.assert_allclose(first @ second - second @ first, bracket, atol=1e-12)
    casimir = sum(generator @ generator for generator in generators)
    np.testing.assert_allclose(casimir, np.trace(casimir) / 4 * np.eye(4), atol=1e-12)
