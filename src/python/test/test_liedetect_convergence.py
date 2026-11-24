"""Regression tests for stationarity and start selection in LieDetect."""

from types import SimpleNamespace

import numpy as np
import pytest

from gudhi.datasets.linear_orbits import sample_from_lie_group_rep
from gudhi.liedetect import OrbitFitter, optimization


def _rotation(angle):
    cosine, sine = np.cos(angle), np.sin(angle)
    return np.array([[cosine, -sine], [sine, cosine]])


def _rotation_problem(*, nonfinite_gradient=False):
    manifold = optimization.pymanopt.manifolds.SpecialOrthogonalGroup(2)

    @optimization.pymanopt.function.autograd(manifold)
    def cost(point):
        return -point[0, 0]

    if nonfinite_gradient:
        @optimization.pymanopt.function.autograd(manifold)
        def gradient(point):
            return np.full(point.shape, np.nan)

        return optimization.pymanopt.Problem(manifold, cost, riemannian_gradient=gradient)
    return optimization.pymanopt.Problem(manifold, cost)


def test_returned_gradient_norm_matches_analytic_gradient(monkeypatch):
    original = optimization.pymanopt.optimizers.SteepestDescent.run
    reported = []

    def run(solver, *args, **kwargs):
        result = original(solver, *args, **kwargs)
        reported.append(result.gradient_norm)
        return result

    monkeypatch.setattr(optimization.pymanopt.optimizers.SteepestDescent, "run", run)
    angle = np.pi / 3
    result = optimization._run_optimizer(
        _rotation_problem(), 2, None, _rotation(angle), {"max_iterations": 1}
    )
    # The angle curve has Frobenius speed sqrt(2), and the objective is -cos(angle).
    assert reported[0] == pytest.approx(np.sin(angle) / np.sqrt(2))
    assert result.gradient_norm == pytest.approx(abs(result.point[1, 0]) / np.sqrt(2))
    assert abs(result.gradient_norm - reported[0]) > 0.1
    assert result.cost < -np.cos(angle)


def test_point_inside_gradient_tolerance_is_left_in_place():
    initial_point = _rotation(1e-3)
    result = optimization._run_optimizer(
        _rotation_problem(), 2, None, initial_point, {"min_gradient_norm": 1e-2}
    )
    np.testing.assert_array_equal(result.point, initial_point)
    assert 0 < result.gradient_norm < 1e-2
    assert result.iterations == 1


@pytest.mark.parametrize("failure", ["off_manifold", "nonfinite_gradient"])
def test_returned_estimate_is_validated(monkeypatch, failure):
    nonfinite_gradient = failure == "nonfinite_gradient"
    problem = _rotation_problem(nonfinite_gradient=nonfinite_gradient)
    point = -np.eye(2) if nonfinite_gradient else 2 * np.eye(2)

    def run(*args, **kwargs):
        return SimpleNamespace(point=point, cost=1.0, gradient_norm=0)

    monkeypatch.setattr(optimization.pymanopt.optimizers.SteepestDescent, "run", run)
    error = optimization._NonfiniteOptimizationError if nonfinite_gradient else optimization._InvalidOptimizationError
    with pytest.raises(error):
        optimization._run_optimizer(problem, 2, None, np.eye(2), None)


def _result(point, cost, gradient_norm):
    return SimpleNamespace(
        point=point, cost=cost, gradient_norm=gradient_norm,
        iterations=1, stopping_criterion="Iteration limit",
    )


@pytest.mark.parametrize("method", ["bottom_lie_pca", "full_lie_pca"])
def test_multistart_selection_and_ties(monkeypatch, method):
    reps = [((1, 2),), ((1, 3),)]
    fitter = OrbitFitter(sample_from_lie_group_rep("torus", reps[0], 32))
    fitter.lie_pca(6, 1)
    initial = np.eye(4)
    scripted = [
        (0.8, 0), (0.05, 0.1), (0.25, 0),
        (0.2, 0), (0.2, 0), (0.3, 0),
        (0.2, 0), (0.2, 0), (0.2, 0),
        (0.3, 0), (0.3, 0), (0.3, 0),
    ]

    def run(require_convergence, outcomes):
        options = iter(outcomes)
        received = []

        def solve(problem, ambient_dim, rng, initial_point, optimizer_options):
            received.append(initial_point)
            cost, gradient = next(options)
            return _result(initial, cost, gradient)

        monkeypatch.setattr(optimization, "_run_optimizer", solve)
        winner, _ = fitter.closest_algebra(
            "torus", reps_to_test=reps, method=method, n_starts=3,
            initial_point=initial, seed=17, require_convergence=require_convergence,
        )
        assert len(received) == 12
        for index, point in enumerate(received):
            if index % 3 == 0:
                np.testing.assert_array_equal(point, initial)
            else:
                assert point is None
        return winner, fitter.fit_diagnostics_

    winner, permissive = run(False, scripted)
    assert winner == reps[0]
    assert permissive["objective_value"] == 0.05
    assert permissive["gradient_norm"] == 0.1
    assert not permissive["converged"]
    assert permissive["runs"][(reps[0], "+1")]["selected_start"] == 1

    winner, strict = run(True, scripted)
    assert winner == reps[0]  # Earlier representation and component win equal-cost ties.
    assert strict["objective_value"] == 0.2
    assert strict["converged"]  # The solver stop reason can still be an iteration limit.
    assert strict["runs"][(reps[0], "+1")]["starts"][1]["status"] == "not_stationary"
    assert strict["runs"][(reps[0], "-1")]["selected_start"] == 0

    with pytest.raises(RuntimeError, match="no eligible run reached the gradient tolerance"):
        run(True, [(cost, 0.1) for cost, _ in scripted])
    assert fitter.fit_diagnostics_ is strict


def test_multistart_seed_reproduces_random_starts(monkeypatch):
    original = optimization.pymanopt.optimizers.SteepestDescent.run
    received = []

    def run(solver, problem, *, initial_point):
        received.append(initial_point.copy())
        return original(solver, problem, initial_point=initial_point)

    monkeypatch.setattr(optimization.pymanopt.optimizers.SteepestDescent, "run", run)
    rep = ((1, 2),)
    fitter = OrbitFitter(sample_from_lie_group_rep("torus", rep, 32))
    fitter.lie_pca(6, 1)
    for _ in range(2):
        fitter.closest_algebra(
            "torus", reps_to_test=[rep], n_starts=3, initial_point=np.eye(4),
            seed=7, optimizer_options={"max_iterations": 1},
        )
    np.testing.assert_array_equal(received[:6], received[6:])
    np.testing.assert_array_equal(received[0], np.eye(4))
    np.testing.assert_array_equal(received[3], np.eye(4))
    assert not np.allclose(received[1], received[2])


def test_joint_normal_form_loss_and_strict_multistart(monkeypatch):
    j = np.array([[0., 1.], [-1., 0.]])
    zero = np.zeros_like(j)
    matrices = [np.block([[2 * j, zero], [zero, j]]), np.block([[j, zero], [zero, 3 * j]])]
    mixed = np.eye(4)
    mixed[np.ix_([0, 2], [0, 2])] = _rotation(0.4)
    points = iter((np.eye(4), mixed))
    block_mask = np.kron(np.eye(2, dtype=bool), np.ones((2, 2), dtype=bool))
    energies = [np.linalg.norm(matrix, "fro") ** 2 for matrix in matrices]

    def run(problem, ambient_dim, rng, initial_point, optimizer_options):
        point = next(points)
        mixed_blocks = [point @ matrix @ point.T for matrix in matrices]
        off_block_energy = sum(np.sum(matrix[~block_mask] ** 2) for matrix in mixed_blocks)
        cost = problem.cost(point)
        assert cost == pytest.approx(off_block_energy + sum(energy ** 2 - energy for energy in energies))
        # Simulate an ineligible lowest-cost result to isolate the start selection policy.
        return _result(point, cost, 0.1 if off_block_energy == 0 else 0)

    monkeypatch.setattr(optimization, "_run_optimizer", run)
    runs = {}
    _, basis = optimization.normal_form_skew_symmetric_matrices(
        matrices, n_starts=2, require_convergence=True, runs=runs
    )
    np.testing.assert_array_equal(basis, mixed.T)
    summary = runs["joint_normal_form"]
    assert summary["selected_start"] == 1
    assert summary["starts"][0]["cost"] < summary["starts"][1]["cost"]
    assert summary["starts"][0]["status"] == "not_stationary"


def test_rank_one_abelian_has_no_iterative_policy(monkeypatch):
    def fail(*args, **kwargs):
        pytest.fail("Rank-one Abelian fitting should use direct Schur decomposition.")

    monkeypatch.setattr(optimization, "_run_optimizer", fail)
    rep = ((1,),)
    fitter = OrbitFitter(sample_from_lie_group_rep("torus", rep, 24))
    fitter.lie_pca(5, 1)
    fitted, _ = fitter.closest_algebra(
        "torus", reps_to_test=[rep], method="abelian", require_convergence=True, n_starts=3
    )
    assert fitted == rep
    assert fitter.fit_diagnostics_["gradient_norm"] is None
    assert fitter.fit_diagnostics_["converged"] is None


@pytest.mark.parametrize("n_starts", [0, -1, 1.5, True])
def test_invalid_start_count(n_starts):
    with pytest.raises(ValueError, match="n_starts"):
        optimization.find_closest_algebra("torus", np.eye(4), n_starts=n_starts)
