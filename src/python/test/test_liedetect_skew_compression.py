"""Independent regressions for raw and skew-compressed Lie-PCA."""

from types import SimpleNamespace

import numpy as np
import pytest

from gudhi.liedetect import OrbitFitter, optimization, orbit_fitter
from gudhi.liedetect.liepca import _compress_lie_pca_to_skew, get_lie_pca_operator


def _matrix_space_projections(dimension):
    symmetric, skew = [], []
    for basis in np.eye(dimension**2).reshape(-1, dimension, dimension):
        symmetric.append(((basis + basis.T) / 2).ravel())
        skew.append(((basis - basis.T) / 2).ravel())
    return np.column_stack(symmetric), np.column_stack(skew)


@pytest.mark.parametrize("method", ["PCA", "covariance"])
def test_operator_defaults_to_raw_and_opt_in_compresses(method):
    points = np.random.default_rng(0).normal(size=(40, 4))
    points /= np.linalg.norm(points, axis=1, keepdims=True)
    symmetric, skew = _matrix_space_projections(4)
    raw = get_lie_pca_operator(points, 8, 2, method, compress_to_skew=False)
    assert np.linalg.norm(symmetric @ raw @ skew) > 1e-3

    np.testing.assert_array_equal(get_lie_pca_operator(points, 8, 2, method), raw)
    fitter = OrbitFitter(points)
    np.testing.assert_array_equal(fitter.lie_pca(8, 2, method), raw)
    compressed = fitter.lie_pca(8, 2, method, compress_to_skew=True)
    np.testing.assert_allclose(compressed, symmetric + skew @ raw @ skew, atol=1e-14)
    np.testing.assert_allclose(compressed @ symmetric, symmetric, atol=1e-14)
    np.testing.assert_allclose(symmetric @ compressed, symmetric, atol=1e-14)
    np.testing.assert_allclose(symmetric @ compressed @ skew, 0, atol=1e-14)
    np.testing.assert_array_equal(compressed, compressed.T)

    if method == "PCA":
        for operator in (raw, compressed):
            values = np.linalg.eigvalsh(operator)
            assert values[0] >= -1e-14
            assert values[-1] <= 1 + 1e-14


def test_old_correction_can_lose_positivity_for_a_noncontractive_operator():
    symmetric, skew = _matrix_space_projections(2)
    sym_direction = (np.eye(2) / np.sqrt(2)).ravel()
    skew_direction = np.array([[0, 1], [-1, 0]]).ravel() / np.sqrt(2)
    vector = 3 * sym_direction + skew_direction
    raw = np.outer(vector, vector)
    old = symmetric + (raw @ skew + skew @ raw) / 2
    assert np.linalg.eigvalsh(old)[0] == pytest.approx(-0.5)
    assert np.linalg.eigvalsh(_compress_lie_pca_to_skew(raw, 2))[0] == pytest.approx(1)


def test_covariance_estimator_can_be_indefinite_on_skew_inputs():
    offsets = np.array([0.0, 0.01, 0.03, 0.08])
    points = np.column_stack((np.ones(4), offsets, np.zeros(4)))
    direction = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]]).ravel()
    expected = np.mean((offsets**2 + 1 - np.sqrt(2)) / (1 + offsets**2))
    assert expected < 0
    for compress in (False, True):
        operator = get_lie_pca_operator(points, 2, 2, "covariance", compress)
        assert direction @ operator @ direction == pytest.approx(expected)


def test_ideal_sphere_kernel_and_perturbation():
    # For a uniform unit sphere, E[x_i x_j x_k x_l] gives this raw operator.
    n = 3
    columns = [
        ((np.trace(basis) * np.eye(n) + basis + basis.T) / (n * (n + 2))).ravel()
        for basis in np.eye(n**2).reshape(-1, n, n)
    ]
    ideal = np.column_stack(columns)
    corrected_ideal = _compress_lie_pca_to_skew(ideal, n)
    symmetric, skew = _matrix_space_projections(n)
    np.testing.assert_allclose(corrected_ideal, symmetric, atol=1e-14)
    np.testing.assert_allclose(corrected_ideal @ skew, 0, atol=1e-14)
    assert np.linalg.matrix_rank(ideal) == np.linalg.matrix_rank(corrected_ideal) == 6

    matrix = np.random.default_rng(3).normal(size=(n**2, n**2))
    empirical = matrix.T @ matrix
    assert np.linalg.norm(symmetric @ empirical @ skew) > 1
    raw_difference = empirical - ideal
    corrected_difference = _compress_lie_pca_to_skew(empirical, n) - corrected_ideal
    np.testing.assert_allclose(corrected_difference, skew @ raw_difference @ skew)
    for norm in (2, "fro"):
        assert np.linalg.norm(corrected_difference, norm) <= np.linalg.norm(raw_difference, norm)


def test_full_fit_uses_the_selected_operator(monkeypatch):
    symmetric, skew = _matrix_space_projections(2)
    # L(A)=E_11 A E_22 retains only A_12. Its raw cross block is nonzero.
    raw = np.diag([0.0, 1.0, 0.0, 0.0])
    legacy = symmetric + (raw @ skew + skew @ raw) / 2

    def assemble(*, compress_to_skew, **kwargs):
        return _compress_lie_pca_to_skew(raw, 2) if compress_to_skew else raw

    def evaluate(problem, ambient_dim, rng, initial_point, optimizer_options):
        return SimpleNamespace(
            point=initial_point,
            cost=problem.cost(initial_point),
            iterations=0,
            stopping_criterion="Fixed test conjugation",
            gradient_norm=0.0,
        )

    monkeypatch.setattr(orbit_fitter, "get_lie_pca_operator", assemble)
    monkeypatch.setattr(optimization, "_run_optimizer", evaluate)
    fitter = OrbitFitter(np.eye(2))
    options = dict(reps_to_test=[((1,),)], method="full_lie_pca", initial_point=np.eye(2))
    fitter.lie_pca(1, 1)
    fitter.closest_algebra("torus", **options)
    assert fitter.fit_diagnostics_["objective_value"] == pytest.approx(1 / 2)

    old_result = optimization.optimization_full_lie_pca(
        "torus", legacy, ((1,),), initial_point=np.eye(2)
    )
    assert old_result.cost == pytest.approx(5 / 16)

    fitter.lie_pca(1, 1, compress_to_skew=True)
    fitter.closest_algebra("torus", **options)
    assert fitter.fit_diagnostics_["objective_value"] == pytest.approx(1 / 4)


@pytest.mark.parametrize("tie", [False, True])
@pytest.mark.parametrize("method", ["bottom_lie_pca", "abelian"])
def test_bottom_methods_can_lose_skew_rank_at_a_tie(tie, method, monkeypatch):
    symmetric, skew = _matrix_space_projections(2)
    operator = symmetric + (1 if tie else 0.5) * skew
    sym_direction = np.array([1.0, 0.0, 0.0, 0.0])
    skew_direction = np.array([0.0, 1.0, -1.0, 0.0]) / np.sqrt(2)
    other_sym = np.array([0.0, 1.0, 1.0, 0.0]) / np.sqrt(2)
    last_sym = np.array([0.0, 0.0, 0.0, 1.0])
    vectors = np.column_stack(
        (sym_direction, skew_direction, other_sym, last_sym)
        if tie else (skew_direction, sym_direction, other_sym, last_sym)
    )
    eigenvalues = np.array([1.0 if tie else 0.5, 1.0, 1.0, 1.0])
    np.testing.assert_allclose(operator @ vectors, vectors * eigenvalues)

    class FrameObserved(Exception):
        pass

    def inspect_frame(frame):
        assert np.linalg.norm(frame[0]) == pytest.approx(0.0 if tie else 1.0)
        raise FrameObserved

    monkeypatch.setattr(optimization.np.linalg, "eigh", lambda matrix: (eigenvalues, vectors))
    monkeypatch.setattr(optimization, "skew_sym_frame_to_projection", inspect_frame)
    with pytest.raises(FrameObserved):
        optimization.find_closest_algebra(
            "torus", operator, reps_to_test=[((1,),)], method=method, verbose=False
        )
