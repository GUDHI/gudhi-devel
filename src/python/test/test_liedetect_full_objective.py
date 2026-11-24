"""Independent mathematical regressions for the full Lie-PCA objective."""

from types import SimpleNamespace

import numpy as np
import pytest

from gudhi.liedetect import optimization
from gudhi.liedetect.utils import sample_from_lie_algebra


def _rotation_generator(weights):
    generator = np.zeros((2 * len(weights), 2 * len(weights)))
    for i, weight in enumerate(weights):
        generator[2 * i, 2 * i + 1] = weight
        generator[2 * i + 1, 2 * i] = -weight
    return generator


def test_normalized_full_objective_reverses_the_raw_ranking():
    aligned = _rotation_generator((1, 3)).ravel() / np.sqrt(20)
    operator = np.eye(16) - 0.8 * np.outer(aligned, aligned)
    candidates = [((1, 2),), ((1, 3),)]
    diagnostics = {}

    representation, algebra = optimization.find_closest_algebra(
        "torus",
        operator,
        reps_to_test=candidates,
        method="full_lie_pca",
        initial_point=np.eye(4),
        verbose=False,
        diagnostics=diagnostics,
    )

    assert diagnostics["runs"][(candidates[0], "+1")]["cost"] == pytest.approx(0.0592)
    assert diagnostics["runs"][(candidates[1], "+1")]["cost"] == pytest.approx(0.04)
    assert diagnostics["objective_value"] == pytest.approx(0.04)
    assert representation == ((1, 3),)
    assert np.linalg.norm(algebra[0]) ** 2 == pytest.approx(20)

    base = np.array([1.0, 0.0, 1.0, 0.0])
    points = sample_from_lie_algebra("torus", representation, algebra, 24, x=base)
    first_plane = points[:, 0] + 1j * points[:, 1]
    second_plane = points[:, 2] + 1j * points[:, 3]
    np.testing.assert_allclose(second_plane, first_plane**3, rtol=0, atol=1e-10)


@pytest.mark.parametrize("determinant", ["+1", "-1"])
def test_full_objective_is_invariant_under_nonorthogonal_rebasing(determinant, monkeypatch):
    weights = ((1, 0, 1), (0, 1, 1))
    rebased = ((1, 2, 3), (0, 1, 1))
    angle = 0.3
    rotation = np.eye(6)
    rotation[np.ix_([0, 2], [0, 2])] = [
        [np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]
    ]
    conjugation = rotation.copy()
    if determinant == "-1":
        conjugation[:, 0] *= -1

    columns = np.column_stack([
        (conjugation @ _rotation_generator(row) @ conjugation.T).ravel()
        for row in weights
    ])
    projector = columns @ np.linalg.solve(columns.T @ columns, columns.T)
    operator = np.random.default_rng(5).normal(size=(36, 36))
    expected = np.trace(operator.T @ operator @ projector)

    def evaluate_at_initial_point(problem, ambient_dim, rng, initial_point, optimizer_options):
        # Inspect the exact objective at a controlled conjugation, without solver motion.
        return SimpleNamespace(point=initial_point, cost=problem.cost(initial_point))

    monkeypatch.setattr(optimization, "_run_optimizer", evaluate_at_initial_point)
    for rows in (weights, rebased):
        result = optimization.optimization_full_lie_pca(
            "torus", operator, rows, determinant, initial_point=rotation
        )
        np.testing.assert_allclose(result.point, conjugation)
        assert result.cost == pytest.approx(expected, rel=1e-12)


def test_frobenius_frame_remains_orthonormal_after_conjugation():
    frame = optimization._frobenius_orthonormal_frame([
        _rotation_generator((1, 0, 1)),
        _rotation_generator((1, 2, 3)),
    ])
    rotation = np.eye(6)
    rotation[:2, :2] = [[0.6, -0.8], [0.8, 0.6]]
    reflection = np.diag([-1, 1, 1, 1, 1, 1])
    for orthogonal in (rotation, rotation @ reflection):
        columns = (orthogonal @ frame @ orthogonal.T).reshape(2, -1)
        np.testing.assert_allclose(columns @ columns.T, np.eye(2), rtol=0, atol=1e-12)


def test_full_objective_rejects_dependent_generators():
    with pytest.raises(ValueError, match="linearly independent candidate generators"):
        optimization.optimization_full_lie_pca(
            "torus", np.eye(36), ((1, 0, 1), (2, 0, 2))
        )
