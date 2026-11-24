"""Regressions for the final LieDetect review findings."""

from pathlib import Path
from textwrap import dedent

import numpy as np
import pytest

from gudhi.datasets.linear_orbits import sample_from_lie_group_rep
from gudhi.liedetect import OrbitFitter, optimization, orbit_fitter
from gudhi.liedetect.algebra import get_canonical_pushforward_algebra


@pytest.mark.parametrize("candidate", [((2, 1),), ((-1, 2),), ((2, -1),)])
def test_rank_one_abelian_reconstructs_the_scored_orbit(candidate, monkeypatch):
    canonical = ((1, 2),)
    weights = np.array(canonical[0])
    angles = 2 * np.pi * np.arange(80)[:, None] / 80
    phases = angles * weights
    points = np.column_stack([
        np.cos(phases[:, 0]), -np.sin(phases[:, 0]),
        np.cos(phases[:, 1]), -np.sin(phases[:, 1]),
    ]) / 2
    generator = np.kron(np.diag(weights), [[0, 1], [-1, 0]])
    unit_generator = generator.ravel() / np.linalg.norm(generator)
    operator = np.eye(16) - np.outer(unit_generator, unit_generator)
    monkeypatch.setattr(orbit_fitter, "get_lie_pca_operator", lambda **kwargs: operator)

    fitter = OrbitFitter(points)
    fitter.lie_pca(nb_neighbors=5, orbit_dim=1)
    fitted, algebra = fitter.closest_algebra(
        "torus", method="abelian", reps_to_test=[candidate, canonical]
    )
    assert fitted == candidate
    assert fitter.fit_diagnostics_["objective_value"] == pytest.approx(0, abs=1e-14)
    fitted_generator = algebra[0].ravel() / np.linalg.norm(algebra[0])
    np.testing.assert_allclose(operator @ fitted_generator, 0, atol=1e-14)
    fitter.sample_orbit(len(points), x=points[0])
    assert fitter.hausdorff_distances_[0] < 1e-12


@pytest.mark.parametrize("method", ["bottom_lie_pca", "abelian"])
@pytest.mark.parametrize("group_dim", [1, 2])
def test_bottom_methods_reject_rank_loss_after_skew_projection(method, group_dim):
    n = 2 * group_dim
    symmetric = np.zeros((n, n))
    symmetric[0, 0] = 1
    kernel_frame = [symmetric.ravel()]
    if group_dim == 2:
        skew = np.zeros((n, n))
        skew[0, 1] = 1 / np.sqrt(2)
        skew[1, 0] = -skew[0, 1]
        kernel_frame.append(skew.ravel())
    frame = np.array(kernel_frame)
    operator = np.eye(n**2) - frame.T @ frame
    weights = tuple(tuple(int(i == j) for j in range(group_dim)) for i in range(group_dim))
    with pytest.raises(
        ValueError, match=f"numerical rank {group_dim - 1} after skew projection, expected {group_dim}"
    ):
        optimization.find_closest_algebra(
            "torus", operator, reps_to_test=[weights], method=method, verbose=False
        )


def test_random_torus_ignores_zero_rows_and_keeps_axis_draw_order():
    count, seed = 10, 4
    weight_table = ((0, 0), (2, 4), (0, 0))
    base = np.array([1.0, 0.0, 1.0, 0.0])
    samples = sample_from_lie_group_rep(
        "torus", weight_table, count, x=base, method="random_uniform", seed=seed
    )
    angle = 2 * np.pi * np.random.default_rng(seed).random((count, 3))[:, 1]
    expected = np.column_stack([np.cos(angle), -np.sin(angle), np.cos(2 * angle), -np.sin(2 * angle)])
    np.testing.assert_allclose(samples, expected, atol=1e-12)
    constant = sample_from_lie_group_rep(
        "torus", ((0, 0),), count, x=base, method="random_uniform", seed=seed
    )
    np.testing.assert_array_equal(constant, np.tile(base, (count, 1)))


@pytest.mark.parametrize("group", ["SO(3)", "SU(2)"])
@pytest.mark.parametrize("method", ["bottom_lie_pca", "full_lie_pca"])
def test_explicit_nonabelian_fitting_requires_a_nontrivial_irrep(group, method):
    fitter = OrbitFitter(np.eye(3))
    fitter.lie_pca(nb_neighbors=2, orbit_dim=1)
    with pytest.raises(ValueError, match="must contain a nontrivial irrep"):
        fitter.closest_algebra(group, reps_to_test=[(3,), (1, 1, 1)], method=method)


def test_explicit_torus_rank_error_identifies_numerical_limit():
    a = 10**8
    weights = ((a, a + 1), (a + 1, a + 2))
    assert a * (a + 2) - (a + 1)**2 == -1
    with pytest.raises(ValueError, match="row rank 2 at working precision"):
        optimization.find_closest_algebra(
            "torus", np.eye(16), reps_to_test=[weights], method="full_lie_pca", verbose=False
        )


def test_linear_orbits_manual_examples_run_independently():
    manual = Path(__file__).resolve().parents[1] / "doc" / "linear_orbits.rst"
    blocks = manual.read_text(encoding="utf-8").split(".. code-block:: python\n\n")[1:]
    assert len(blocks) == 3
    for block in blocks:
        scope = {}
        code = dedent(block.split("\n.. autofunction::", 1)[0])
        exec(compile(code, str(manual), "exec"), scope)
        assert scope["gen_points"].shape == (500, 4)
        assert np.isfinite(scope["gen_points"]).all()


@pytest.mark.parametrize(
    ("group", "dimension"),
    [
        ("SU(2)", d) for d in (1, 3, 4, 5, 7, 8, 9, 11, 12, 15, 16, 20, 24, 32)
    ] + [
        ("SO(3)", d) for d in (1, 3, 5, 7, 9, 11, 15)
    ],
)
def test_canonical_irrep_brackets_and_casimir(group, dimension):
    a, b, c = get_canonical_pushforward_algebra(group, (dimension,))
    for generator in (a, b, c):
        np.testing.assert_allclose(generator + generator.T, 0, atol=1e-12)
    for left, right, bracket in ((a, b, c), (b, c, a), (c, a, b)):
        np.testing.assert_allclose(left @ right - right @ left, bracket, rtol=0, atol=1e-12)
    spin = (dimension - 1) / 2 if dimension % 2 else (dimension - 2) / 4
    expected = -spin * (spin + 1) * np.eye(dimension)
    np.testing.assert_allclose(a @ a + b @ b + c @ c, expected, rtol=0, atol=1e-12)
