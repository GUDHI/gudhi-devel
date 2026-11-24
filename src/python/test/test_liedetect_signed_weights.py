"""Regression tests for complete signed torus weight searches."""

from itertools import product

import numpy as np
import pytest

from gudhi.liedetect import OrbitFitter, are_representations_equivalent
from gudhi.liedetect.algebra import get_lattices


def test_mixed_sign_four_block_orbit_is_in_the_bound_one_search():
    weights = ((1, 0, 1, 1), (0, 1, 1, -1))
    angles = np.array(list(product(np.arange(4) * np.pi / 2, repeat=2)))
    phases = angles @ np.array(weights)
    grid = np.stack((np.cos(phases), np.sin(phases)), axis=-1).reshape(-1, 8)
    np.testing.assert_allclose(grid.mean(axis=0), 0, atol=1e-15)
    np.testing.assert_allclose(grid.T @ grid / len(grid), np.eye(8) / 2, atol=1e-15)

    candidates = get_lattices(2, 4, 1, "orbit-equivalence")
    assert len(candidates) == 1
    assert are_representations_equivalent("torus", candidates[0], weights)

    # Independent random angle samples avoid the distance ties of the proof grid.
    angles = np.random.default_rng(17).uniform(0, 2 * np.pi, (240, 2))
    phases = angles @ np.array(weights)
    cloud = np.stack((np.cos(phases), np.sin(phases)), axis=-1).reshape(-1, 8)
    fitter = OrbitFitter(cloud)
    fitter.lie_pca(15, 2)
    inferred, _ = fitter.closest_algebra("torus", group_dim=2, frequency_max=1, seed=10)
    assert are_representations_equivalent("torus", inferred, weights)


def _character_classes(weights):
    return {max(w, tuple(-entry for entry in w)) for w in zip(*weights)}


@pytest.mark.parametrize("method", ["span-equivalence", "orbit-equivalence"])
def test_signed_search_contracts_for_distinct_and_repeated_characters(method):
    full_span = get_lattices(2, 3, 1, method)
    unrestricted = get_lattices(2, 3, 1, method, span_ambient_space=False)
    assert full_span and unrestricted
    assert all(len(_character_classes(weights)) == 3 for weights in full_span)
    assert any(len(_character_classes(weights)) < 3 for weights in unrestricted)
    for weights in full_span + unrestricted:
        assert np.max(np.abs(weights)) <= 1
        assert np.linalg.matrix_rank(weights) == 2
        assert all(any(column) for column in zip(*weights))
    if method == "span-equivalence":
        # The target rows (1, 0, 1) and (0, 1, -1) span this plane.
        normal = np.array([-1, 1, 1])
        assert any(np.all(np.asarray(weights) @ normal == 0) for weights in full_span)
