"""Compare the Phase 1 speedups with the previous mathematical calculations."""

import itertools
from math import gcd

import numpy as np
import pytest
from scipy.linalg import eigh, expm
from scipy.stats import qmc, special_ortho_group
from sklearn.neighbors import KDTree

from gudhi.liedetect.algebra import (
    are_representations_equivalent,
    get_canonical_pushforward_algebra,
    get_lattices,
    invariant_of_lattices,
)
from gudhi.liedetect.liepca import get_lie_pca_operator
from gudhi.liedetect.utils import sample_orbit_from_algebra_su2


@pytest.mark.parametrize("method", ["PCA", "covariance"])
@pytest.mark.parametrize("correction", [False, True])
def test_lie_pca_assembly_matches_dense_reference(method, correction):
    points = np.random.default_rng(16).normal(size=(20, 4))
    neighbors = KDTree(points).query(points, 6, return_distance=False)[:, 1:]
    normal_projectors = []
    line_projectors = []
    for point, indices in zip(points, neighbors):
        differences = point - points[indices]
        covariance = differences.T @ differences / len(indices)
        if method == "PCA":
            _, vectors = eigh(covariance)
            normal = vectors[:, :2] @ vectors[:, :2].T
        else:
            normal = np.eye(4) - covariance * np.sqrt(2) / np.linalg.norm(covariance)
        normal_projectors.append(normal)
        line_projectors.append(np.outer(point, point) / (point @ point))

    columns = []
    for basis in np.eye(16).reshape(16, 4, 4):
        symmetric = (basis + basis.T) / 2
        argument = (basis - basis.T) / 2 if correction else basis
        image = np.mean(
            [normal @ argument @ line for normal, line in zip(normal_projectors, line_projectors)],
            axis=0,
        )
        columns.append((image + symmetric if correction else image).ravel())
    expected = np.array(columns).T
    if correction:
        expected = (expected + expected.T) / 2

    computed = get_lie_pca_operator(points, 5, 2, method=method, correction=correction)
    np.testing.assert_allclose(computed, expected, atol=1e-13)


def _legacy_invariant(weights, method):
    if len(weights) == 1:
        return (tuple(sorted(map(abs, weights[0]))),)

    def projection(rows):
        rows = rows / np.linalg.norm(rows, axis=1, keepdims=True)
        basis, _ = np.linalg.qr(rows.T)
        return tuple(np.round(np.sum([np.outer(v, v) for v in basis.T], axis=0).ravel(), 5))

    weights = np.asarray(weights)
    if method == "span-equivalence":
        return projection(weights)
    width = weights.shape[1]
    projections = [
        projection(weights[:, permutation] * signs)
        for permutation in itertools.permutations(range(width))
        for signs in itertools.product([-1, 1], repeat=width)
    ]
    return min(projections)


@pytest.mark.parametrize("rank", [1, 2])
@pytest.mark.parametrize("span", [False, True])
def test_lattice_streaming_preserves_representatives_and_order(rank, span):
    width, bound = 3, 2 if rank == 1 else 1
    characters = list(itertools.product(range(bound + 1), repeat=rank))[1:]
    blocks = itertools.combinations(characters, width) if span else itertools.product(characters, repeat=width)
    candidates = [tuple(zip(*candidate)) for candidate in blocks]
    if rank == 1:
        candidates = [candidate for candidate in candidates if gcd(*candidate[0]) == 1]
    representatives = {}
    for candidate in candidates:
        if np.linalg.matrix_rank(candidate) == rank:
            representatives.setdefault(_legacy_invariant(candidate, "span-equivalence"), candidate)
    assert get_lattices(rank, width, bound, "span-equivalence", span) == list(representatives.values())
    orbits = {}
    for candidate in representatives.values():
        expected = _legacy_invariant(candidate, "orbit-equivalence")
        assert invariant_of_lattices(candidate, "orbit-equivalence") == expected
        orbits.setdefault(expected, candidate)
    assert get_lattices(rank, width, bound, "orbit-equivalence", span) == list(orbits.values())


def test_rank_two_torus_equivalence_distinguishes_three_column_images():
    weights = ((1, 0, 1), (0, 1, 1))
    assert are_representations_equivalent("torus", weights, ((-1, 1, 0), (-1, 0, 1)))
    assert not are_representations_equivalent("torus", weights, ((1, 0, 1), (0, 1, 2)))


@pytest.mark.parametrize(("group", "rep_type"), [("SO(3)", (3, 5)), ("SU(2)", (4,))])
@pytest.mark.parametrize("method", ["evenly_spaced", "random_uniform"])
def test_spectral_sampling_matches_dense_rotations(group, rep_type, method):
    rng = np.random.default_rng(3)
    dimension = sum(rep_type)
    orthogonal = special_ortho_group.rvs(dimension, random_state=rng)
    algebra = [
        orthogonal @ generator @ orthogonal.T
        for generator in get_canonical_pushforward_algebra(group, rep_type)
    ]
    _, y, z = algebra
    x = rng.normal(size=dimension)
    count = 24
    if method == "evenly_spaced":
        parameters = qmc.Halton(d=3, scramble=False).random(count)
    else:
        parameters = np.random.default_rng(17).random((count, 3))
    alpha = 2 * np.pi * parameters[:, 0]
    beta = np.arccos(2 * parameters[:, 1] - 1)
    gamma = (4 if group == "SU(2)" else 2) * np.pi * parameters[:, 2]
    expected = np.array([expm(a * z) @ expm(b * y) @ expm(c * z) @ x for a, b, c in zip(alpha, beta, gamma)])
    computed = sample_orbit_from_algebra_su2(rep_type, algebra, x, count, method=method, seed=17)
    np.testing.assert_allclose(computed, expected, atol=1e-12)
    np.testing.assert_allclose(np.linalg.norm(computed, axis=1), np.linalg.norm(x), atol=1e-12)
