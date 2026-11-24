"""Mathematical regressions for random torus weight presentations."""

import itertools

import numpy as np
import pytest

from gudhi.datasets.linear_orbits import sample_from_lie_group, sample_from_lie_group_rep
from gudhi.liedetect import algebra, are_representations_equivalent


@pytest.mark.parametrize(("rank", "bound"), [(1, 4), (2, 1), (2, 2), (3, 1)])
def test_balanced_radix_indexes_exactly_the_nonzero_sign_classes(rank, bound):
    signed_box = itertools.product(range(-bound, bound + 1), repeat=rank)
    canonical = [weight for weight in signed_box if weight > (0,) * rank]
    assert [
        algebra._decode_weight_class(index, rank, bound)
        for index in range(1, len(canonical) + 1)
    ] == canonical


def test_sparse_draw_enumerates_each_ordered_pair_once():
    class ScriptedRng:
        def __init__(self, indices):
            self.indices = iter(indices)

        def integers(self, stop):
            index = next(self.indices)
            assert 0 <= index < stop
            return index

    outcomes = [
        tuple(algebra._sample_weight_class_indices(ScriptedRng(indices), 4, 2))
        for indices in itertools.product(range(4), range(3))
    ]
    assert sorted(outcomes) == sorted(itertools.permutations(range(1, 5), 2))


def test_sparse_draw_supports_pools_beyond_int64():
    class ByteOnlyRng:
        def __init__(self):
            self.rng = np.random.default_rng(4)
            self.requested_sizes = []

        def bytes(self, size):
            self.requested_sizes.append(size)
            return self.rng.bytes(size)

    stop = 10**40
    rng = ByteOnlyRng()
    indices = list(algebra._sample_weight_class_indices(rng, stop, 6))
    assert len(set(indices)) == 6
    assert all(type(index) is int and 1 <= index <= stop for index in indices)
    assert rng.requested_sizes and max(rng.requested_sizes) <= 17
    assert indices == list(algebra._sample_weight_class_indices(ByteOnlyRng(), stop, 6))


def test_full_span_sampler_accepts_nonprimitive_circle_weights(monkeypatch):
    def indices(rng, population, count):
        assert (population, count) == (4, 2)
        return iter((2, 4))

    monkeypatch.setattr(algebra, "_sample_weight_class_indices", indices)
    weights = algebra.get_random_lattice(1, 2, 4, span_ambient_space=True, seed=3)
    assert tuple(map(abs, weights[0])) == (2, 4)
    assert are_representations_equivalent("torus", weights, ((1, 2),))

    angles = 2 * np.pi * np.arange(17) / 17
    phases = angles[:, None] * np.array(weights[0])
    points = np.stack((np.cos(phases), np.sin(phases)), axis=-1).reshape(-1, 4)
    np.testing.assert_allclose(points.mean(axis=0), 0, atol=1e-15)
    np.testing.assert_allclose(points.T @ points / len(points), np.eye(4) / 2, atol=1e-15)


def test_rank_two_full_span_sampling_boundary():
    weights = algebra.get_random_lattice(2, 4, 1, span_ambient_space=True, seed=19)
    columns = list(zip(*weights))
    assert {max(w, tuple(-value for value in w)) for w in columns} == {
        (0, 1), (1, -1), (1, 0), (1, 1)
    }
    assert np.linalg.matrix_rank(weights) == 2
    with pytest.raises(ValueError, match="distinct nonzero frequency"):
        algebra.get_random_lattice(2, 5, 1, span_ambient_space=True, seed=19)


def test_exact_rank_accepts_a_poorly_conditioned_unimodular_table(monkeypatch):
    k = 10**9
    weights = ((k, k + 1), (k - 1, k))
    assert weights[0][0] * weights[1][1] - weights[0][1] * weights[1][0] == 1

    class ScriptedRng:
        def integers(self, low, high, *, size, endpoint):
            assert size == (2, 2)
            return np.array(weights, dtype=np.int64)

    monkeypatch.setattr(algebra, "default_rng", lambda seed: ScriptedRng())
    assert algebra.get_random_lattice(2, 2, k + 1) == weights


def test_full_span_rejection_has_a_bound(monkeypatch):
    calls = 0

    def collinear_indices(rng, population, count):
        nonlocal calls
        calls += 1
        assert (population, count) == (12, 2)
        # In balanced base five these decode to (1, 0) and (2, 0).
        return iter((5, 10))

    monkeypatch.setattr(algebra, "_sample_weight_class_indices", collinear_indices)
    with pytest.raises(RuntimeError, match="100 attempts"):
        algebra.get_random_lattice(2, 2, 2, span_ambient_space=True, seed=2)
    assert calls == 100


def test_public_sampler_validates_the_base_point_and_is_reproducible():
    options = dict(frequency_max=1, group_dim=2, method="random_uniform", seed=9)
    points, weights = sample_from_lie_group("torus", 8, 48, span_ambient_space=True, **options)
    again, again_weights = sample_from_lie_group("torus", 8, 48, span_ambient_space=True, **options)
    np.testing.assert_array_equal(points, again)
    assert weights == again_weights
    assert np.linalg.matrix_rank(points - points.mean(axis=0)) == 8

    with_zero_plane = np.array([1, 0, 1, 0, 1, 0, 0, 0])
    with pytest.raises(ValueError, match="nonzero in every rotation plane"):
        sample_from_lie_group("torus", 8, 48, x=with_zero_plane, span_ambient_space=True, **options)
    lower_span, _ = sample_from_lie_group("torus", 8, 48, x=with_zero_plane, **options)
    assert lower_span.shape == (48, 8)
    explicit = sample_from_lie_group_rep(
        "torus", weights, 48, x=with_zero_plane, method="random_uniform", seed=9
    )
    assert np.linalg.matrix_rank(explicit - explicit.mean(axis=0)) == 6
