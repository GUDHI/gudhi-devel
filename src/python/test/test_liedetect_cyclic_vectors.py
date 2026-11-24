"""Independent regressions for centered span and non-Abelian multiplicities."""

import itertools

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from gudhi.datasets.linear_orbits import sample_from_lie_group, sample_from_lie_group_rep
from gudhi.liedetect import OrbitFitter
from gudhi.liedetect.algebra import get_constrained_partitions


def _signed_permutation_rotations():
    for permutation in itertools.permutations(range(3)):
        for signs in itertools.product((-1, 1), repeat=3):
            rotation = np.eye(3)[:, permutation] * signs
            if np.linalg.det(rotation) > 0:
                yield rotation


def _centered_rank(points):
    return np.linalg.matrix_rank(points - points.mean(axis=0))


def test_so3_stiefel_orbit_spans_both_copies():
    rotations = np.array(list(_signed_permutation_rotations()))
    independent = np.column_stack((rotations[:, :, 0], rotations[:, :, 1]))
    proportional = np.column_stack((rotations[:, :, 0], 2 * rotations[:, :, 0]))
    assert np.linalg.matrix_rank(independent) == _centered_rank(independent) == 6
    assert np.linalg.matrix_rank(proportional) == _centered_rank(proportional) == 3

    # Use an independently generated continuous-group sample for tangent fitting.
    rotations = Rotation.random(240, random_state=11).as_matrix()
    points = np.column_stack((rotations[:, :, 0], rotations[:, :, 1]))
    fitter = OrbitFitter(points)
    fitter.lie_pca(nb_neighbors=20, orbit_dim=3)
    rep_type, _ = fitter.closest_algebra("SO(3)", seed=16)
    assert rep_type == (3, 3)


@pytest.mark.parametrize(
    ("group", "partition", "has_cyclic_vector"),
    [
        ("SO(3)", (3, 3, 3), True),
        ("SO(3)", (3, 3, 3, 3), False),
        ("SU(2)", (4, 4), False),
        ("SU(2)", (8, 8), True),
    ],
)
def test_multiplicity_bound(group, partition, has_cyclic_vector):
    dimension = sum(partition)
    assert partition in get_constrained_partitions(group, dimension)
    restricted = get_constrained_partitions(group, dimension, span_ambient_space=True)
    assert (partition in restricted) == has_cyclic_vector


def test_trivial_component_can_increase_raw_span_only():
    points = np.array([np.r_[1, q[:, 0]] for q in _signed_permutation_rotations()])
    assert np.linalg.matrix_rank(points) == 4
    assert _centered_rank(points) == 3
    assert (1, 3) in get_constrained_partitions("SO(3)", 4)
    with pytest.raises(ValueError, match="No partitions"):
        get_constrained_partitions("SO(3)", 4, span_ambient_space=True)


def test_random_full_span_sampler_uses_cyclic_base_point():
    def sample(x=None):
        return sample_from_lie_group(
            "SO(3)", 6, 48, x=x, span_ambient_space=True, seed=42
        )

    points, rep_type = sample()
    assert rep_type == (3, 3)
    assert _centered_rank(points) == 6
    np.testing.assert_array_equal(points, sample()[0])

    independent = np.r_[1, 0, 0, 0, 1, 0]
    assert _centered_rank(sample(independent)[0]) == 6

    proportional = np.r_[1, 0, 0, 2, 0, 0]
    with pytest.raises(ValueError, match="must be cyclic"):
        sample(proportional)
    points = sample_from_lie_group_rep("SO(3)", (3, 3), 48, x=proportional)
    assert _centered_rank(points) == 3
