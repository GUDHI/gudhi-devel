"""Regression tests for the public LieDetect sampling and preprocessing APIs."""

import numpy as np
import pytest

from gudhi.datasets.linear_orbits import sample_from_lie_group_rep
from gudhi.liedetect import Orthonormalize


@pytest.mark.parametrize(("group", "rep_type"), [("SU(2)", (4,)), ("SO(3)", (3,))])
def test_nonabelian_sampling_method(group, rep_type):
    random_sample = sample_from_lie_group_rep(
        group, rep_type, 16, method="random_uniform", seed=1
    )
    np.testing.assert_array_equal(
        random_sample,
        sample_from_lie_group_rep(group, rep_type, 16, method="random_uniform", seed=1),
    )
    assert not np.allclose(
        random_sample,
        sample_from_lie_group_rep(group, rep_type, 16, method="random_uniform", seed=2),
    )
    np.testing.assert_array_equal(
        sample_from_lie_group_rep(group, rep_type, 16, method="evenly_spaced", seed=1),
        sample_from_lie_group_rep(group, rep_type, 16, method="evenly_spaced", seed=2),
    )
    with pytest.raises(ValueError, match="Method .* not recognized"):
        sample_from_lie_group_rep(group, rep_type, 16, method="invalid")


def test_inverse_transform_uses_transformed_training_points():
    angles = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    points = np.column_stack((10 + 3 * np.cos(angles), -4 + 7 * np.sin(angles)))
    normalizer = Orthonormalize(points)
    normalizer.fit()

    np.testing.assert_allclose(normalizer.inverse_transform(), points, atol=1e-12)
    np.testing.assert_allclose(
        normalizer.inverse_transform(add_mean=False),
        points - points.mean(axis=0),
        atol=1e-12,
    )

    transformed = normalizer.transform(points[:4])
    restored = normalizer.inverse_transform(transformed, in_place=True)
    np.testing.assert_allclose(restored, points[:4], atol=1e-12)
    np.testing.assert_array_equal(normalizer.points, restored)
