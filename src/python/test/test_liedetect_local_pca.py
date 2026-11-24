"""Independent regressions for the neighborhood SVD and numerical rank check."""

import numpy as np
import pytest

from gudhi.liedetect.liepca import get_lie_pca_operator


def test_direct_pca_retains_a_direction_lost_in_the_formed_gram_matrix():
    t = 1e-9
    differences = np.array([[1.0, 1, 0], [t, -t, 0]])
    assert np.linalg.matrix_rank(differences.T @ differences) == 1
    assert np.linalg.matrix_rank(differences) == 2

    origin = np.array([0.0, 0, 1])
    points = np.vstack([origin, origin - differences])
    computed = get_lie_pca_operator(points, nb_neighbors=2, orbit_dim=2)

    # Each neighborhood consists of the two other vertices of the same
    # noncollinear triangle in z=1. Its exact normal projector is E_33.
    normal = np.diag([0.0, 0, 1])
    lines = [np.outer(point, point) / (point @ point) for point in points]
    expected = np.column_stack([
        np.mean([normal @ basis @ line for line in lines], axis=0).ravel()
        for basis in np.eye(9).reshape(9, 3, 3)
    ])
    np.testing.assert_allclose(computed, expected, atol=1e-14)


def test_pca_rejects_a_collinear_neighborhood_but_covariance_remains_available():
    points = np.array([[1, 0.0, 0], [1, 0.01, 0], [1, 0.03, 0], [1, 0.08, 0]])
    with pytest.raises(ValueError, match="point 0 has numerical rank 1, below orbit_dim=2"):
        get_lie_pca_operator(points, nb_neighbors=2, orbit_dim=2, method="PCA")

    covariance = get_lie_pca_operator(points, nb_neighbors=2, orbit_dim=2, method="covariance")
    assert np.linalg.eigvalsh(covariance)[0] < 0


def test_pca_rejects_an_algebraically_nonzero_direction_below_the_numerical_floor():
    origin = np.array([0.0, 0, 1])
    points = np.vstack([origin, origin - [1, 1, 0], origin - [1e-16, -1e-16, 0]])
    with pytest.raises(ValueError, match="point 0 has numerical rank 1, below orbit_dim=2"):
        get_lie_pca_operator(points, nb_neighbors=2, orbit_dim=2)
