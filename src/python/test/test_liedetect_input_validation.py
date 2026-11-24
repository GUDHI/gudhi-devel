"""Regression tests for LieDetect preprocessing and basic input validation."""

import numpy as np
import pytest
from scipy.linalg import sqrtm

from gudhi.liedetect import Orthonormalize
from gudhi.liedetect.liepca import get_lie_pca_operator


def test_whitening_integer_points_matches_symmetric_inverse_square_root():
    points = np.array([[5, 7], [9, 7], [6, 12], [8, 15]])
    centered = points - points.mean(axis=0)
    covariance = centered.T @ centered / len(points)
    reference = centered @ sqrtm(np.linalg.inv(covariance))
    reference /= np.linalg.norm(reference, axis=1).mean()

    normalizer = Orthonormalize(points)
    np.testing.assert_allclose(normalizer.fit_transform(), reference, atol=1e-12)
    np.testing.assert_allclose(normalizer.inverse_transform(), points, atol=1e-12)


def test_whitening_reports_singular_covariance():
    points = np.array([[0, 0], [1, 0], [2, 0]])
    with pytest.raises(ValueError, match="full-rank centered covariance"):
        Orthonormalize(points).fit()


def test_whitening_rejects_complex_coordinates():
    with pytest.raises(ValueError, match="real coordinates"):
        Orthonormalize(np.array([[1 + 2j, 0], [1, 3], [2, 4]]))


@pytest.mark.parametrize(
    ("points", "nb_neighbors", "orbit_dim", "message"),
    [
        (np.empty((0, 2)), 1, 1, "at least two points"),
        ([[1, 0], [0, 1]], 2, 1, "nb_neighbors"),
        ([[1, 0], [0, 1]], 1, 0, "orbit_dim"),
        ([[1, 0], [0, 1]], 1, 1.5, "orbit_dim"),
        ([[1, 0], [0, np.inf]], 1, 1, "coordinates must be finite"),
        ([[1, 0], [0, 0]], 1, 1, "nonzero squared norms"),
        ([[1, 0], [1, 0]], 1, 1, "Local covariance must have"),
    ],
)
def test_lie_pca_invalid_input(points, nb_neighbors, orbit_dim, message):
    with pytest.raises(ValueError, match=message):
        get_lie_pca_operator(points, nb_neighbors, orbit_dim)


def test_lie_pca_rejects_ambient_dimension_as_tangent_dimension():
    points = np.array([[1, 0], [0, 1], [-1, 0]])
    with pytest.raises(ValueError, match="strictly smaller than the ambient dimension"):
        get_lie_pca_operator(points, nb_neighbors=1, orbit_dim=2)
