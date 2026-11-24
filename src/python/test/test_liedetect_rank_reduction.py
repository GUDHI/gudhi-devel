"""Independent regressions for affine rank reduction and geometric input checks."""

import numpy as np
import pytest

from gudhi.liedetect import OrbitFitter, Orthonormalize
from gudhi.liedetect.liepca import get_lie_pca_operator, orthonormalize


def _embedded_ellipse():
    angles = np.linspace(0, 2 * np.pi, 96, endpoint=False)
    a = 3 * np.cos(angles)
    b = 7 * np.sin(angles)
    return np.column_stack((a, b, a + b)) + [10, -4, 2]


def test_oblique_affine_ellipse_reduces_and_fits():
    points = _embedded_ellipse()
    normalizer = Orthonormalize(points, rank_tol=0)
    reduced = normalizer.fit_transform()
    assert reduced.shape == (len(points), 2)
    np.testing.assert_allclose(reduced.mean(axis=0), 0, atol=1e-14)
    np.testing.assert_allclose(reduced.T @ reduced / len(points), np.eye(2) / 2, atol=1e-14)
    assert np.linalg.norm(reduced, axis=1).mean() == pytest.approx(1)
    np.testing.assert_allclose(normalizer.transform(points), reduced, atol=1e-14)
    np.testing.assert_allclose(normalizer.inverse_transform(), points, atol=1e-13)

    fitter = OrbitFitter(reduced)
    fitter.lie_pca(8, 1)
    rep_type, _ = fitter.closest_algebra("torus", group_dim=1, frequency_max=2, seed=1)
    assert rep_type == ((1,),)
    reconstructed = normalizer.inverse_transform(fitter.sample_orbit(96))
    assert reconstructed.shape == points.shape
    np.testing.assert_allclose(
        (reconstructed - points.mean(axis=0)) @ np.array([1, 1, -1]),
        0,
        atol=1e-13,
    )


def test_relative_threshold_and_affine_projection():
    # Centered covariance is diag(1, 0.01, 0) independent of translation.
    points = np.array([[1, 0.1, 0], [1, -0.1, 0], [-1, 0.1, 0], [-1, -0.1, 0]])
    mean = np.array([3.0, 4, 5])
    points += mean
    for threshold, dimension in ((0.009, 2), (0.011, 1)):
        reduced, forward = orthonormalize(points, rank_tol=threshold)
        assert reduced.shape == (4, dimension)
        assert forward.shape == (dimension, 3)

    normalizer = Orthonormalize(points, rank_tol=0.011)
    normalizer.fit()
    query = np.array([[9.0, 8, 7]])
    transformed = normalizer.transform(query)
    np.testing.assert_allclose(normalizer.inverse_transform(transformed), [[9, 4, 5]])
    np.testing.assert_allclose(
        normalizer.inverse_transform(transformed, add_mean=False), [[6, 0, 0]]
    )
    with pytest.raises(ValueError, match="1 fitted coordinates"):
        normalizer.inverse_transform(query)


def test_underdetermined_whitening_and_zero_rank_error():
    points = np.array([[1, 2, 0, 0, 0], [2, 3, 0, 0, 0], [3, 2, 0, 0, 0]])
    normalizer = Orthonormalize(points, rank_tol=0)
    assert normalizer.fit_transform().shape == (3, 2)
    np.testing.assert_allclose(normalizer.inverse_transform(), points, atol=1e-13)
    with pytest.raises(ValueError, match="no retained directions"):
        Orthonormalize(np.ones((4, 3)), rank_tol=0).fit()


def test_refit_after_in_place_lift_uses_new_mean():
    points = _embedded_ellipse()
    normalizer = Orthonormalize(points, rank_tol=0)
    reduced = normalizer.fit_transform()
    lifted = normalizer.inverse_transform(reduced + [2, -1], in_place=True)
    np.testing.assert_array_equal(normalizer.points, lifted)
    assert not np.allclose(lifted.mean(axis=0), points.mean(axis=0))
    # The previous fitted transform remains usable until an explicit refit.
    np.testing.assert_allclose(normalizer.inverse_transform(reduced), points, atol=1e-13)
    refitted = normalizer.fit_transform()
    np.testing.assert_allclose(refitted.mean(axis=0), 0, atol=1e-14)
    np.testing.assert_allclose(normalizer.inverse_transform(), lifted, atol=1e-13)


@pytest.mark.parametrize("rank_tol", [-0.01, 1, np.nan])
def test_rank_tol_must_be_a_relative_threshold(rank_tol):
    with pytest.raises(ValueError, match="relative covariance threshold"):
        Orthonormalize(_embedded_ellipse(), rank_tol=rank_tol).fit()


@pytest.mark.parametrize("method", ["PCA", "covariance"])
def test_neighbor_count_must_support_orbit_dimension(method):
    with pytest.raises(ValueError, match="'nb_neighbors' must be at least 'orbit_dim'"):
        get_lie_pca_operator(np.eye(3), nb_neighbors=1, orbit_dim=2, method=method)


def test_tied_duplicate_neighborhood_matches_direct_normal_projection(monkeypatch):
    points = np.array([[1.0, 0], [1, 0], [0, 1]])
    # The first two queries both report the duplicate with index 1 first.
    class TiedTree:
        def __init__(self, *args, **kwargs):
            pass

        def query(self, *args, **kwargs):
            return np.array([[1, 0, 2], [1, 0, 2], [2, 1, 0]])

    monkeypatch.setattr("sklearn.neighbors.KDTree", TiedTree)
    actual = get_lie_pca_operator(points, nb_neighbors=2, orbit_dim=1)
    # All nonzero differences are parallel to (1,-1).
    normal = np.ones((2, 2)) / 2
    line_projectors = [np.outer(point, point) for point in points]
    expected = np.column_stack([
        np.mean([normal @ basis @ line for line in line_projectors], axis=0).ravel()
        for basis in np.eye(4).reshape(4, 2, 2)
    ])
    np.testing.assert_allclose(actual, expected, atol=1e-14)
