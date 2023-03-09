""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Theo Lacombe, Marc Glisse

    Copyright (C) 2019 Inria

    Modification(s):
      - 2020/07 Théo Lacombe: Added tests about handling essential parts in diagrams.
      - YYYY/MM Author: Description of the modification
"""

from gudhi.wasserstein.wasserstein import _proj_on_diag, _finite_part, _handle_essential_parts, _get_essential_parts
from gudhi.wasserstein.wasserstein import _warn_infty
from gudhi.wasserstein import wasserstein_distance as pot
from gudhi.hera import wasserstein_distance as hera
import numpy as np
import pytest


__author__ = "Theo Lacombe"
__copyright__ = "Copyright (C) 2019 Inria"
__license__ = "MIT"


def test_proj_on_diag():
    dgm = np.array([[1., 1.], [1., 2.], [3., 5.]])
    assert np.array_equal(_proj_on_diag(dgm), [[1., 1.], [1.5, 1.5], [4., 4.]])
    empty = np.empty((0, 2))
    assert np.array_equal(_proj_on_diag(empty), empty)


def test_finite_part():
    diag = np.array([[0, 1], [3, 5], [2, np.inf], [3, np.inf], [-np.inf, 8], [-np.inf, 12], [-np.inf, -np.inf],
                     [np.inf, np.inf], [-np.inf, np.inf], [-np.inf, np.inf]])
    assert np.array_equal(_finite_part(diag), [[0, 1], [3, 5]])


def test_handle_essential_parts():
    diag1 = np.array([[0, 1], [3, 5],
                      [2, np.inf], [3, np.inf],
                      [-np.inf, 8], [-np.inf, 12],
                      [-np.inf, -np.inf],
                      [np.inf, np.inf],
                      [-np.inf, np.inf], [-np.inf, np.inf]])

    diag2 = np.array([[0, 2], [3, 5],
                      [2, np.inf], [4, np.inf],
                      [-np.inf, 8], [-np.inf, 11],
                      [-np.inf, -np.inf],
                      [np.inf, np.inf],
                      [-np.inf, np.inf], [-np.inf, np.inf]])

    diag3 = np.array([[0, 2], [3, 5],
                      [2, np.inf], [4, np.inf], [6, np.inf],
                      [-np.inf, 8], [-np.inf, 11],
                      [-np.inf, -np.inf],
                      [np.inf, np.inf],
                      [-np.inf, np.inf], [-np.inf, np.inf]])

    c, m = _handle_essential_parts(diag1, diag2, order=1)
    assert c == pytest.approx(2, 0.0001)  # Note: here c is only the cost due to essential part (thus 2, not 3)
    # Similarly, the matching only corresponds to essential parts.
    # Note that (-inf,-inf) and (+inf,+inf) coordinates are matched to the diagonal.
    assert np.array_equal(m, [[4, 4], [5, 5], [2, 2], [3, 3], [8, 8], [9, 9], [6, -1], [7, -1], [-1, 6], [-1, 7]])

    c, m = _handle_essential_parts(diag1, diag3, order=1)
    assert c == np.inf
    assert (m is None)


def test_get_essential_parts():
    diag1 = np.array([[0, 1], [3, 5], [2, np.inf], [3, np.inf], [-np.inf, 8], [-np.inf, 12], [-np.inf, -np.inf],
                     [np.inf, np.inf], [-np.inf, np.inf], [-np.inf, np.inf]])

    diag2 = np.array([[0, 1], [3, 5], [2, np.inf], [3, np.inf]])

    res  = _get_essential_parts(diag1)
    res2 = _get_essential_parts(diag2)
    assert np.array_equal(res[0], [4, 5])
    assert np.array_equal(res[1], [2, 3])
    assert np.array_equal(res[2], [8, 9])
    assert np.array_equal(res[3], [6]   )
    assert np.array_equal(res[4], [7]   )

    assert np.array_equal(res2[0], []    )
    assert np.array_equal(res2[1], [2, 3])
    assert np.array_equal(res2[2], []    )
    assert np.array_equal(res2[3], []    )
    assert np.array_equal(res2[4], []    )


def test_warn_infty():
    with pytest.warns(UserWarning):
        assert _warn_infty(matching=False)==np.inf
        c, m = _warn_infty(matching=True)
        assert (c == np.inf)
        assert (m is None)


def _to_set(X):
    return { (i, j) for i, j in X }

def _same_permuted(X, Y):
    return _to_set(X) == _to_set(Y)


def _basic_wasserstein(wasserstein_distance, delta, test_infinity=True, test_matching=True):
    diag1 = np.array([[2.7, 3.7], [9.6, 14.0], [34.2, 34.974]])
    diag2 = np.array([[2.8, 4.45], [9.5, 14.1]])
    diag3 = np.array([[0, 2], [4, 6]])
    diag4 = np.array([[0, 3], [4, 8]])
    emptydiag = np.array([])

    # We just need to handle positive numbers here
    def approx(x):
        return pytest.approx(x, rel=delta)

    assert wasserstein_distance(emptydiag, emptydiag, internal_p=2.,     order=1.) == 0.
    assert wasserstein_distance(emptydiag, emptydiag, internal_p=np.inf, order=1.) == 0.
    assert wasserstein_distance(emptydiag, emptydiag, internal_p=np.inf, order=2.) == 0.
    assert wasserstein_distance(emptydiag, emptydiag, internal_p=2.,     order=2.) == 0.

    assert wasserstein_distance(diag3, emptydiag, internal_p=np.inf,     order=1.) == approx(2.)
    assert wasserstein_distance(diag3, emptydiag, internal_p=1.,         order=1.) == approx(4.)

    assert wasserstein_distance(diag4, emptydiag, internal_p=1.,     order=2.) == approx(5.)  # thank you Pythagorician triplets
    assert wasserstein_distance(diag4, emptydiag, internal_p=np.inf, order=2.) == approx(2.5)
    assert wasserstein_distance(diag4, emptydiag, internal_p=2.,     order=2.) == approx(3.5355339059327378)

    assert wasserstein_distance(diag1, diag2, internal_p=2.,   order=1.)   == approx(1.4453593023967701)
    assert wasserstein_distance(diag1, diag2, internal_p=2.35, order=1.74) == approx(0.9772734057168739)

    assert wasserstein_distance(diag1, emptydiag, internal_p=2.35, order=1.7863) == approx(3.141592214572228)

    assert wasserstein_distance(diag3, diag4, internal_p=1.,     order=1.) == approx(3.)
    assert wasserstein_distance(diag3, diag4, internal_p=np.inf, order=1.) == approx(3.)  # no diag matching here
    assert wasserstein_distance(diag3, diag4, internal_p=np.inf, order=2.) == approx(np.sqrt(5))
    assert wasserstein_distance(diag3, diag4, internal_p=1.,     order=2.) == approx(np.sqrt(5))
    assert wasserstein_distance(diag3, diag4, internal_p=4.5,    order=2.) == approx(np.sqrt(5))

    if test_infinity:
        diag5 = np.array([[0, 3], [4, np.inf]])
        diag6 = np.array([[7, 8], [4, 6], [3, np.inf]])

        assert wasserstein_distance(diag4, diag5) == np.inf
        assert wasserstein_distance(diag5, diag6, order=1, internal_p=np.inf) == approx(4.)
        assert wasserstein_distance(diag5, emptydiag) == np.inf

    if test_matching:
        match = wasserstein_distance(emptydiag, emptydiag, matching=True, internal_p=1., order=2)[1]
        # Accept [] or np.array of shape (2, 0)
        assert len(match) == 0
        match = wasserstein_distance(emptydiag, emptydiag, matching=True, internal_p=np.inf, order=2.24)[1]
        assert len(match) == 0
        match = wasserstein_distance(emptydiag, diag2, matching=True, internal_p=np.inf, order=2.)[1]
        assert _same_permuted(match, [[-1, 0], [-1, 1]])
        match = wasserstein_distance(diag2, emptydiag, matching=True, internal_p=np.inf, order=2.24)[1]
        assert _same_permuted(match, [[0, -1], [1, -1]])
        match = wasserstein_distance(diag1, diag2, matching=True, internal_p=2., order=2.)[1]
        assert _same_permuted(match, [[0, 0], [1, 1], [2, -1]])

    if test_matching and test_infinity:
        diag7 = np.array([[0, 3], [4, np.inf], [5, np.inf]])
        diag8 = np.array([[0,1], [0, np.inf], [-np.inf, -np.inf], [np.inf, np.inf]])
        diag9 = np.array([[-np.inf, -np.inf], [np.inf, np.inf]])
        diag10 = np.array([[0,1], [-np.inf, -np.inf], [np.inf, np.inf]])

        match = wasserstein_distance(diag5, diag6, matching=True, internal_p=2., order=2.)[1]
        assert _same_permuted(match, [[0, -1], [-1,0], [-1, 1], [1, 2]])
        match = wasserstein_distance(diag5, diag7, matching=True, internal_p=2., order=2.)[1]
        assert (match is None)
        cost, match = wasserstein_distance(diag7, emptydiag, matching=True, internal_p=2., order=2.3)
        assert (cost == np.inf)
        assert (match is None)
        cost, match = wasserstein_distance(emptydiag, diag7, matching=True, internal_p=2.42, order=2.)
        assert (cost == np.inf)
        assert (match is None)
        cost, match = wasserstein_distance(diag8, diag9, matching=True, internal_p=2., order=2.)
        assert (cost == np.inf)
        assert (match is None)
        cost, match = wasserstein_distance(diag9, diag10, matching=True, internal_p=1., order=1.)
        assert (cost == 1)
        assert _same_permuted(match, [[0, -1],[1, -1],[-1, 0], [-1, 1], [-1, 2]]) # type 4 and 5 are match to the diag anyway.
        cost, match = wasserstein_distance(diag9, emptydiag, matching=True, internal_p=2., order=2.)
        assert (cost == 0.)
        assert _same_permuted(match, [[0, -1], [1, -1]])


def hera_wrap(**extra):
    def fun(*kargs,**kwargs):
        return hera(*kargs,**kwargs,**extra)
    return fun


def pot_wrap(**extra):
    def fun(*kargs,**kwargs):
        return pot(*kargs,**kwargs,**extra)
    return fun


def test_wasserstein_distance_pot():
    _basic_wasserstein(pot, 1e-15, test_infinity=False, test_matching=True)  # pot with its standard args
    _basic_wasserstein(pot_wrap(enable_autodiff=True, keep_essential_parts=False), 1e-15, test_infinity=False, test_matching=False)


def test_wasserstein_distance_hera():
    _basic_wasserstein(hera_wrap(delta=1e-12), 1e-12, test_matching=True)
    _basic_wasserstein(hera_wrap(delta=.1), .1, test_matching=True)

