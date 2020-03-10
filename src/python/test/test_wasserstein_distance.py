""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Theo Lacombe, Marc Glisse

    Copyright (C) 2019 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.wasserstein import wasserstein_distance as pot
from gudhi.hera import wasserstein_distance as hera
import numpy as np
import pytest

__author__ = "Theo Lacombe"
__copyright__ = "Copyright (C) 2019 Inria"
__license__ = "MIT"

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


    if test_matching:
        match = wasserstein_distance(emptydiag, emptydiag, matching=True, internal_p=1., order=2)[1]
        assert np.array_equal(match, np.array([]))
        match = wasserstein_distance(emptydiag, emptydiag, matching=True, internal_p=np.inf, order=2.24)[1]
        assert np.array_equal(match, np.array([]))
        match = wasserstein_distance(emptydiag, diag2, matching=True, internal_p=np.inf, order=2.)[1]
        assert np.array_equal(match , np.array([[-1, 0], [-1, 1]]))
        match = wasserstein_distance(diag2, emptydiag, matching=True, internal_p=np.inf, order=2.24)[1]
        assert np.array_equal(match , np.array([[0, -1], [1, -1]]))
        match = wasserstein_distance(diag1, diag2, matching=True, internal_p=2., order=2.)[1]
        assert np.array_equal(match, np.array([[0, 0], [1, 1], [2, -1]]))
        


def hera_wrap(delta):
    def fun(*kargs,**kwargs):
        return hera(*kargs,**kwargs,delta=delta)
    return fun

def test_wasserstein_distance_pot():
    _basic_wasserstein(pot, 1e-15, test_infinity=False, test_matching=True)

def test_wasserstein_distance_hera():
    _basic_wasserstein(hera_wrap(1e-12), 1e-12, test_matching=False)
    _basic_wasserstein(hera_wrap(.1), .1, test_matching=False)
