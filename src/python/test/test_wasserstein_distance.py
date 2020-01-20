""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Theo Lacombe

    Copyright (C) 2019 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.wasserstein import wasserstein_distance as pot
from gudhi.hera import wasserstein_distance as hera
import numpy as np

__author__ = "Theo Lacombe"
__copyright__ = "Copyright (C) 2019 Inria"
__license__ = "MIT"

def _basic_wasserstein(wasserstein_distance, delta, test_infinity=True):
    diag1 = np.array([[2.7, 3.7], [9.6, 14.0], [34.2, 34.974]])
    diag2 = np.array([[2.8, 4.45], [9.5, 14.1]])
    diag3 = np.array([[0, 2], [4, 6]])
    diag4 = np.array([[0, 3], [4, 8]])
    emptydiag = np.array([])

    # We just need to handle positive numbers here
    def approx(a, b):
        f = 1 + delta
        return a <= b*f and b <= a*f

    assert wasserstein_distance(emptydiag, emptydiag, internal_p=2.,     order=1.) == 0.
    assert wasserstein_distance(emptydiag, emptydiag, internal_p=np.inf, order=1.) == 0.
    assert wasserstein_distance(emptydiag, emptydiag, internal_p=np.inf, order=2.) == 0.
    assert wasserstein_distance(emptydiag, emptydiag, internal_p=2.,     order=2.) == 0.

    assert approx(wasserstein_distance(diag3, emptydiag, internal_p=np.inf,     order=1.), 2.)
    assert approx(wasserstein_distance(diag3, emptydiag, internal_p=1.,         order=1.), 4.)

    assert approx(wasserstein_distance(diag4, emptydiag, internal_p=1.,     order=2.), 5.)  # thank you Pythagorician triplets
    assert approx(wasserstein_distance(diag4, emptydiag, internal_p=np.inf, order=2.), 2.5)
    assert approx(wasserstein_distance(diag4, emptydiag, internal_p=2.,     order=2.), 3.5355339059327378)

    assert approx(wasserstein_distance(diag1, diag2, internal_p=2.,   order=1.)  , 1.4453593023967701)
    assert approx(wasserstein_distance(diag1, diag2, internal_p=2.35, order=1.74), 0.9772734057168739)

    assert approx(wasserstein_distance(diag1, emptydiag, internal_p=2.35, order=1.7863), 3.141592214572228)

    assert approx(wasserstein_distance(diag3, diag4, internal_p=1.,     order=1.), 3.)
    assert approx(wasserstein_distance(diag3, diag4, internal_p=np.inf, order=1.), 3.)  # no diag matching here
    assert approx(wasserstein_distance(diag3, diag4, internal_p=np.inf, order=2.), np.sqrt(5))
    assert approx(wasserstein_distance(diag3, diag4, internal_p=1.,     order=2.), np.sqrt(5))
    assert approx(wasserstein_distance(diag3, diag4, internal_p=4.5,    order=2.), np.sqrt(5))

    if(not test_infinity):
        return

    diag5 = np.array([[0, 3], [4, np.inf]])
    diag6 = np.array([[7, 8], [4, 6], [3, np.inf]])

    assert wasserstein_distance(diag4, diag5) == np.inf
    assert approx(wasserstein_distance(diag5, diag6, order=1, internal_p=np.inf), 4.)

def hera_wrap(delta):
    def fun(*kargs,**kwargs):
        return hera(*kargs,**kwargs,delta=delta)
    return fun

def test_wasserstein_distance_pot():
    _basic_wasserstein(pot, 1e-15, False)

def test_wasserstein_distance_hera():
    _basic_wasserstein(hera_wrap(1e-12), 1e-12)
    _basic_wasserstein(hera_wrap(.1), .1)
