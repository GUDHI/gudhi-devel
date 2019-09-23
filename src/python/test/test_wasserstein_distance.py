import gudhi
import numpy as np

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Theo Lacombe

    Copyright (C) 2019 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__author__ = "Theo Lacombe"
__copyright__ = "Copyright (C) 2019 Inria"
__license__ = "MIT"


def test_basic_wasserstein():
    diag1 = np.array([[2.7, 3.7], [9.6, 14.0], [34.2, 34.974]])
    diag2 = np.array([[2.8, 4.45], [9.5, 14.1]])
    diag3 = np.array([[0, 2], [4, 6]])
    diag4 = np.array([[0, 3], [4, 8]])
    emptydiag = np.array([[]])

    assert gudhi.wasserstein_distance(emptydiag, emptydiag, q=2., p=1.) == 0.
    assert gudhi.wasserstein_distance(emptydiag, emptydiag, q=np.inf, p=1.) == 0.
    assert gudhi.wasserstein_distance(emptydiag, emptydiag, q=np.inf, p=2.) == 0.
    assert gudhi.wasserstein_distance(emptydiag, emptydiag, q=2., p=2.) == 0.

    assert gudhi.wasserstein_distance(diag3, emptydiag, q=np.inf, p=1.) == 2.
    assert gudhi.wasserstein_distance(diag3, emptydiag, q=1., p=1.) == 4.

    assert gudhi.wasserstein_distance(diag4, emptydiag, q=1., p=2.) == 5.  # thank you Pythagorician triplets
    assert gudhi.wasserstein_distance(diag4, emptydiag, q=np.inf, p=2.) == 2.5
    assert gudhi.wasserstein_distance(diag4, emptydiag, q=2., p=2.) == 3.5355339059327378

    assert gudhi.wasserstein_distance(diag1, diag2, q=2., p=1.) == 1.4453593023967701
    assert gudhi.wasserstein_distance(diag1, diag2, q=2.35, p=1.74) == 0.9772734057168739

    assert gudhi.wasserstein_distance(diag1, emptydiag, q=2.35, p=1.7863) == 3.141592214572228

    assert gudhi.wasserstein_distance(diag3, diag4, q=1., p=1.) == 3.
    assert gudhi.wasserstein_distance(diag3, diag4, q=np.inf, p=1.) == 3.  # no diag matching here
    assert gudhi.wasserstein_distance(diag3, diag4, q=np.inf, p=2.) == np.sqrt(5)
    assert gudhi.wasserstein_distance(diag3, diag4, q=1., p=2.) == np.sqrt(5)
    assert gudhi.wasserstein_distance(diag3, diag4, q=4.5, p=2.) == np.sqrt(5)

    

