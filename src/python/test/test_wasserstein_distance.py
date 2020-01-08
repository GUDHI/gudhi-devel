from gudhi.wasserstein import wasserstein_distance
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

    assert wasserstein_distance(emptydiag, emptydiag, internal_p=2., q=1.) == 0.
    assert wasserstein_distance(emptydiag, emptydiag, internal_p=np.inf, q=1.) == 0.
    assert wasserstein_distance(emptydiag, emptydiag, internal_p=np.inf, q=2.) == 0.
    assert wasserstein_distance(emptydiag, emptydiag, internal_p=2., q=2.) == 0.

    assert wasserstein_distance(diag3, emptydiag, internal_p=np.inf, q=1.) == 2.
    assert wasserstein_distance(diag3, emptydiag, internal_p=1., q=1.) == 4.

    assert wasserstein_distance(diag4, emptydiag, internal_p=1., q=2.) == 5.  # thank you Pythagorician triplets
    assert wasserstein_distance(diag4, emptydiag, internal_p=np.inf, q=2.) == 2.5
    assert wasserstein_distance(diag4, emptydiag, internal_p=2., q=2.) == 3.5355339059327378

    assert wasserstein_distance(diag1, diag2, internal_p=2., q=1.) == 1.4453593023967701
    assert wasserstein_distance(diag1, diag2, internal_p=2.35, q=1.74) == 0.9772734057168739

    assert wasserstein_distance(diag1, emptydiag, internal_p=2.35, q=1.7863) == 3.141592214572228

    assert wasserstein_distance(diag3, diag4, internal_p=1., q=1.) == 3.
    assert wasserstein_distance(diag3, diag4, internal_p=np.inf, q=1.) == 3.  # no diag matching here
    assert wasserstein_distance(diag3, diag4, internal_p=np.inf, q=2.) == np.sqrt(5)
    assert wasserstein_distance(diag3, diag4, internal_p=1., q=2.) == np.sqrt(5)
    assert wasserstein_distance(diag3, diag4, internal_p=4.5, q=2.) == np.sqrt(5)

