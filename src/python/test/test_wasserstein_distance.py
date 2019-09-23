import gudhi

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Theo Lacombe

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__author__ = "Theo Lacombe"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"


def test_basic_bottleneck():
    diag1 = np.array([[2.7, 3.7], [9.6, 14.0], [34.2, 34.974]])
    diag2 = np.array([[2.8, 4.45], [9.5, 14.1]])

    assert gudhi.wasserstein_distance(diag1, diag2) == 1.4453593023967701
