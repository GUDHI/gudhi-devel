""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

import gudhi
import gudhi.hera
import pytest

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"


def test_basic_bottleneck():
    diag1 = [[2.7, 3.7], [9.6, 14.0], [34.2, 34.974], [3.0, float("Inf")]]
    diag2 = [[2.8, 4.45], [9.5, 14.1], [3.2, float("Inf")]]

    assert gudhi.bottleneck_distance(diag1, diag2) == 0.75
    assert gudhi.bottleneck_distance(diag1, diag2, 0.1) == pytest.approx(0.75, abs=0.1)
    assert gudhi.hera.bottleneck_distance(diag1, diag2, 0) == 0.75
    assert gudhi.hera.bottleneck_distance(diag1, diag2, 0.1) == pytest.approx(0.75, rel=0.1)

    import numpy as np

    # Translating both diagrams along the diagonal should not affect the distance,
    # checks that negative numbers are not an issue
    diag1 = np.array(diag1) - 100
    diag2 = np.array(diag2) - 100

    assert gudhi.bottleneck_distance(diag1, diag2) == 0.75
    assert gudhi.bottleneck_distance(diag1, diag2, 0.1) == pytest.approx(0.75, abs=0.1)
    assert gudhi.hera.bottleneck_distance(diag1, diag2, 0) == 0.75
    assert gudhi.hera.bottleneck_distance(diag1, diag2, 0.1) == pytest.approx(0.75, rel=0.1)
