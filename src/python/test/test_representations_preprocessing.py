""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2021 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.representations.preprocessing import DimensionSelector
import numpy as np
import pytest

H0_0 = np.array([0.0, 0.0])
H1_0 = np.array([1.0, 0.0])
H0_1 = np.array([0.0, 1.0])
H1_1 = np.array([1.0, 1.0])
H0_2 = np.array([0.0, 2.0])
H1_2 = np.array([1.0, 2.0])


def test_dimension_selector():
    X = [[H0_0, H1_0], [H0_1, H1_1], [H0_2, H1_2]]
    ds = DimensionSelector(index=0)
    h0 = ds.fit_transform(X)
    np.testing.assert_array_equal(h0[0], H0_0)
    np.testing.assert_array_equal(h0[1], H0_1)
    np.testing.assert_array_equal(h0[2], H0_2)

    ds = DimensionSelector(index=1)
    h1 = ds.fit_transform(X)
    np.testing.assert_array_equal(h1[0], H1_0)
    np.testing.assert_array_equal(h1[1], H1_1)
    np.testing.assert_array_equal(h1[2], H1_2)

    ds = DimensionSelector(index=2)
    with pytest.raises(IndexError):
        h2 = ds.fit_transform([[H0_0, H1_0], [H0_1, H1_1], [H0_2, H1_2]])
