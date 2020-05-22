""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Yuichi Ike

    Copyright (C) 2020 Inria, Copyright (C) 2020 FUjitsu Laboratories Ltd.

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.dtm_rips_complex import DTMRipsComplex
import numpy as np
from math import sqrt
import pytest

def test_dtm_rips_complex():
    pts = np.array([[2.0, 2.0], [0.0, 1.0], [3.0, 4.0]])
    dtm_rips = DTMRipsComplex(points=pts, k=2)
    st = dtm_rips.create_simplex_tree(max_dimension=2)
    st.persistence()
    persistence_intervals0 = st.persistence_intervals_in_dimension(0)
    assert persistence_intervals0 == pytest.approx(np.array([[3.16227766, 5.39834564],[3.16227766, 5.39834564], [3.16227766, float("inf")]]))
    
def test_compatibility_with_rips():
    distance_matrix = np.array([[0, 1, 1, sqrt(2)], [1, 0, sqrt(2), 1], [1, sqrt(2), 0, 1], [sqrt(2), 1, 1, 0]])
    dtm_rips = DTMRipsComplex(distance_matrix=distance_matrix, max_filtration=42)
    st = dtm_rips.create_simplex_tree(max_dimension=1)
    assert list(st.get_filtration()) == [
        ([0], 0.0),
        ([1], 0.0),
        ([2], 0.0),
        ([3], 0.0),
        ([0, 1], 1.0),
        ([0, 2], 1.0),
        ([1, 3], 1.0),
        ([2, 3], 1.0),
        ([1, 2], sqrt(2)),
        ([0, 3], sqrt(2)),
    ]

