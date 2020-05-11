""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Yuichi Ike and Masatoshi Takenouchi

    Copyright (C) 2020 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.weighted_rips_complex import WeightedRipsComplex
from gudhi.point_cloud.dtm import DistanceToMeasure
import numpy as np
from math import sqrt
from scipy.spatial.distance import cdist
import pytest

def test_non_dtm_rips_complex():
    dist = [[], [1]]
    weights = [1, 100]
    w_rips = WeightedRipsComplex(distance_matrix=dist, weights=weights)
    st = w_rips.create_simplex_tree(max_dimension=2)
    assert st.filtration([0,1]) == pytest.approx(200.0)

def test_compatibility_with_rips():
    distance_matrix = [[0], [1, 0], [1, sqrt(2), 0], [sqrt(2), 1, 1, 0]]
    w_rips = WeightedRipsComplex(distance_matrix=distance_matrix,max_filtration=42)
    st = w_rips.create_simplex_tree(max_dimension=1)
    assert list(st.get_filtration()) == [
        ([0], 0.0),
        ([1], 0.0),
        ([2], 0.0),
        ([3], 0.0),
        ([0, 1], 1.0),
        ([0, 2], 1.0),
        ([1, 3], 1.0),
        ([2, 3], 1.0),
        ([1, 2], 1.4142135623730951),
        ([0, 3], 1.4142135623730951),
    ]

def test_compatibility_with_filtered_rips():
    distance_matrix = [[0], [1, 0], [1, sqrt(2), 0], [sqrt(2), 1, 1, 0]]
    w_rips = WeightedRipsComplex(distance_matrix=distance_matrix,max_filtration=1.0)
    st = w_rips.create_simplex_tree(max_dimension=1)
    
    assert st.__is_defined() == True
    assert st.__is_persistence_defined() == False

    assert st.num_simplices() == 8
    assert st.num_vertices() == 4

def test_dtm_rips_complex():
    pts = np.array([[2.0, 2.0], [0.0, 1.0], [3.0, 4.0]])
    dist = cdist(pts,pts)
    dtm = DistanceToMeasure(2, q=2, metric="precomputed")
    r = dtm.fit_transform(dist)
    w_rips = WeightedRipsComplex(distance_matrix=dist, weights=r)
    st = w_rips.create_simplex_tree(max_dimension=2)
    st.persistence()
    persistence_intervals0 = st.persistence_intervals_in_dimension(0)
    assert persistence_intervals0 == pytest.approx(np.array([[3.16227766, 5.39834564],[3.16227766, 5.39834564], [3.16227766, float("inf")]]))
    
