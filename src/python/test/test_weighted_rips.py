""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Yuichi Ike and Masatoshi Takenouchi

    Copyright (C) 2020 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.weighted_rips_complex import WeightedRipsComplex
from gudhi.point_cloud.dtm import DTM
import numpy as np
from scipy.spatial.distance import cdist
import pytest

def test_dtm_rips_complex():
    pts = np.array([[2.0, 2], [0, 1], [3, 4]])
    dist = cdist(pts,pts)
    dtm = DTM(2, q=2, metric="precomputed")    
    r = dtm.fit_transform(dist)
    w_rips = WeightedRipsComplex(distance_mattix=dist, weights=r)
    st = w_rips.create_simplex_tree(max_dimension=2)
    persistence_intervals0 = st.persistence_intervals_in_dimension(0)
    assert persistence_intervals0 == pytest.approx(np.array([[1.58113883, 2.69917282],[1.58113883, 2.69917282], [1.58113883, float("inf")]]))
    
