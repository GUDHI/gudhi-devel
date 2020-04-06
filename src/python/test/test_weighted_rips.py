""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Yuichi Ike

    Copyright (C) 2020 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.weighted_rips_complex import WeightedRipsComplex
from gudhi.point_cloud.dtm import DTM
import numpy
from scipy.spatial.distance import cdist
import pytest

def test_dtm_rips_complex():
    pts = numpy.array([[2.0, 2], [0, 1], [3, 4]])
    dist = cdist(pts,pts)
    dtm = DTM(2, q=2, metric="precomputed")    
    r = dtm.fit_transform(dist)
    w_rips = WeightedRipsComplex(distance_mattix=dist, filtration_values=r)
    st = w_rips.create_simplex_tree(max_dimension=2)
    diag = st.persistence()
    assert diag == [(0, (1.5811388300841898, float("inf"))), (0, (1.5811388300841898, 2.699172818834085)), (0, (1.5811388300841898, 2.699172818834085))]
    
    