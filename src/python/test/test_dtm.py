""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Marc Glisse

    Copyright (C) 2020 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.point_cloud.dtm import DTM
import numpy


def test_dtm_euclidean():
    pts = numpy.random.rand(1000,4)
    k = 3
    dtm = DTM(k,implementation="ckdtree")
    print(dtm.fit_transform(pts))
    dtm = DTM(k,implementation="sklearn")
    print(dtm.fit_transform(pts))
    dtm = DTM(k,implementation="sklearn",algorithm="brute")
    print(dtm.fit_transform(pts))
    dtm = DTM(k,implementation="hnsw")
    print(dtm.fit_transform(pts))
    from scipy.spatial.distance import cdist
    d = cdist(pts,pts)
    dtm = DTM(k,metric="precomputed")
    print(dtm.fit_transform(d))
    dtm = DTM(k,implementation="keops")
    print(dtm.fit_transform(pts))

