""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Marc Glisse

    Copyright (C) 2020 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.point_cloud.dtm import DistanceToMeasure
import numpy
import pytest
import torch


def test_dtm_compare_euclidean():
    pts = numpy.random.rand(1000, 4)
    k = 3
    dtm = DistanceToMeasure(k, implementation="ckdtree")
    r0 = dtm.fit_transform(pts)
    dtm = DistanceToMeasure(k, implementation="sklearn")
    r1 = dtm.fit_transform(pts)
    assert r1 == pytest.approx(r0)
    dtm = DistanceToMeasure(k, implementation="sklearn", algorithm="brute")
    r2 = dtm.fit_transform(pts)
    assert r2 == pytest.approx(r0)
    dtm = DistanceToMeasure(k, implementation="hnsw")
    r3 = dtm.fit_transform(pts)
    assert r3 == pytest.approx(r0)
    from scipy.spatial.distance import cdist

    d = cdist(pts, pts)
    dtm = DistanceToMeasure(k, metric="precomputed")
    r4 = dtm.fit_transform(d)
    assert r4 == pytest.approx(r0)
    dtm = DistanceToMeasure(k, metric="precomputed", n_jobs=2)
    r4b = dtm.fit_transform(d)
    assert r4b == pytest.approx(r0)
    dtm = DistanceToMeasure(k, implementation="keops")
    r5 = dtm.fit_transform(pts)
    assert r5 == pytest.approx(r0)
    pts2 = torch.tensor(pts, requires_grad=True)
    assert pts2.grad is None
    dtm = DistanceToMeasure(k, implementation="keops", enable_autodiff=True)
    r6 = dtm.fit_transform(pts2)
    assert r6.detach().numpy() == pytest.approx(r0)
    r6.sum().backward()
    assert pts2.grad is not None


def test_dtm_precomputed():
    dist = numpy.array([[1.0, 3, 8], [1, 5, 5], [0, 2, 3]])
    dtm = DistanceToMeasure(2, q=1, metric="neighbors")
    r = dtm.fit_transform(dist)
    assert r == pytest.approx([2.0, 3, 1])

    dist = numpy.array([[2.0, 2], [0, 1], [3, 4]])
    dtm = DistanceToMeasure(2, q=2, metric="neighbors")
    r = dtm.fit_transform(dist)
    assert r == pytest.approx([2.0, 0.707, 3.5355], rel=0.01)
