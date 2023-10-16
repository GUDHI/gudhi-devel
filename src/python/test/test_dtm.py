""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Marc Glisse

    Copyright (C) 2020 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.point_cloud.dtm import DistanceToMeasure, DTMDensity
import numpy
import pytest
import torch
import math
import warnings


def test_dtm_compare_euclidean():
    pts = numpy.random.rand(1000, 4)
    k = 6
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
    assert r3 == pytest.approx(r0, rel=0.1)
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
    assert not torch.isnan(pts2.grad).any()
    pts2 = torch.tensor(pts, requires_grad=True)
    assert pts2.grad is None
    dtm = DistanceToMeasure(k, implementation="ckdtree", enable_autodiff=True)
    r7 = dtm.fit_transform(pts2)
    assert r7.detach().numpy() == pytest.approx(r0)
    r7.sum().backward()
    assert not torch.isnan(pts2.grad).any()


def test_dtm_precomputed():
    dist = numpy.array([[1.0, 3, 8], [1, 5, 5], [0, 2, 3]])
    dtm = DistanceToMeasure(2, q=1, metric="neighbors")
    r = dtm.fit_transform(dist)
    assert r == pytest.approx([2.0, 3, 1])

    dist = numpy.array([[2.0, 2], [0, 1], [3, 4]])
    dtm = DistanceToMeasure(2, q=2, metric="neighbors")
    r = dtm.fit_transform(dist)
    assert r == pytest.approx([2.0, 0.707, 3.5355], rel=0.01)


def test_density_normalized():
    sample = numpy.random.normal(0, 1, (1000000, 2))
    queries = numpy.array([[0.0, 0.0], [-0.5, 0.7], [0.4, 1.7]])
    expected = numpy.exp(-(queries ** 2).sum(-1) / 2) / (2 * math.pi)
    estimated = DTMDensity(k=150, normalize=True).fit(sample).transform(queries)
    assert estimated == pytest.approx(expected, rel=0.4)


def test_density():
    distances = [[0, 1, 10], [2, 0, 30], [1, 3, 5]]
    density = DTMDensity(k=2, metric="neighbors", dim=1).fit_transform(distances)
    expected = numpy.array([2.0, 1.0, 0.5])
    assert density == pytest.approx(expected)
    distances = [[0, 1], [2, 0], [1, 3]]
    density = DTMDensity(metric="neighbors", dim=1).fit_transform(distances)
    assert density == pytest.approx(expected)
    density = DTMDensity(weights=[0.5, 0.5], metric="neighbors", dim=1).fit_transform(distances)
    assert density == pytest.approx(expected)

def test_dtm_overflow_warnings():
    pts = numpy.array([[10., 100000000000000000000000000000.], [1000., 100000000000000000000000000.]])
    impl_warn = ["keops", "hnsw"]
    for impl in impl_warn:
        with warnings.catch_warnings(record=True) as w:
            dtm = DistanceToMeasure(2, implementation=impl)
            r = dtm.fit_transform(pts)
            assert len(w) == 1
            assert issubclass(w[0].category, RuntimeWarning)
            assert "Overflow" in str(w[0].message)
