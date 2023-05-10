""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Marc Glisse

    Copyright (C) 2020 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.point_cloud.knn import KNearestNeighbors
import numpy as np
import pytest


def test_knn_explicit():
    base = np.array([[1.0, 1], [1, 2], [4, 2], [4, 3]])
    query = np.array([[1.0, 1], [2, 2], [4, 4]])
    knn = KNearestNeighbors(2, metric="manhattan", return_distance=True, return_index=True)
    knn.fit(base)
    r = knn.transform(query)
    assert r[0] == pytest.approx(np.array([[0, 1], [1, 0], [3, 2]]))
    assert r[1] == pytest.approx(np.array([[0.0, 1], [1, 2], [1, 2]]))

    knn = KNearestNeighbors(2, metric="chebyshev", return_distance=True, return_index=False)
    knn.fit(base)
    r = knn.transform(query)
    assert r == pytest.approx(np.array([[0.0, 1], [1, 1], [1, 2]]))
    r = (
        KNearestNeighbors(2, metric="chebyshev", return_distance=True, return_index=False, implementation="keops")
        .fit(base)
        .transform(query)
    )
    assert r == pytest.approx(np.array([[0.0, 1], [1, 1], [1, 2]]))
    r = (
        KNearestNeighbors(2, metric="chebyshev", return_distance=True, return_index=False, implementation="keops", enable_autodiff=True)
        .fit(base)
        .transform(query)
    )
    assert r == pytest.approx(np.array([[0.0, 1], [1, 1], [1, 2]]))

    knn = KNearestNeighbors(2, metric="minkowski", p=3, return_distance=False, return_index=True)
    knn.fit(base)
    r = knn.transform(query)
    assert np.array_equal(r, [[0, 1], [1, 0], [3, 2]])
    r = (
        KNearestNeighbors(2, metric="minkowski", p=3, return_distance=False, return_index=True, implementation="keops")
        .fit(base)
        .transform(query)
    )
    assert np.array_equal(r, [[0, 1], [1, 0], [3, 2]])

    dist = np.array([[0.0, 3, 8], [1, 0, 5], [1, 2, 0]])
    knn = KNearestNeighbors(2, metric="precomputed", return_index=True, return_distance=False)
    r = knn.fit_transform(dist)
    assert np.array_equal(r, [[0, 1], [1, 0], [2, 0]])
    knn = KNearestNeighbors(2, metric="precomputed", return_index=True, return_distance=True, sort_results=True)
    r = knn.fit_transform(dist)
    assert np.array_equal(r[0], [[0, 1], [1, 0], [2, 0]])
    assert np.array_equal(r[1], [[0, 3], [0, 1], [0, 1]])
    # Second time in parallel
    knn = KNearestNeighbors(2, metric="precomputed", return_index=True, return_distance=False, n_jobs=2, sort_results=True)
    r = knn.fit_transform(dist)
    assert np.array_equal(r, [[0, 1], [1, 0], [2, 0]])
    knn = KNearestNeighbors(2, metric="precomputed", return_index=True, return_distance=True, n_jobs=2)
    r = knn.fit_transform(dist)
    assert np.array_equal(r[0], [[0, 1], [1, 0], [2, 0]])
    assert np.array_equal(r[1], [[0, 3], [0, 1], [0, 1]])


def test_knn_compare():
    base = np.array([[1.0, 1], [1, 2], [4, 2], [4, 3]])
    query = np.array([[1.0, 1], [2, 2], [4, 4]])
    r0 = (
        KNearestNeighbors(2, implementation="ckdtree", return_index=True, return_distance=False)
        .fit(base)
        .transform(query)
    )
    r1 = (
        KNearestNeighbors(2, implementation="sklearn", return_index=True, return_distance=False)
        .fit(base)
        .transform(query)
    )
    r2 = (
        KNearestNeighbors(2, implementation="hnsw", return_index=True, return_distance=False).fit(base).transform(query)
    )
    r3 = (
        KNearestNeighbors(2, implementation="keops", return_index=True, return_distance=False)
        .fit(base)
        .transform(query)
    )
    assert np.array_equal(r0, r1) and np.array_equal(r0, r2) and np.array_equal(r0, r3)

    r0 = (
        KNearestNeighbors(2, implementation="ckdtree", return_index=True, return_distance=True)
        .fit(base)
        .transform(query)
    )
    r1 = (
        KNearestNeighbors(2, implementation="sklearn", return_index=True, return_distance=True)
        .fit(base)
        .transform(query)
    )
    r2 = KNearestNeighbors(2, implementation="hnsw", return_index=True, return_distance=True).fit(base).transform(query)
    r3 = (
        KNearestNeighbors(2, implementation="keops", return_index=True, return_distance=True).fit(base).transform(query)
    )
    assert np.array_equal(r0[0], r1[0]) and np.array_equal(r0[0], r2[0]) and np.array_equal(r0[0], r3[0])
    d0 = pytest.approx(r0[1])
    assert r1[1] == d0 and r2[1] == d0 and r3[1] == d0


def test_knn_nop():
    # This doesn't look super useful...
    p = np.array([[0.0]])
    assert None is KNearestNeighbors(
        k=1, return_index=False, return_distance=False, implementation="sklearn"
    ).fit_transform(p)
    assert None is KNearestNeighbors(
        k=1, return_index=False, return_distance=False, implementation="ckdtree"
    ).fit_transform(p)
    assert None is KNearestNeighbors(
        k=1, return_index=False, return_distance=False, implementation="hnsw", ef=5
    ).fit_transform(p)
    assert None is KNearestNeighbors(
        k=1, return_index=False, return_distance=False, implementation="keops"
    ).fit_transform(p)
    assert None is KNearestNeighbors(
        k=1, return_index=False, return_distance=False, metric="precomputed"
    ).fit_transform(p)
