""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2024 Inria

    Modification(s):
      - 2024/07 Vincent Rouvreau: Separate Delaunay and Alpha tests in a different python file
      - YYYY/MM Author: Description of the modification
"""

from gudhi import AlphaComplex, DelaunayComplex, DelaunayCechComplex
import math
import numpy as np
import pytest
import random

try:
    # python3
    from itertools import zip_longest
except ImportError:
    # python2
    from itertools import izip_longest as zip_longest


def _empty_complex(simplicial_complex, precision):
    cplx = simplicial_complex(precision=precision)


def _one_2d_point_complex(simplicial_complex, precision):
    cplx = simplicial_complex(points=[[0, 0]], precision=precision)


def test_empty_complex():
    for simplicial_complex in [AlphaComplex, DelaunayComplex, DelaunayCechComplex]:
        for precision in ["fast", "safe", "exact"]:
            _empty_complex(simplicial_complex, precision)
            _one_2d_point_complex(simplicial_complex, precision)


def _infinite_threshold(simplicial_complex, precision):
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    cplx = simplicial_complex(points=point_list, precision=precision)

    simplex_tree = cplx.create_simplex_tree()
    assert simplex_tree._is_persistence_defined() == False

    assert simplex_tree.num_simplices() == 11
    assert simplex_tree.num_vertices() == 4

    assert list(simplex_tree.get_filtration()) == [
        ([0], 0.0),
        ([1], 0.0),
        ([2], 0.0),
        ([3], 0.0),
        ([0, 1], 0.25),
        ([0, 2], 0.25),
        ([1, 3], 0.25),
        ([2, 3], 0.25),
        ([1, 2], 0.5),
        ([0, 1, 2], 0.5),
        ([1, 2, 3], 0.5),
    ]

    assert simplex_tree.get_star([0]) == [
        ([0], 0.0),
        ([0, 1], 0.25),
        ([0, 1, 2], 0.5),
        ([0, 2], 0.25),
    ]
    assert simplex_tree.get_cofaces([0], 1) == [([0, 1], 0.25), ([0, 2], 0.25)]


def test_infinite_threshold():
    for simplicial_complex in [AlphaComplex, DelaunayCechComplex]:
        for precision in ["fast", "safe", "exact"]:
            _infinite_threshold(simplicial_complex, precision)


def _finite_threshold(simplicial_complex, precision):
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    cplx = simplicial_complex(points=point_list, precision=precision)

    simplex_tree = cplx.create_simplex_tree(max_alpha_square=0.25)

    assert simplex_tree.num_simplices() == 8
    assert simplex_tree.num_vertices() == 4

    assert list(simplex_tree.get_filtration()) == [
        ([0], 0.0),
        ([1], 0.0),
        ([2], 0.0),
        ([3], 0.0),
        ([0, 1], 0.25),
        ([0, 2], 0.25),
        ([1, 3], 0.25),
        ([2, 3], 0.25),
    ]
    assert simplex_tree.get_star([0]) == [([0], 0.0), ([0, 1], 0.25), ([0, 2], 0.25)]
    assert simplex_tree.get_cofaces([0], 1) == [([0, 1], 0.25), ([0, 2], 0.25)]


def test_filtered_complex():
    for simplicial_complex in [AlphaComplex, DelaunayCechComplex]:
        for precision in ["fast", "safe", "exact"]:
            _finite_threshold(simplicial_complex, precision)


def _safe_persistence_comparison(simplicial_complex, precision):
    # generate periodic signal
    time = np.arange(0, 10, 1)
    signal = [math.sin(x) for x in time]
    delta = math.pi
    delayed = [math.sin(x + delta) for x in time]

    # construct embedding
    embedding1 = [[signal[i], -signal[i]] for i in range(len(time))]
    embedding2 = [[signal[i], delayed[i]] for i in range(len(time))]

    # build simplicial_complex and simplex tree
    cplx1 = simplicial_complex(points=embedding1, precision=precision)
    simplex_tree1 = cplx1.create_simplex_tree()

    cplx2 = simplicial_complex(points=embedding2, precision=precision)
    simplex_tree2 = cplx2.create_simplex_tree()

    diag1 = simplex_tree1.persistence()
    diag2 = simplex_tree2.persistence()

    for first_p, second_p in zip_longest(diag1, diag2):
        assert first_p[0] == pytest.approx(second_p[0])
        assert first_p[1] == pytest.approx(second_p[1])


def test_safe_persistence_comparison():
    for simplicial_complex in [AlphaComplex, DelaunayCechComplex]:
        # Won't work for 'fast' version
        _safe_persistence_comparison(simplicial_complex, "safe")
        _safe_persistence_comparison(simplicial_complex, "exact")


def _delaunay_complex(precision):
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    cplx = DelaunayComplex(points=point_list, precision=precision)

    simplex_tree = cplx.create_simplex_tree()

    assert simplex_tree.num_simplices() == 11
    assert simplex_tree.num_vertices() == 4

    for filtered_value in simplex_tree.get_filtration():
        assert math.isnan(filtered_value[1])
    for filtered_value in simplex_tree.get_star([0]):
        assert math.isnan(filtered_value[1])
    for filtered_value in simplex_tree.get_cofaces([0], 1):
        assert math.isnan(filtered_value[1])


def test_delaunay_complex():
    for precision in ["fast", "safe", "exact"]:
        _delaunay_complex(precision)


def _3d_points_on_a_plane(simplicial_complex, precision):
    cplx = simplicial_complex(
        points=[
            [1.0, 1.0, 0.0],
            [7.0, 0.0, 0.0],
            [4.0, 6.0, 0.0],
            [9.0, 6.0, 0.0],
            [0.0, 14.0, 0.0],
            [2.0, 19.0, 0.0],
            [9.0, 17.0, 0.0],
        ],
        precision=precision,
    )

    simplex_tree = cplx.create_simplex_tree()
    assert simplex_tree.dimension() == 2
    assert simplex_tree.num_vertices() == 7
    assert simplex_tree.num_simplices() == 25


def test_3d_points_on_a_plane():
    for simplicial_complex in [AlphaComplex, DelaunayComplex, DelaunayCechComplex]:
        for precision in ["fast", "safe", "exact"]:
            _3d_points_on_a_plane(simplicial_complex, precision)


def _duplicated_2d_points_on_a_plane(simplicial_complex, precision):
    cplx = simplicial_complex(
        points=[
            [1.0, 1.0],
            [7.0, 0.0],  # This point is duplicate
            [4.0, 6.0],
            [9.0, 6.0],
            [0.0, 14.0],
            [2.0, 19.0],
            [7.0, 0.0],  # This point is duplicate
            [9.0, 17.0],
        ],
        precision=precision,
    )

    simplex_tree = cplx.create_simplex_tree()
    assert simplex_tree.dimension() == 2
    assert simplex_tree.num_vertices() == 7
    assert simplex_tree.num_simplices() == 25


def test_duplicated_2d_points_on_a_plane():
    for simplicial_complex in [AlphaComplex, DelaunayComplex, DelaunayCechComplex]:
        for precision in ["fast", "safe", "exact"]:
            _duplicated_2d_points_on_a_plane(simplicial_complex, precision)

def test_weighted_exceptions():
    points=[[1, 1], [7, 0], [4, 6], [9, 6], [0, 14], [2, 19], [9, 17]]
    nb_wgts = len(points)
    while nb_wgts == len(points):
      nb_wgts = random.randint(2, 2*len(points))
    weights = nb_wgts * [1.]
    print(f"nb points = {len(points)} vs nb weights = {len(weights)}")
    with pytest.raises(ValueError):
        # When number of points does not correspond to the number of weights
        dc = DelaunayComplex(points=points, weights=weights)

    weights = len(points) * [1.]
    dc = DelaunayComplex(points=points, weights=weights)
    # No Weighted Delaunay-Cech available
    with pytest.raises(ValueError):
        dc.create_simplex_tree(filtration='cech')

    # Weighted Alpha complex cannot output square root values
    with pytest.raises(ValueError):
        dc.create_simplex_tree(filtration='alpha', output_squared_values=False)
