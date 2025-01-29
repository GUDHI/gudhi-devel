""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2024 Inria

    Modification(s):
      - 2024/07 Vincent Rouvreau: Separate Delaunay and Alpha tests in a different python file
      - YYYY/MM Author: Description of the modification
"""

from gudhi import weighted_alpha_complex, alpha_complex, delaunay_complex, delaunay_cech_complex
import math
import numpy as np
import pytest

try:
    # python3
    from itertools import zip_longest
except ImportError:
    # python2
    from itertools import izip_longest as zip_longest



def test_empty_complex():
    # Specific for delaunay_complex as no precision
    stree = delaunay_complex()
    assert stree.is_empty()

    stree = delaunay_complex(points=[[0, 0]])
    assert stree.num_vertices() == 1
    assert stree.num_simplices() == 1

    for precision in ['fast', 'safe', 'exact']:
        # Specific for delaunay_complex as it requires weights
        stree = weighted_alpha_complex(precision = precision)
        assert stree.is_empty()

        stree = weighted_alpha_complex(points=[[0, 0]], weights=[0], precision = precision)
        assert stree.num_vertices() == 1
        assert stree.num_simplices() == 1

        for simplicial_complex_helper in [alpha_complex, delaunay_cech_complex]:
            stree = simplicial_complex_helper(precision = precision)
            assert stree.is_empty()

            stree = simplicial_complex_helper(points=[[0, 0]], precision = precision)
            assert stree.num_vertices() == 1
            assert stree.num_simplices() == 1


def _infinite_threshold(simplicial_complex_helper, precision):
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    simplex_tree = simplicial_complex_helper(points=point_list, precision=precision)

    assert simplex_tree._is_persistence_defined() == False

    assert simplex_tree.num_simplices() == 11
    assert simplex_tree.num_vertices() == 4

    diag_filt = 1. / math.sqrt(2.)
    simplices = [filt[0] for filt in simplex_tree.get_filtration()]
    assert simplices ==  [
        [0],
        [1],
        [2],
        [3],
        [0, 1],
        [0, 2],
        [1, 3],
        [2, 3],
        [1, 2],
        [0, 1, 2],
        [1, 2, 3],
    ]
    filtrations = [filt[1] for filt in simplex_tree.get_filtration()]
    np.testing.assert_array_almost_equal(filtrations, [
        0.0,
        0.0,
        0.0,
        0.0,
        0.5,
        0.5,
        0.5,
        0.5,
        diag_filt,
        diag_filt,
        diag_filt,
    ])

    simplices = [filt[0] for filt in simplex_tree.get_star([0])]
    assert simplices == [
        [0],
        [0, 1],
        [0, 1, 2],
        [0, 2],
    ]
    filtrations = [filt[1] for filt in simplex_tree.get_star([0])]
    np.testing.assert_array_almost_equal(filtrations, [
        0.0,
        0.5,
        diag_filt,
        0.5,
    ])

    assert simplex_tree.get_cofaces([0], 1) == [([0, 1], 0.5), ([0, 2], 0.5)]

def test_infinite_threshold():
    for simplicial_complex_helper in [alpha_complex, delaunay_cech_complex]:
        for precision in ['fast', 'safe', 'exact']:
            _infinite_threshold(simplicial_complex_helper, precision)

def _finite_threshold(simplicial_complex_helper, precision):
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    simplex_tree = simplicial_complex_helper(points=point_list, max_alpha=0.5, precision=precision)

    assert simplex_tree.num_simplices() == 8
    assert simplex_tree.num_vertices() == 4

    assert list(simplex_tree.get_filtration()) == [
        ([0], 0.0),
        ([1], 0.0),
        ([2], 0.0),
        ([3], 0.0),
        ([0, 1], 0.5),
        ([0, 2], 0.5),
        ([1, 3], 0.5),
        ([2, 3], 0.5),
    ]
    assert simplex_tree.get_star([0]) == [([0], 0.0), ([0, 1], 0.5), ([0, 2], 0.5)]
    assert simplex_tree.get_cofaces([0], 1) == [([0, 1], 0.5), ([0, 2], 0.5)]

def test_filtered_complex():
    for simplicial_complex_helper in [alpha_complex, delaunay_cech_complex]:
        for precision in ['fast', 'safe', 'exact']:
            _finite_threshold(simplicial_complex_helper, precision)

def _safe_persistence_comparison(simplicial_complex_helper, precision):
    #generate periodic signal
    time = np.arange(0, 10, 1)
    signal = [math.sin(x) for x in time]
    delta = math.pi
    delayed = [math.sin(x + delta) for x in time]

    #construct embedding
    embedding1 = [[signal[i], -signal[i]] for i in range(len(time))]
    embedding2 = [[signal[i], delayed[i]] for i in range(len(time))]

    #build simplicial_complex_helper and simplex tree
    simplex_tree1 = simplicial_complex_helper(points=embedding1, precision = precision)

    simplex_tree2 = simplicial_complex_helper(points=embedding2, precision = precision)

    diag1 = simplex_tree1.persistence()
    diag2 = simplex_tree2.persistence()

    for (first_p, second_p) in zip_longest(diag1, diag2):
        assert first_p[0] == pytest.approx(second_p[0])
        assert first_p[1] == pytest.approx(second_p[1])


def test_safe_persistence_comparison():
    for simplicial_complex_helper in [alpha_complex, delaunay_cech_complex]:
        # Won't work for 'fast' version
        _safe_persistence_comparison(simplicial_complex_helper, 'safe')
        _safe_persistence_comparison(simplicial_complex_helper, 'exact')

def test_delaunay_complex():
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    simplex_tree = delaunay_complex(points=point_list)

    assert simplex_tree.num_simplices() == 11
    assert simplex_tree.num_vertices() == 4

    for filtered_value in simplex_tree.get_filtration():
        assert math.isnan(filtered_value[1])
    for filtered_value in simplex_tree.get_star([0]):
        assert math.isnan(filtered_value[1])
    for filtered_value in simplex_tree.get_cofaces([0], 1):
        assert math.isnan(filtered_value[1])

def _3d_points_on_a_plane(simplicial_complex_helper):
    simplex_tree = simplicial_complex_helper(points = [[1.0, 1.0 , 0.0],
                                                       [7.0, 0.0 , 0.0],
                                                       [4.0, 6.0 , 0.0],
                                                       [9.0, 6.0 , 0.0],
                                                       [0.0, 14.0, 0.0],
                                                       [2.0, 19.0, 0.0],
                                                       [9.0, 17.0, 0.0]])

    assert simplex_tree.dimension() == 2
    assert simplex_tree.num_vertices() == 7
    assert simplex_tree.num_simplices() == 25

def test_3d_points_on_a_plane():
    for simplicial_complex_helper in [alpha_complex, delaunay_complex, delaunay_cech_complex]:
        for precision in ['fast', 'safe', 'exact']:
            _3d_points_on_a_plane(simplicial_complex_helper)

def _duplicated_2d_points_on_a_plane(simplicial_complex_helper):
    simplex_tree = simplicial_complex_helper(points = [[1.0, 1.0 ],
                                                       [7.0, 0.0 ], # This point is duplicate
                                                       [4.0, 6.0 ],
                                                       [9.0, 6.0 ],
                                                       [0.0, 14.0],
                                                       [2.0, 19.0],
                                                       [7.0, 0.0 ], # This point is duplicate
                                                       [9.0, 17.0]])

    assert simplex_tree.dimension() == 2
    assert simplex_tree.num_vertices() == 7
    assert simplex_tree.num_simplices() == 25

def test_duplicated_2d_points_on_a_plane():
    for simplicial_complex_helper in [alpha_complex, delaunay_complex, delaunay_cech_complex]:
        _duplicated_2d_points_on_a_plane(simplicial_complex_helper)

def test_output_squared_values():
    for simplicial_complex_helper in [alpha_complex, delaunay_cech_complex]:
        for precision in ['fast', 'safe', 'exact']:
            for max_alpha in [float('inf'), math.sqrt(20.)]:
                pts=[[1, 1], [7, 0], [4, 6], [9, 6], [0, 14], [2, 19], [9, 17]]
                stree = simplicial_complex_helper(points=pts, precision=precision,
                                                  output_squared_values=True, max_alpha=max_alpha)
                stree_sqrt = simplicial_complex_helper(points=pts, precision=precision,
                                                       output_squared_values=False, max_alpha=max_alpha)
                assert stree.num_simplices() == stree_sqrt.num_simplices()
                for simplex, filt in stree_sqrt.get_filtration():
                    # np.testing.assert_almost_equal(float('nan'), float('nan')) is ok
                    # while float('nan') == float('nan') is False
                    np.testing.assert_almost_equal(filt, math.sqrt(stree.filtration(simplex)))
