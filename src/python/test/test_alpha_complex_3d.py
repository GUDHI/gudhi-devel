""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2021 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi import AlphaComplex3D
import pytest
import numpy as np

try:
    # python3
    from itertools import zip_longest
except ImportError:
    # python2
    from itertools import izip_longest as zip_longest



def _empty_alpha(precision):
    alpha_complex = AlphaComplex3D(precision = precision)
    assert alpha_complex.__is_defined() == True

def _one_3d_point_alpha(precision):
    alpha_complex = AlphaComplex3D(points=[[0, 0, 0]], precision = precision)
    assert alpha_complex.__is_defined() == True

def test_empty_alpha():
    for precision in ['fast', 'safe', 'exact']:
        _empty_alpha(precision)
        _one_3d_point_alpha(precision)

def _infinite_alpha(precision):
    point_list = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]]
    alpha_complex = AlphaComplex3D(points=point_list, precision = precision)
    assert alpha_complex.__is_defined() == True

    stree = alpha_complex.create_simplex_tree()
    assert stree.__is_persistence_defined() == False

    assert stree.num_simplices() == 51
    assert stree.num_vertices() == len(point_list)

    for filtration in stree.get_filtration():
        if len(filtration[0]) == 1:
            assert filtration[1] == 0.
        if len(filtration[0]) == 4:
            assert filtration[1] == 0.75

    for idx in range(len(point_list)):
        pt_idx = point_list.index(alpha_complex.get_point(idx))
        assert pt_idx >= 0
        assert pt_idx < len(point_list)

    with pytest.raises(IndexError):
        alpha_complex.get_point(len(point_list))

def test_infinite_alpha():
    for precision in ['fast', 'safe', 'exact']:
        _infinite_alpha(precision)

def _filtered_alpha(precision):
    point_list = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]]
    filtered_alpha = AlphaComplex3D(points=point_list, precision = precision)

    stree = filtered_alpha.create_simplex_tree(max_alpha_square=0.25)

    assert stree.num_simplices() == 20
    assert stree.num_vertices() == len(point_list)

    for filtration in stree.get_filtration():
        if len(filtration[0]) == 1:
            assert filtration[1] == 0.
        elif len(filtration[0]) == 2:
            assert filtration[1] == 0.25
        else:
            assert False

    for idx in range(len(point_list)):
        pt_idx = point_list.index(filtered_alpha.get_point(idx))
        assert pt_idx >= 0
        assert pt_idx < len(point_list)

    with pytest.raises(IndexError):
        filtered_alpha.get_point(len(point_list))

def test_filtered_alpha():
    for precision in ['fast', 'safe', 'exact']:
        _filtered_alpha(precision)

def _3d_points_on_a_plane(precision):
    alpha = AlphaComplex3D(points = [[1.0, 1.0 , 0.0],
                                     [7.0, 0.0 , 0.0],
                                     [4.0, 6.0 , 0.0],
                                     [9.0, 6.0 , 0.0],
                                     [0.0, 14.0, 0.0],
                                     [2.0, 19.0, 0.0],
                                     [9.0, 17.0, 0.0]], precision = precision)

    with pytest.raises(ValueError):
        stree = alpha.create_simplex_tree()

def test_3d_points_on_a_plane():
    for precision in ['fast', 'safe', 'exact']:
        _3d_points_on_a_plane(precision)

def test_inconsistency_points_and_weights():
    points = [[1.0, 1.0 , 1.0],
              [7.0, 0.0 , 2.0],
              [4.0, 6.0 , 0.0],
              [9.0, 6.0 , 1.0],
              [0.0, 14.0, 2.0],
              [2.0, 19.0, 0.0],
              [9.0, 17.0, 1.0]]
    with pytest.raises(ValueError):
        # 7 points, 8 weights, on purpose
        alpha = AlphaComplex3D(points = points,
                                weights = [1., 2., 3., 4., 5., 6., 7., 8.])

    with pytest.raises(ValueError):
        # 7 points, 6 weights, on purpose
        alpha = AlphaComplex3D(points = points,
                                weights = [1., 2., 3., 4., 5., 6.])

def _weighted_doc_example(precision):
    pts = [[ 1., -1., -1.],
           [-1.,  1., -1.],
           [-1., -1.,  1.],
           [ 1.,  1.,  1.],
           [ 2.,  2.,  2.]]
    wgts = [4., 4., 4., 4., 1.]
    alpha = AlphaComplex3D(points = pts, weights = wgts, precision = precision)
    stree = alpha.create_simplex_tree()

    # Needs to retrieve points as points are shuffled
    get_idx = lambda idx: pts.index(alpha.get_point(idx))
    indices = [get_idx(x) for x in range(len(pts))]

    assert stree.filtration([indices[x] for x in [0, 1, 2, 3]]) == pytest.approx(-1.)
    assert stree.filtration([indices[x] for x in [0, 1, 3, 4]]) == pytest.approx(95.)
    assert stree.filtration([indices[x] for x in [0, 2, 3, 4]]) == pytest.approx(95.)
    assert stree.filtration([indices[x] for x in [1, 2, 3, 4]]) == pytest.approx(95.)

def test_weighted_doc_example():
    for precision in ['fast', 'safe', 'exact']:
        _weighted_doc_example(precision)

def test_points_not_in_3d():
    with pytest.raises(ValueError):
        alpha = AlphaComplex3D(points = np.random.rand(6,2))
    with pytest.raises(ValueError):
        alpha = AlphaComplex3D(points = np.random.rand(6,4))

    alpha = AlphaComplex3D(points = np.random.rand(6,3))
    assert alpha.__is_defined() == True