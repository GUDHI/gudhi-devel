""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

import gudhi
import math
import numpy as np
import pytest
try:
    # python3
    from itertools import zip_longest
except ImportError:
    # python2
    from itertools import izip_longest as zip_longest

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"


def _empty_alpha(precision):
    alpha_complex = gudhi.AlphaComplex(points=[[0, 0]], precision = precision)
    assert alpha_complex.__is_defined() == True

def test_empty_alpha():
    for precision in ['fast', 'safe', 'exact']:
        _empty_alpha(precision)

def _infinite_alpha(precision):
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    alpha_complex = gudhi.AlphaComplex(points=point_list, precision = precision)
    assert alpha_complex.__is_defined() == True

    simplex_tree = alpha_complex.create_simplex_tree()
    assert simplex_tree.__is_persistence_defined() == False

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

    assert point_list[0] == alpha_complex.get_point(0)
    assert point_list[1] == alpha_complex.get_point(1)
    assert point_list[2] == alpha_complex.get_point(2)
    assert point_list[3] == alpha_complex.get_point(3)
    try:
        alpha_complex.get_point(4) == []
    except IndexError:
        pass
    else:
        assert False
    try:
        alpha_complex.get_point(125) == []
    except IndexError:
        pass
    else:
        assert False

def test_infinite_alpha():
    for precision in ['fast', 'safe', 'exact']:
        _infinite_alpha(precision)

def _filtered_alpha(precision):
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    filtered_alpha = gudhi.AlphaComplex(points=point_list, precision = precision)

    simplex_tree = filtered_alpha.create_simplex_tree(max_alpha_square=0.25)

    assert simplex_tree.num_simplices() == 8
    assert simplex_tree.num_vertices() == 4

    assert point_list[0] == filtered_alpha.get_point(0)
    assert point_list[1] == filtered_alpha.get_point(1)
    assert point_list[2] == filtered_alpha.get_point(2)
    assert point_list[3] == filtered_alpha.get_point(3)
    try:
        filtered_alpha.get_point(4) == []
    except IndexError:
        pass
    else:
        assert False
    try:
        filtered_alpha.get_point(125) == []
    except IndexError:
        pass
    else:
        assert False

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

def test_filtered_alpha():
    for precision in ['fast', 'safe', 'exact']:
        _filtered_alpha(precision)

def _safe_alpha_persistence_comparison(precision):
    #generate periodic signal
    time = np.arange(0, 10, 1)
    signal = [math.sin(x) for x in time]
    delta = math.pi
    delayed = [math.sin(x + delta) for x in time]
    
    #construct embedding
    embedding1 = [[signal[i], -signal[i]] for i in range(len(time))]
    embedding2 = [[signal[i], delayed[i]] for i in range(len(time))]
    
    #build alpha complex and simplex tree
    alpha_complex1 = gudhi.AlphaComplex(points=embedding1, precision = precision)
    simplex_tree1 = alpha_complex1.create_simplex_tree()
    
    alpha_complex2 = gudhi.AlphaComplex(points=embedding2, precision = precision)
    simplex_tree2 = alpha_complex2.create_simplex_tree()
    
    diag1 = simplex_tree1.persistence()
    diag2 = simplex_tree2.persistence()

    for (first_p, second_p) in zip_longest(diag1, diag2):
        assert first_p[0] == pytest.approx(second_p[0])
        assert first_p[1] == pytest.approx(second_p[1])


def test_safe_alpha_persistence_comparison():
    # Won't work for 'fast' version
    _safe_alpha_persistence_comparison('safe')
    _safe_alpha_persistence_comparison('exact')

def _delaunay_complex(precision):
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    filtered_alpha = gudhi.AlphaComplex(points=point_list, precision = precision)

    simplex_tree = filtered_alpha.create_simplex_tree(default_filtration_value = True)

    assert simplex_tree.num_simplices() == 11
    assert simplex_tree.num_vertices() == 4

    assert point_list[0] == filtered_alpha.get_point(0)
    assert point_list[1] == filtered_alpha.get_point(1)
    assert point_list[2] == filtered_alpha.get_point(2)
    assert point_list[3] == filtered_alpha.get_point(3)
    try:
        filtered_alpha.get_point(4) == []
    except IndexError:
        pass
    else:
        assert False
    try:
        filtered_alpha.get_point(125) == []
    except IndexError:
        pass
    else:
        assert False

    for filtered_value in simplex_tree.get_filtration():
        assert math.isnan(filtered_value[1])
    for filtered_value in simplex_tree.get_star([0]):
        assert math.isnan(filtered_value[1])
    for filtered_value in simplex_tree.get_cofaces([0], 1):
        assert math.isnan(filtered_value[1])

def test_delaunay_complex():
    for precision in ['fast', 'safe', 'exact']:
        _delaunay_complex(precision)

def _3d_points_on_a_plane(precision, default_filtration_value):
    alpha = gudhi.AlphaComplex(off_file=gudhi.__root_source_dir__ + '/data/points/alphacomplexdoc.off',
                         precision = precision)

    simplex_tree = alpha.create_simplex_tree(default_filtration_value = default_filtration_value)
    assert simplex_tree.dimension() == 2
    assert simplex_tree.num_vertices() == 7
    assert simplex_tree.num_simplices() == 25

def test_3d_points_on_a_plane():
    for default_filtration_value in [True, False]:
        for precision in ['fast', 'safe', 'exact']:
            _3d_points_on_a_plane(precision, default_filtration_value)