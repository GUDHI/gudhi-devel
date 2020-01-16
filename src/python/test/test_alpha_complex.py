""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi import AlphaComplex, SimplexTree
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


def test_empty_alpha():
    alpha_complex = AlphaComplex(points=[[0, 0]])
    assert alpha_complex.__is_defined() == True


def test_infinite_alpha():
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    alpha_complex = AlphaComplex(points=point_list)
    assert alpha_complex.__is_defined() == True

    simplex_tree = alpha_complex.create_simplex_tree()
    assert simplex_tree.__is_persistence_defined() == False

    assert simplex_tree.num_simplices() == 11
    assert simplex_tree.num_vertices() == 4

    assert simplex_tree.get_filtration() == [
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
    assert alpha_complex.get_point(4) == []
    assert alpha_complex.get_point(125) == []


def test_filtered_alpha():
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    filtered_alpha = AlphaComplex(points=point_list)

    simplex_tree = filtered_alpha.create_simplex_tree(max_alpha_square=0.25)

    assert simplex_tree.num_simplices() == 8
    assert simplex_tree.num_vertices() == 4

    assert point_list[0] == filtered_alpha.get_point(0)
    assert point_list[1] == filtered_alpha.get_point(1)
    assert point_list[2] == filtered_alpha.get_point(2)
    assert point_list[3] == filtered_alpha.get_point(3)
    assert filtered_alpha.get_point(4) == []
    assert filtered_alpha.get_point(125) == []

    assert simplex_tree.get_filtration() == [
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

def test_safe_alpha_persistence_comparison():
    #generate periodic signal
    time = np.arange(0, 10, 1)
    signal = [math.sin(x) for x in time]
    delta = math.pi
    delayed = [math.sin(x + delta) for x in time]
    
    #construct embedding
    embedding1 = [[signal[i], -signal[i]] for i in range(len(time))]
    embedding2 = [[signal[i], delayed[i]] for i in range(len(time))]
    
    #build alpha complex and simplex tree
    alpha_complex1 = AlphaComplex(points=embedding1)
    simplex_tree1 = alpha_complex1.create_simplex_tree()
    
    alpha_complex2 = AlphaComplex(points=embedding2)
    simplex_tree2 = alpha_complex2.create_simplex_tree()
    
    diag1 = simplex_tree1.persistence()
    diag2 = simplex_tree2.persistence()

    for (first_p, second_p) in zip_longest(diag1, diag2):
        assert first_p[0] == pytest.approx(second_p[0])
        assert first_p[1] == pytest.approx(second_p[1])
