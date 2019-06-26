from gudhi import AlphaComplex, SimplexTree

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"


def test_empty_alpha():
    alpha_complex = AlphaComplex(points=[[0,0]])
    assert alpha_complex.__is_defined() == True

def test_infinite_alpha():
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    alpha_complex = AlphaComplex(points=point_list)
    assert alpha_complex.__is_defined() == True

    simplex_tree = alpha_complex.create_simplex_tree()
    assert simplex_tree.__is_persistence_defined() == False

    assert simplex_tree.num_simplices() == 11
    assert simplex_tree.num_vertices() == 4
 
    assert simplex_tree.get_filtration() == \
           [([0], 0.0), ([1], 0.0), ([2], 0.0), ([3], 0.0),
            ([0, 1], 0.25), ([0, 2], 0.25), ([1, 3], 0.25),
            ([2, 3], 0.25), ([1, 2], 0.5), ([0, 1, 2], 0.5),
            ([1, 2, 3], 0.5)]
    assert simplex_tree.get_star([0]) == \
           [([0], 0.0), ([0, 1], 0.25), ([0, 1, 2], 0.5),
           ([0, 2], 0.25)]
    assert simplex_tree.get_cofaces([0], 1) == \
           [([0, 1], 0.25), ([0, 2], 0.25)]
 
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

    assert simplex_tree.get_filtration() == \
           [([0], 0.0), ([1], 0.0), ([2], 0.0), ([3], 0.0),
            ([0, 1], 0.25), ([0, 2], 0.25), ([1, 3], 0.25),
            ([2, 3], 0.25)]
    assert simplex_tree.get_star([0]) == \
           [([0], 0.0), ([0, 1], 0.25), ([0, 2], 0.25)]
    assert simplex_tree.get_cofaces([0], 1) == \
           [([0, 1], 0.25), ([0, 2], 0.25)]
