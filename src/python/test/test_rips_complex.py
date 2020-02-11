""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi import RipsComplex
from math import sqrt
import pytest

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"


def test_empty_rips():
    rips_complex = RipsComplex()


def test_rips_from_points():
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    rips_complex = RipsComplex(points=point_list, max_edge_length=42)

    simplex_tree = rips_complex.create_simplex_tree(max_dimension=1)

    assert simplex_tree.__is_defined() == True
    assert simplex_tree.__is_persistence_defined() == False

    assert simplex_tree.num_simplices() == 10
    assert simplex_tree.num_vertices() == 4

    assert list(simplex_tree.get_filtration()) == [
        ([0], 0.0),
        ([1], 0.0),
        ([2], 0.0),
        ([3], 0.0),
        ([0, 1], 1.0),
        ([0, 2], 1.0),
        ([1, 3], 1.0),
        ([2, 3], 1.0),
        ([1, 2], 1.4142135623730951),
        ([0, 3], 1.4142135623730951),
    ]

    assert simplex_tree.get_star([0]) == [
        ([0], 0.0),
        ([0, 1], 1.0),
        ([0, 2], 1.0),
        ([0, 3], 1.4142135623730951),
    ]
    assert simplex_tree.get_cofaces([0], 1) == [
        ([0, 1], 1.0),
        ([0, 2], 1.0),
        ([0, 3], 1.4142135623730951),
    ]


def test_filtered_rips_from_points():
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    filtered_rips = RipsComplex(points=point_list, max_edge_length=1.0)

    simplex_tree = filtered_rips.create_simplex_tree(max_dimension=1)

    assert simplex_tree.__is_defined() == True
    assert simplex_tree.__is_persistence_defined() == False

    assert simplex_tree.num_simplices() == 8
    assert simplex_tree.num_vertices() == 4


def test_sparse_filtered_rips_from_points():
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    filtered_rips = RipsComplex(points=point_list, max_edge_length=1.0, sparse=0.001)

    simplex_tree = filtered_rips.create_simplex_tree(max_dimension=1)

    assert simplex_tree.__is_defined() == True
    assert simplex_tree.__is_persistence_defined() == False

    assert simplex_tree.num_simplices() == 8
    assert simplex_tree.num_vertices() == 4


def test_rips_from_distance_matrix():
    distance_matrix = [[0], [1, 0], [1, sqrt(2), 0], [sqrt(2), 1, 1, 0]]
    rips_complex = RipsComplex(distance_matrix=distance_matrix, max_edge_length=42)

    simplex_tree = rips_complex.create_simplex_tree(max_dimension=1)

    assert simplex_tree.__is_defined() == True
    assert simplex_tree.__is_persistence_defined() == False

    assert simplex_tree.num_simplices() == 10
    assert simplex_tree.num_vertices() == 4

    assert list(simplex_tree.get_filtration()) == [
        ([0], 0.0),
        ([1], 0.0),
        ([2], 0.0),
        ([3], 0.0),
        ([0, 1], 1.0),
        ([0, 2], 1.0),
        ([1, 3], 1.0),
        ([2, 3], 1.0),
        ([1, 2], 1.4142135623730951),
        ([0, 3], 1.4142135623730951),
    ]

    assert simplex_tree.get_star([0]) == [
        ([0], 0.0),
        ([0, 1], 1.0),
        ([0, 2], 1.0),
        ([0, 3], 1.4142135623730951),
    ]
    assert simplex_tree.get_cofaces([0], 1) == [
        ([0, 1], 1.0),
        ([0, 2], 1.0),
        ([0, 3], 1.4142135623730951),
    ]


def test_filtered_rips_from_distance_matrix():
    distance_matrix = [[0], [1, 0], [1, sqrt(2), 0], [sqrt(2), 1, 1, 0]]
    filtered_rips = RipsComplex(distance_matrix=distance_matrix, max_edge_length=1.0)

    simplex_tree = filtered_rips.create_simplex_tree(max_dimension=1)

    assert simplex_tree.__is_defined() == True
    assert simplex_tree.__is_persistence_defined() == False

    assert simplex_tree.num_simplices() == 8
    assert simplex_tree.num_vertices() == 4
