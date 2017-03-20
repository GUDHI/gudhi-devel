from gudhi import RipsComplex
from math import sqrt

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Vincent Rouvreau

   Copyright (C) 2016 INRIA

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 INRIA"
__license__ = "GPL v3"


def test_empty_rips():
    rips_complex = RipsComplex()
    assert rips_complex.__is_defined() == True

def test_rips_from_points():
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    rips_complex = RipsComplex(points=point_list, max_edge_length=42)

    simplex_tree = rips_complex.create_simplex_tree(max_dimension=1)

    assert simplex_tree.__is_defined() == True
    assert simplex_tree.__is_persistence_defined() == False

    assert simplex_tree.num_simplices() == 10
    assert simplex_tree.num_vertices() == 4

    assert simplex_tree.get_filtered_tree() == \
           [([0], 0.0), ([1], 0.0), ([2], 0.0), ([3], 0.0),
            ([0, 1], 1.0), ([0, 2], 1.0), ([1, 3], 1.0),
            ([2, 3], 1.0), ([1, 2], 1.4142135623730951),
            ([0, 3], 1.4142135623730951)]
    assert simplex_tree.get_star_tree([0]) == \
           [([0], 0.0), ([0, 1], 1.0), ([0, 2], 1.0),
            ([0, 3], 1.4142135623730951)]
    assert simplex_tree.get_cofaces([0], 1) == \
           [([0, 1], 1.0), ([0, 2], 1.0),
            ([0, 3], 1.4142135623730951)]

def test_filtered_rips_from_points():
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    filtered_rips = RipsComplex(points=point_list, max_edge_length=1.0)

    simplex_tree = filtered_rips.create_simplex_tree(max_dimension=1)

    assert simplex_tree.__is_defined() == True
    assert simplex_tree.__is_persistence_defined() == False

    assert simplex_tree.num_simplices() == 8
    assert simplex_tree.num_vertices() == 4

def test_rips_from_distance_matrix():
    distance_matrix = [[0],
                       [1, 0],
                       [1, sqrt(2), 0],
                       [sqrt(2), 1, 1, 0]]
    rips_complex = RipsComplex(distance_matrix=distance_matrix, max_edge_length=42)

    simplex_tree = rips_complex.create_simplex_tree(max_dimension=1)

    assert simplex_tree.__is_defined() == True
    assert simplex_tree.__is_persistence_defined() == False

    assert simplex_tree.num_simplices() == 10
    assert simplex_tree.num_vertices() == 4

    assert simplex_tree.get_filtered_tree() == \
           [([0], 0.0), ([1], 0.0), ([2], 0.0), ([3], 0.0),
            ([0, 1], 1.0), ([0, 2], 1.0), ([1, 3], 1.0),
            ([2, 3], 1.0), ([1, 2], 1.4142135623730951),
            ([0, 3], 1.4142135623730951)]
    assert simplex_tree.get_star_tree([0]) == \
           [([0], 0.0), ([0, 1], 1.0), ([0, 2], 1.0),
            ([0, 3], 1.4142135623730951)]
    assert simplex_tree.get_cofaces([0], 1) == \
           [([0, 1], 1.0), ([0, 2], 1.0),
            ([0, 3], 1.4142135623730951)]

def test_filtered_rips_from_distance_matrix():
    distance_matrix = [[0],
                       [1, 0],
                       [1, sqrt(2), 0],
                       [sqrt(2), 1, 1, 0]]
    filtered_rips = RipsComplex(distance_matrix=distance_matrix, max_edge_length=1.0)

    simplex_tree = filtered_rips.create_simplex_tree(max_dimension=1)

    assert simplex_tree.__is_defined() == True
    assert simplex_tree.__is_persistence_defined() == False

    assert simplex_tree.num_simplices() == 8
    assert simplex_tree.num_vertices() == 4
