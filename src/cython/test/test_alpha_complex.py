from gudhi import AlphaComplex

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


def test_empty_alpha():
    alpha_complex = AlphaComplex()
    assert alpha_complex.__is_defined() == True
    assert alpha_complex.__is_persistence_defined() == False
    assert alpha_complex.num_simplices() == 0
    assert alpha_complex.num_vertices() == 0

def test_infinite_alpha():
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    alpha_complex = AlphaComplex(points=point_list)
    assert alpha_complex.__is_defined() == True
    assert alpha_complex.__is_persistence_defined() == False

    assert alpha_complex.num_simplices() == 11
    assert alpha_complex.num_vertices() == 4
 
    assert alpha_complex.get_filtered_tree() == \
           [([0], 0.0), ([1], 0.0), ([2], 0.0), ([3], 0.0),
            ([0, 1], 0.25), ([0, 2], 0.25), ([1, 3], 0.25),
            ([2, 3], 0.25), ([1, 2], 0.5), ([0, 1, 2], 0.5),
            ([1, 2, 3], 0.5)]
    assert alpha_complex.get_star_tree([0]) == \
           [([0], 0.0), ([0, 1], 0.25), ([0, 1, 2], 0.5),
           ([0, 2], 0.25)]
    assert alpha_complex.get_coface_tree([0], 1) == \
           [([0, 1], 0.25), ([0, 2], 0.25)]
 
    assert point_list[0] == alpha_complex.get_point(0)
    assert point_list[1] == alpha_complex.get_point(1)
    assert point_list[2] == alpha_complex.get_point(2)
    assert point_list[3] == alpha_complex.get_point(3)
    assert alpha_complex.get_point(4) == []
    assert alpha_complex.get_point(125) == []

def test_filtered_alpha():
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    filtered_alpha = AlphaComplex(points=point_list,
                                  max_alpha_square=0.25)

    assert filtered_alpha.num_simplices() == 8
    assert filtered_alpha.num_vertices() == 4

    assert point_list[0] == filtered_alpha.get_point(0)
    assert point_list[1] == filtered_alpha.get_point(1)
    assert point_list[2] == filtered_alpha.get_point(2)
    assert point_list[3] == filtered_alpha.get_point(3)
    assert filtered_alpha.get_point(4) == []
    assert filtered_alpha.get_point(125) == []

    assert filtered_alpha.get_filtered_tree() == \
           [([0], 0.0), ([1], 0.0), ([2], 0.0), ([3], 0.0),
            ([0, 1], 0.25), ([0, 2], 0.25), ([1, 3], 0.25),
            ([2, 3], 0.25)]
    assert filtered_alpha.get_star_tree([0]) == \
           [([0], 0.0), ([0, 1], 0.25), ([0, 2], 0.25)]
    assert filtered_alpha.get_coface_tree([0], 1) == \
           [([0, 1], 0.25), ([0, 2], 0.25)]
