from gudhi import WitnessComplex, StrongWitnessComplex, SimplexTree

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


def test_empty_witness_complex():
    witness = WitnessComplex()
    assert witness.__is_defined() == False

def test_witness_complex():
    nearest_landmark_table = [[[0, 0], [1, 1], [2, 2], [3, 3], [4, 4]],
                              [[1, 0], [2, 1], [3, 2], [4, 3], [0, 4]],
                              [[2, 0], [3, 1], [4, 2], [0, 3], [1, 4]],
                              [[3, 0], [4, 1], [0, 2], [1, 3], [2, 4]],
                              [[4, 0], [0, 1], [1, 2], [2, 3], [3, 4]]]

    witness_complex = WitnessComplex(nearest_landmark_table=nearest_landmark_table)
    simplex_tree = witness_complex.create_simplex_tree(max_alpha_square=4.1)
    assert simplex_tree.num_vertices() == 5
    assert simplex_tree.num_simplices() == 31
    simplex_tree = witness_complex.create_simplex_tree(max_alpha_square=4.1, limit_dimension=2)
    assert simplex_tree.num_vertices() == 5
    assert simplex_tree.num_simplices() == 25

def test_strong_witness_complex():
    nearest_landmark_table = [[[0, 0], [1, 1], [2, 2], [3, 3], [4, 4]],
                              [[1, 0], [2, 1], [3, 2], [4, 3], [0, 4]],
                              [[2, 0], [3, 1], [4, 2], [0, 3], [1, 4]],
                              [[3, 0], [4, 1], [0, 2], [1, 3], [2, 4]],
                              [[4, 0], [0, 1], [1, 2], [2, 3], [3, 4]]]

    strong_witness_complex = StrongWitnessComplex(nearest_landmark_table=nearest_landmark_table)
    simplex_tree = strong_witness_complex.create_simplex_tree(max_alpha_square=4.1)
    assert simplex_tree.num_vertices() == 5
    assert simplex_tree.num_simplices() == 31
    simplex_tree = strong_witness_complex.create_simplex_tree(max_alpha_square=4.1, limit_dimension=2)
    assert simplex_tree.num_vertices() == 5
    assert simplex_tree.num_simplices() == 25
