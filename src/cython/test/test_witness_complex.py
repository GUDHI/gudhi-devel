from gudhi import WitnessComplex

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
    point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
    witness = WitnessComplex(points=point_list, number_of_landmarks=10)
    assert witness.__is_defined() == True

    # FIXME: Remove this line
    witness.set_dimension(2)

    assert witness.num_simplices() == 13
    assert witness.num_vertices() == 10
    witness.initialize_filtration()
"""
    assert witness.get_filtered_tree() == \
        [([0], 0.0), ([1], 0.0), ([0, 1], 0.0), ([2], 0.0), ([1, 2], 0.0),
         ([3], 0.0), ([2, 3], 0.0), ([4], 0.0), ([5], 0.0), ([6], 0.0),
         ([7], 0.0), ([8], 0.0), ([9], 0.0)]


    assert witness.get_coface_tree([2], 1) == [([1, 2], 0.0), ([2, 3], 0.0)]
    assert witness.get_star_tree([2]) == \
        [([1, 2], 0.0), ([2], 0.0), ([2, 3], 0.0)]
"""
