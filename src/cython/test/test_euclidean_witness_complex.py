import gudhi

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Vincent Rouvreau

   Copyright (C) 2016 Inria

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
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "GPL v3"


def test_empty_euclidean_witness_complex():
    euclidean_witness = gudhi.EuclideanWitnessComplex()
    assert euclidean_witness.__is_defined() == False

def test_witness_complex():
    point_cloud = [[1.0, 1.0], [7.0, 0.0], [4.0, 6.0], [9.0, 6.0],
        [0.0, 14.0], [2.0, 19.0], [9.0, 17.0]]
    landmarks = [[1.0, 1.0], [7.0, 0.0], [4.0, 6.0]]
    euclidean_witness_complex = gudhi.EuclideanWitnessComplex(landmarks=landmarks, witnesses = point_cloud)
    simplex_tree = euclidean_witness_complex.create_simplex_tree(max_alpha_square=4.1)

    assert landmarks[0] == euclidean_witness_complex.get_point(0)
    assert landmarks[1] == euclidean_witness_complex.get_point(1)
    assert landmarks[2] == euclidean_witness_complex.get_point(2)

    assert simplex_tree.get_filtration() == [([0], 0.0), ([1], 0.0),
        ([0, 1], 0.0), ([2], 0.0), ([0, 2], 0.0), ([1, 2], 0.0),
        ([0, 1, 2], 0.0)]

def test_empty_euclidean_strong_witness_complex():
    euclidean_strong_witness = gudhi.EuclideanStrongWitnessComplex()
    assert euclidean_strong_witness.__is_defined() == False

def test_strong_witness_complex():
    point_cloud = [[1.0, 1.0], [7.0, 0.0], [4.0, 6.0], [9.0, 6.0],
        [0.0, 14.0], [2.0, 19.0], [9.0, 17.0]]
    landmarks = [[1.0, 1.0], [7.0, 0.0], [4.0, 6.0]]
    euclidean_strong_witness_complex = gudhi.EuclideanStrongWitnessComplex(landmarks=landmarks, witnesses = point_cloud)
    simplex_tree = euclidean_strong_witness_complex.create_simplex_tree(max_alpha_square=14.9)

    assert landmarks[0] == euclidean_strong_witness_complex.get_point(0)
    assert landmarks[1] == euclidean_strong_witness_complex.get_point(1)
    assert landmarks[2] == euclidean_strong_witness_complex.get_point(2)

    assert simplex_tree.get_filtration() == [([0], 0.0), ([1], 0.0), ([2], 0.0)]

    simplex_tree = euclidean_strong_witness_complex.create_simplex_tree(max_alpha_square=100.0)

    assert simplex_tree.get_filtration() == [([0], 0.0), ([1], 0.0),
    ([2], 0.0), ([1, 2], 15.0), ([0, 2], 34.0), ([0, 1], 37.0),
    ([0, 1, 2], 37.0)]

