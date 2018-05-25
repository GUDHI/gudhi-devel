import gudhi

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Vincent Rouvreau

   Copyright (C) 2017 Inria

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
__copyright__ = "Copyright (C) 2017 Inria"
__license__ = "GPL v3"


def test_non_existing_csv_file():
    # Try to open a non existing file
    matrix = gudhi.read_lower_triangular_matrix_from_csv_file(csv_file='pouetpouettralala.toubiloubabdou')
    assert matrix == []

def test_full_square_distance_matrix_csv_file():
    # Create test file
    test_file = open('full_square_distance_matrix.csv', 'w')
    test_file.write('0;1;2;3;\n1;0;4;5;\n2;4;0;6;\n3;5;6;0;')
    test_file.close()
    matrix = gudhi.read_lower_triangular_matrix_from_csv_file(csv_file="full_square_distance_matrix.csv")
    assert matrix == [[], [1.0], [2.0, 4.0], [3.0, 5.0, 6.0]]

def test_lower_triangular_distance_matrix_csv_file():
    # Create test file
    test_file = open('lower_triangular_distance_matrix.csv', 'w')
    test_file.write('\n1,\n2,3,\n4,5,6,\n7,8,9,10,')
    test_file.close()
    matrix = gudhi.read_lower_triangular_matrix_from_csv_file(csv_file="lower_triangular_distance_matrix.csv", separator=",")
    assert matrix == [[], [1.0], [2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0, 10.0]]

def test_non_existing_persistence_file():
    # Try to open a non existing file
    persistence = gudhi.read_persistence_intervals_grouped_by_dimension(persistence_file='pouetpouettralala.toubiloubabdou')
    assert persistence == []
    persistence = gudhi.read_persistence_intervals_in_dimension(persistence_file='pouetpouettralala.toubiloubabdou', only_this_dim=1)
    assert persistence == []

def test_read_persistence_intervals_without_dimension():
    # Create test file
    test_file = open('persistence_intervals_without_dimension.pers', 'w')
    test_file.write('# Simple persistence diagram without dimension\n2.7 3.7\n9.6 14.\n34.2 34.974\n3. inf')
    test_file.close()
    persistence = gudhi.read_persistence_intervals_in_dimension(persistence_file='persistence_intervals_without_dimension.pers')
    assert persistence == [(2.7, 3.7), (9.6, 14.), (34.2, 34.974), (3., float('Inf'))]
    persistence = gudhi.read_persistence_intervals_in_dimension(persistence_file='persistence_intervals_without_dimension.pers', only_this_dim=0)
    assert persistence == []
    persistence = gudhi.read_persistence_intervals_in_dimension(persistence_file='persistence_intervals_without_dimension.pers', only_this_dim=1)
    assert persistence == []
    persistence = gudhi.read_persistence_intervals_grouped_by_dimension(persistence_file='persistence_intervals_without_dimension.pers')
    assert persistence == {-1: [(2.7, 3.7), (9.6, 14.0), (34.2, 34.974), (3.0, float('Inf'))]}

def test_read_persistence_intervals_with_dimension():
    # Create test file
    test_file = open('persistence_intervals_with_dimension.pers', 'w')
    test_file.write('# Simple persistence diagram with dimension\n0 2.7 3.7\n1 9.6 14.\n3 34.2 34.974\n1 3. inf')
    test_file.close()
    persistence = gudhi.read_persistence_intervals_in_dimension(persistence_file='persistence_intervals_with_dimension.pers')
    assert persistence == [(2.7, 3.7), (9.6, 14.), (34.2, 34.974), (3., float('Inf'))]
    persistence = gudhi.read_persistence_intervals_in_dimension(persistence_file='persistence_intervals_with_dimension.pers', only_this_dim=0)
    assert persistence == [(2.7, 3.7)]
    persistence = gudhi.read_persistence_intervals_in_dimension(persistence_file='persistence_intervals_with_dimension.pers', only_this_dim=1)
    assert persistence == [(9.6, 14.), (3., float('Inf'))]
    persistence = gudhi.read_persistence_intervals_in_dimension(persistence_file='persistence_intervals_with_dimension.pers', only_this_dim=2)
    assert persistence == []
    persistence = gudhi.read_persistence_intervals_in_dimension(persistence_file='persistence_intervals_with_dimension.pers', only_this_dim=3)
    assert persistence == [(34.2, 34.974)]
    persistence = gudhi.read_persistence_intervals_grouped_by_dimension(persistence_file='persistence_intervals_with_dimension.pers')
    assert persistence == {0: [(2.7, 3.7)], 1: [(9.6, 14.0), (3.0, float('Inf'))], 3: [(34.2, 34.974)]}
