from gudhi import CubicalComplex

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


def test_empty_constructor():
    # Try to create an empty CubicalComplex
    cub = CubicalComplex()
    assert cub.__is_defined() == False
    assert cub.__is_persistence_defined() == False

def test_non_existing_perseus_file_constructor():
    # Try to open a non existing file
    cub = CubicalComplex(perseus_file='pouetpouettralala.toubiloubabdou')
    assert cub.__is_defined() == False
    assert cub.__is_persistence_defined() == False

def test_dimension_or_perseus_file_constructor():
    # Create test file
    test_file = open('CubicalOneSphere.txt', 'w')
    test_file.write('2\n3\n3\n0\n0\n0\n0\n100\n0\n0\n0\n0\n')
    test_file.close()
    # CubicalComplex can be constructed from dimensions and
    # top_dimensional_cells OR from a perseus file style name.
    cub = CubicalComplex(dimensions=[3, 3],
                         top_dimensional_cells = [1,2,3,4,5,6,7,8,9],
                         perseus_file='CubicalOneSphere.txt')
    assert cub.__is_defined() == False
    assert cub.__is_persistence_defined() == False

    cub = CubicalComplex(top_dimensional_cells = [1,2,3,4,5,6,7,8,9],
                         perseus_file='CubicalOneSphere.txt')
    assert cub.__is_defined() == False
    assert cub.__is_persistence_defined() == False

    cub = CubicalComplex(dimensions=[3, 3],
                         perseus_file='CubicalOneSphere.txt')
    assert cub.__is_defined() == False
    assert cub.__is_persistence_defined() == False

def test_dimension_constructor():
    cub = CubicalComplex(dimensions=[3, 3],
                         top_dimensional_cells = [1,2,3,4,5,6,7,8,9])
    assert cub.__is_defined() == True
    assert cub.__is_persistence_defined() == False
    assert cub.persistence() == [(1, (0.0, 100.0)), (0, (0.0, 1.8446744073709552e+19))]
    assert cub.__is_persistence_defined() == True
    assert cub.betti_numbers() == [1, 0]
    assert cub.persistent_betti_numbers(0, 1000) == [0, 0]

def test_dimension_constructor():
    # Create test file
    test_file = open('CubicalOneSphere.txt', 'w')
    test_file.write('2\n3\n3\n0\n0\n0\n0\n100\n0\n0\n0\n0\n')
    test_file.close()
    cub = CubicalComplex(perseus_file='CubicalOneSphere.txt')
    assert cub.__is_defined() == True
    assert cub.__is_persistence_defined() == False
    assert cub.persistence() == [(1, (0.0, 100.0)), (0, (0.0, 1.8446744073709552e+19))]
    assert cub.__is_persistence_defined() == True
    assert cub.betti_numbers() == [1, 0]
    assert cub.persistent_betti_numbers(0, 1000) == [1, 0]
