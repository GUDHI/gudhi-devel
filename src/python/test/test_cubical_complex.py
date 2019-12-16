from gudhi import CubicalComplex
import numpy as np

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


def test_empty_constructor():
    # Try to create an empty CubicalComplex
    cub = CubicalComplex()
    assert cub.__is_defined() == False
    assert cub.__is_persistence_defined() == False


def test_non_existing_perseus_file_constructor():
    # Try to open a non existing file
    cub = CubicalComplex(perseus_file="pouetpouettralala.toubiloubabdou")
    assert cub.__is_defined() == False
    assert cub.__is_persistence_defined() == False


def test_dimension_or_perseus_file_constructor():
    # Create test file
    test_file = open("CubicalOneSphere.txt", "w")
    test_file.write("2\n3\n3\n0\n0\n0\n0\n100\n0\n0\n0\n0\n")
    test_file.close()
    # CubicalComplex can be constructed from dimensions and
    # top_dimensional_cells OR from a Perseus-style file name.
    cub = CubicalComplex(
        dimensions=[3, 3],
        top_dimensional_cells=[1, 2, 3, 4, 5, 6, 7, 8, 9],
        perseus_file="CubicalOneSphere.txt",
    )
    assert cub.__is_defined() == False
    assert cub.__is_persistence_defined() == False

    cub = CubicalComplex(
        top_dimensional_cells=[1, 2, 3, 4, 5, 6, 7, 8, 9],
        perseus_file="CubicalOneSphere.txt",
    )
    assert cub.__is_defined() == False
    assert cub.__is_persistence_defined() == False

    cub = CubicalComplex(dimensions=[3, 3], perseus_file="CubicalOneSphere.txt")
    assert cub.__is_defined() == False
    assert cub.__is_persistence_defined() == False


def simple_constructor(cub):
    cub = CubicalComplex(
        dimensions=[3, 3], top_dimensional_cells=[1, 2, 3, 4, 5, 6, 7, 8, 9]
    )
    assert cub.__is_defined() == True
    assert cub.__is_persistence_defined() == False
    assert cub.persistence() == [(0, (1.0, float("inf")))]
    assert cub.__is_persistence_defined() == True
    assert cub.betti_numbers() == [1, 0, 0]
    assert cub.persistent_betti_numbers(0, 1000) == [0, 0, 0]

def test_simple_constructor_from_top_cells():
    cub = CubicalComplex(
        dimensions=[3, 3],
        top_dimensional_cells=[1, 2, 3, 4, 5, 6, 7, 8, 9],
    )
    simple_constructor(cub)

def test_simple_constructor_from_numpy_array():
    cub = CubicalComplex(
        numpy_array=np.array([[1, 2, 3],
                              [4, 5, 6],
                              [7, 8, 9]])
    )
    simple_constructor(cub)

def user_case_simple_constructor(cub):
    assert cub.__is_defined() == True
    assert cub.__is_persistence_defined() == False
    assert cub.persistence() == [(1, (0.0, 1.0)), (0, (0.0, float("inf")))]
    assert cub.__is_persistence_defined() == True
    other_cub = CubicalComplex(
        dimensions=[3, 3],
        top_dimensional_cells=[1000.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    )
    assert other_cub.persistence() == [(1, (0.0, 1.0)), (0, (0.0, float("inf")))]

def test_user_case_simple_constructor_from_top_cells():
    cub = CubicalComplex(
        dimensions=[3, 3],
        top_dimensional_cells=[float("inf"), 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    )
    user_case_simple_constructor(cub)

def test_user_case_simple_constructor_from_numpy_array():
    cub = CubicalComplex(
        numpy_array=np.array([[float("inf"), 0.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 0.0, 0.0]])
    )
    user_case_simple_constructor(cub)

def test_dimension_file_constructor():
    # Create test file
    test_file = open("CubicalOneSphere.txt", "w")
    test_file.write("2\n3\n3\n0\n0\n0\n0\n100\n0\n0\n0\n0\n")
    test_file.close()
    cub = CubicalComplex(perseus_file="CubicalOneSphere.txt")
    assert cub.__is_defined() == True
    assert cub.__is_persistence_defined() == False
    assert cub.persistence() == [(1, (0.0, 100.0)), (0, (0.0, float("inf")))]
    assert cub.__is_persistence_defined() == True
    assert cub.betti_numbers() == [1, 0, 0]
    assert cub.persistent_betti_numbers(0, 1000) == [1, 0, 0]
