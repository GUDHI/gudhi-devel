""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi import CubicalComplex, PeriodicCubicalComplex
import numpy as np
import pytest

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
    with pytest.raises(FileNotFoundError):
        cub = CubicalComplex(perseus_file="pouetpouettralala.toubiloubabdou")


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
        top_dimensional_cells=np.array([[1, 2, 3],
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
        top_dimensional_cells=np.array([[float("inf"), 0.0, 0.0],
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

def test_connected_sublevel_sets():
    array_cells = np.array([[3, 3], [2, 2], [4, 4]])
    linear_cells = [3, 3, 2, 2, 4, 4]
    dimensions = [2, 3]
    periodic_dimensions = [False, False]
    # with a numpy array version
    cub = CubicalComplex(top_dimensional_cells = array_cells)
    assert cub.persistence() == [(0, (2.0, float("inf")))]
    assert cub.betti_numbers() == [1, 0, 0]
    # with vector of dimensions
    cub = CubicalComplex(dimensions = dimensions,
                         top_dimensional_cells = linear_cells)
    assert cub.persistence() == [(0, (2.0, float("inf")))]
    assert cub.betti_numbers() == [1, 0, 0]
    # periodic with a numpy array version
    cub = PeriodicCubicalComplex(top_dimensional_cells = array_cells,
                                periodic_dimensions = periodic_dimensions)
    assert cub.persistence() == [(0, (2.0, float("inf")))]
    assert cub.betti_numbers() == [1, 0, 0]
    # periodic with vector of dimensions
    cub = PeriodicCubicalComplex(dimensions = dimensions,
                                 top_dimensional_cells = linear_cells,
                                 periodic_dimensions = periodic_dimensions)
    assert cub.persistence() == [(0, (2.0, float("inf")))]
    assert cub.betti_numbers() == [1, 0, 0]

def test_cubical_generators():
    cub = CubicalComplex(top_dimensional_cells = [[0, 0, 0], [0, 1, 0], [0, 0, 0]])
    cub.persistence()
    assert np.array_equal(cub.cofaces_of_persistence_pairs(), np.array([[1, 7, 4], [0, 8, -1]]))
