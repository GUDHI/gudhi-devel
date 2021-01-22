""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi import CubicalComplex, PeriodicCubicalComplex, CubicalComplexBaseClass
import numpy as np
import pytest

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"


def test_empty_constructor():
    with pytest.raises(TypeError):
        cub = CubicalComplexBaseClass()
    with pytest.raises(TypeError):
        cub = CubicalComplex()
    with pytest.raises(TypeError):
        cub = PeriodicCubicalComplex()


def test_non_existing_perseus_file_constructor():
    # Try to open a non existing file
    with pytest.raises(FileNotFoundError):
        cub = CubicalComplex(perseus_file="pouetpouettralala.toubiloubabdou")
    with pytest.raises(FileNotFoundError):
        cub = PeriodicCubicalComplex(perseus_file="pouetpouettralala.toubiloubabdou")


def test_dimension_or_perseus_file_constructor():
    # Create test file
    test_file = open("CubicalOneSphere.txt", "w")
    test_file.write("2\n3\n3\n0\n0\n0\n0\n100\n0\n0\n0\n0\n")
    test_file.close()
    # CubicalComplex can be constructed from dimensions and
    # top_dimensional_cells OR from a Perseus-style file name.
    with pytest.raises(TypeError):
        cub = CubicalComplex(
            dimensions=[3, 3],
            top_dimensional_cells=[1, 2, 3, 4, 5, 6, 7, 8, 9],
            perseus_file="CubicalOneSphere.txt",
        )

    with pytest.raises(TypeError):
        cub = CubicalComplex(
            top_dimensional_cells=[1, 2, 3, 4, 5, 6, 7, 8, 9],
            perseus_file="CubicalOneSphere.txt",
        )

    with pytest.raises(TypeError):
        cub = CubicalComplex(dimensions=[3, 3], perseus_file="CubicalOneSphere.txt")


def simple_constructor(cub):
    assert cub.dimension() == 2
    assert cub.num_simplices() == 49
    assert cub.persistence() == [(0, (1.0, float("inf")))]
    assert cub.betti_numbers() == [1, 0, 0]
    assert cub.persistent_betti_numbers(0, 1000) == [0, 0, 0]

def test_simple_constructor_from_top_cells():
    cub = CubicalComplex(
        dimensions=[3, 3],
        top_dimensional_cells=[1, 2, 3, 4, 5, 6, 7, 8, 9])
    simple_constructor(cub)

def test_simple_constructor_from_numpy_array():
    cub = CubicalComplex(
        top_dimensional_cells=np.array([[1, 2, 3],
                                        [4, 5, 6],
                                        [7, 8, 9]])
    )
    simple_constructor(cub)

def user_case_simple_constructor(cub):
    assert cub.dimension() == 2
    assert cub.num_simplices() == 49
    assert cub.persistence() == [(1, (0.0, 1.0)), (0, (0.0, float("inf")))]

def test_user_case_simple_constructor_from_top_cells():
    cub = CubicalComplex(
        dimensions=[3, 3],
        top_dimensional_cells=[float("inf"), 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0])
    user_case_simple_constructor(cub)
    other_cub = CubicalComplex(
        dimensions=[3, 3],
        top_dimensional_cells=[1000.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0])
    user_case_simple_constructor(other_cub)

def test_user_case_simple_constructor_from_numpy_array():
    cub = CubicalComplex(
        top_dimensional_cells=np.array([[float("inf"), 0.0, 0.0],
                                        [0.0, 1.0, 0.0],
                                        [0.0, 0.0, 0.0]])
    )
    user_case_simple_constructor(cub)
    other_cub = CubicalComplex(
        top_dimensional_cells=np.array([[1000.0, 0.0, 0.0],
                                        [0.0, 1.0, 0.0],
                                        [0.0, 0.0, 0.0]])
    )
    user_case_simple_constructor(other_cub)

def test_dimension_file_constructor():
    # Create test file
    test_file = open("CubicalOneSphere.txt", "w")
    test_file.write("2\n3\n3\n0\n0\n0\n0\n100\n0\n0\n0\n0\n")
    test_file.close()
    cub = CubicalComplex(perseus_file="CubicalOneSphere.txt")
    assert cub.persistence() == [(1, (0.0, 100.0)), (0, (0.0, float("inf")))]
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
    g = cub.cofaces_of_persistence_pairs()
    assert len(g[0]) == 2
    assert len(g[1]) == 1
    assert np.array_equal(g[0][0], np.empty(shape=[0,2]))
    assert np.array_equal(g[0][1], np.array([[7, 4]]))
    assert np.array_equal(g[1][0], np.array([8]))

def test_cubical_cofaces_of_persistence_pairs_when_pd_has_no_paired_birth_and_death():
    cubCpx = CubicalComplex(dimensions=[1,2], top_dimensional_cells=[0.0, 1.0])
    Diag = cubCpx.persistence(homology_coeff_field=2, min_persistence=0)
    pairs = cubCpx.cofaces_of_persistence_pairs()
    assert pairs[0] == []
    assert np.array_equal(pairs[1][0], np.array([0]))

def test_periodic_cofaces_of_persistence_pairs_when_pd_has_no_paired_birth_and_death():
    perCubCpx = PeriodicCubicalComplex(dimensions=[1,2], top_dimensional_cells=[0.0, 1.0],
                                       periodic_dimensions=[True, True])
    Diag = perCubCpx.persistence(homology_coeff_field=2, min_persistence=0)
    pairs = perCubCpx.cofaces_of_persistence_pairs()
    assert pairs[0] == []
    assert np.array_equal(pairs[1][0], np.array([0]))
    assert np.array_equal(pairs[1][1], np.array([0, 1]))
    assert np.array_equal(pairs[1][2], np.array([1]))

def test_periodic_dimension_file_constructor():
    # Create test file
    test_file = open("2-torus.txt", "w")
    test_file.write("2\n-3\n-3\n0\n0\n0\n0\n10\n0\n0\n0\n0\n")
    test_file.close()
    cub = PeriodicCubicalComplex(perseus_file="2-torus.txt")
    assert cub.dimension() == 2
    assert cub.num_simplices() == 36
    assert cub.persistence() == [(2, (10.0, float("inf"))), (1, (0.0, float("inf"))),
                                 (1, (0.0, float("inf"))), (0, (0.0, float("inf")))]
    assert cub.betti_numbers() == [1, 2, 1]
    assert cub.persistent_betti_numbers(0, 1000) == [1, 2, 0]
    assert np.array_equal(cub.persistence_intervals_in_dimension(1), np.array([[ 0., float("inf")],
                                                                               [ 0., float("inf")]]))
