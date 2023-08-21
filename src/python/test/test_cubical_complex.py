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
    with pytest.raises(ValueError):
        cub = CubicalComplex()

def test_non_existing_perseus_file_constructor():
    # Try to open a non existing file
    with pytest.raises(OSError):
        cub = CubicalComplex(perseus_file="pouetpouettralala.toubiloubabdou")

def test_dimensions_one_cells_type_or_perseus_file_constructor():
    # Create test file
    test_file = open("CubicalOneSphere.txt", "w")
    test_file.write("2\n3\n3\n0\n0\n0\n0\n100\n0\n0\n0\n0\n")
    test_file.close()
    # CubicalComplex can be constructed from dimensions and
    # cells OR from a Perseus-style file name.
    with pytest.raises(ValueError):
        cub = CubicalComplex(
            dimensions=[3, 3],
            top_dimensional_cells=[1, 2, 3, 4, 5, 6, 7, 8, 9],
            perseus_file="CubicalOneSphere.txt",
        )
    with pytest.raises(ValueError):
        cub = CubicalComplex(
            top_dimensional_cells=[1, 2, 3, 4, 5, 6, 7, 8, 9],
            perseus_file="CubicalOneSphere.txt",
        )
    with pytest.raises(ValueError):
        cub = CubicalComplex(dimensions=[3, 3],
                             perseus_file="CubicalOneSphere.txt")
    with pytest.raises(ValueError):
        cub = CubicalComplex(
            dimensions=[3, 3],
            top_dimensional_cells=[1, 2, 3, 4],
            vertices=[1, 2, 3, 4, 5, 6, 7, 8, 9],
            perseus_file="CubicalOneSphere.txt",
        )
    with pytest.raises(ValueError):
        cub = CubicalComplex(
            top_dimensional_cells=[1, 2, 3, 4],
            vertices=[1, 2, 3, 4, 5, 6, 7, 8, 9],
            perseus_file="CubicalOneSphere.txt",
        )
    with pytest.raises(ValueError):
        cub = CubicalComplex(
            dimensions=[3, 3],
            top_dimensional_cells=[1, 2, 3, 4],
            vertices=[1, 2, 3, 4, 5, 6, 7, 8, 9],
        )
    with pytest.raises(ValueError):
        cub = CubicalComplex(
            top_dimensional_cells=[1, 2, 3, 4],
            vertices=[1, 2, 3, 4, 5, 6, 7, 8, 9],
        )
    with pytest.raises(ValueError):
        cub = CubicalComplex(
            dimensions=[3, 3],
            vertices=[1, 2, 3, 4, 5, 6, 7, 8, 9],
            perseus_file="CubicalOneSphere.txt",
        )
    with pytest.raises(ValueError):
        cub = CubicalComplex(
            vertices=[1, 2, 3, 4, 5, 6, 7, 8, 9],
            perseus_file="CubicalOneSphere.txt",
        )

def simple_constructor(cub):
    assert cub._is_defined() == True
    assert cub._is_persistence_defined() == False
    assert cub.persistence() == [(0, (1.0, float("inf")))]
    assert cub._is_persistence_defined() == True
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

def test_constructor_from_vertices():
    cub = CubicalComplex(
        dimensions=[4, 4],
        vertices=[1, 1, 2, 3, 1, 1, 2, 3, 4, 4, 5, 6, 7, 7, 8, 9])
    simple_constructor(cub)
    a = np.random.rand(10, 20, 30)
    cub = PeriodicCubicalComplex(vertices=a, periodic_dimensions=(True, False, True))
    assert np.array_equal(a, cub.vertices())

def user_case_simple_constructor(cub):
    assert cub._is_defined() == True
    assert cub._is_persistence_defined() == False
    assert cub.persistence() == [(1, (0.0, 1.0)), (0, (0.0, float("inf")))]
    assert cub._is_persistence_defined() == True
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
    assert cub._is_defined() == True
    assert cub._is_persistence_defined() == False
    assert cub.persistence() == [(1, (0.0, 100.0)), (0, (0.0, float("inf")))]
    assert cub._is_persistence_defined() == True
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

def test_cubical_generators_from_top_dimensional_cells():
    cub = CubicalComplex(top_dimensional_cells = [[0, 0, 0], [0, 1, 0], [0, 0, 0]])
    cub.persistence()
    g = cub.cofaces_of_persistence_pairs()
    assert len(g[0]) == 2
    assert len(g[1]) == 1
    assert np.array_equal(g[0][0], np.empty(shape=[0,2]))
    assert np.array_equal(g[0][1], np.array([[7, 4]]))
    assert np.array_equal(g[1][0], np.array([8]))

def test_cubical_generators_from_vertices():
    cub = CubicalComplex(vertices = [[0, 0, 0, 0], [0, 1, 1, 0], [0, 1, 1, 0], [0, 0, 0, 0]])
    cub.persistence()
    g = cub.vertices_of_persistence_pairs()
    assert len(g[0]) == 2
    assert len(g[1]) == 1
    assert np.array_equal(g[0][0], np.empty(shape=[0,2]))
    assert np.array_equal(g[0][1], np.array([[14, 10]]))
    assert np.array_equal(g[1][0], np.array([14]))

    cub = CubicalComplex(vertices=[[2, 1, 6, 8],
                                   [5, 9, 4, 7],
                                   [4, 7, 3, 9]])
    assert set(cub.persistence()) == {(1, (7.0, 9.0)), (0, (1.0, np.inf)), (0, (3.0, 6.0)), (0, (4.0, 5.0))}
    g = cub.vertices_of_persistence_pairs()
    assert g[1] == [[3]]
    # The order [[8, 6], [2, 1]] would also be fine
    assert len(g[0]) == 2 and np.array_equal(g[0][0], [[2, 1], [8, 6]]) and np.array_equal(g[0][1], [[5, 4]])

def test_periodic_cubical_generators_from_vertices():
    cub = PeriodicCubicalComplex(vertices = [[0, 0, 0, 0], [0, 1, 1, 0], [0, 1, 1, 0], [0, 0, 0, 0]], periodic_dimensions=[False, False])
    cub.persistence()
    g = cub.vertices_of_persistence_pairs()
    assert len(g[0]) == 2
    assert len(g[1]) == 1
    assert np.array_equal(g[0][0], np.empty(shape=[0,2]))
    assert np.array_equal(g[0][1], np.array([[15, 10]]))
    assert np.array_equal(g[1][0], np.array([0]))

    v = np.array([[2, 1, 6, 8],
                  [5, 8, 4, 0],
                  [4, 7, 3, 9]], dtype=float)
    cub = PeriodicCubicalComplex(vertices=v, periodic_dimensions=(True, True))
    assert np.array_equal(
        cub.all_cells(),
        np.array(
            [
                [2.0, 2.0, 1.0, 6.0, 6.0, 8.0, 8.0, 8.0],
                [5.0, 8.0, 8.0, 8.0, 6.0, 8.0, 8.0, 8.0],
                [5.0, 8.0, 8.0, 8.0, 4.0, 4.0, 0.0, 5.0],
                [5.0, 8.0, 8.0, 8.0, 4.0, 9.0, 9.0, 9.0],
                [4.0, 7.0, 7.0, 7.0, 3.0, 9.0, 9.0, 9.0],
                [4.0, 7.0, 7.0, 7.0, 6.0, 9.0, 9.0, 9.0],
            ]
        )
    )
    p = cub.persistence()
    g = cub.vertices_of_persistence_pairs()
    assert set(p) == {
        (2, (9.0, np.inf)),
        (1, (5.0, np.inf)),
        (1, (6.0, np.inf)),
        (1, (6.0, 8.0)),
        (0, (0.0, np.inf)),
        (0, (1.0, 5.0)),
        (0, (3.0, 4.0)),
    }
    vr = np.ravel(v, order='F')
    diag = set()
    for (i,a) in enumerate(g[0]):
        for (b,d) in a:
            diag.add((i, (vr[b], vr[d])))
    for (i,a) in enumerate(g[1]):
        for b in a:
            diag.add((i, (vr[b], np.inf)))
    assert set(p) == diag

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

def test_cubical_persistence_intervals_in_dimension():
    cub = CubicalComplex(
        dimensions=[3, 3],
        top_dimensional_cells=[1, 2, 3, 4, 5, 6, 7, 8, 9],
    )
    cub.compute_persistence()
    H0 = cub.persistence_intervals_in_dimension(0)
    assert np.array_equal(H0, np.array([[ 1., float("inf")]]))
    assert cub.persistence_intervals_in_dimension(1).shape == (0, 2)

def test_periodic_cubical_persistence_intervals_in_dimension():
    cub = PeriodicCubicalComplex(
        dimensions=[3, 3],
        top_dimensional_cells=[1, 2, 3, 4, 5, 6, 7, 8, 9],
        periodic_dimensions = [True, True]
    )
    cub.compute_persistence()
    H0 = cub.persistence_intervals_in_dimension(0)
    assert np.array_equal(H0, np.array([[ 1., float("inf")]]))
    H1 = cub.persistence_intervals_in_dimension(1)
    assert np.array_equal(H1, np.array([[ 3., float("inf")], [ 7., float("inf")]]))
    H2 = cub.persistence_intervals_in_dimension(2)
    assert np.array_equal(H2, np.array([[ 9., float("inf")]]))
    assert cub.persistence_intervals_in_dimension(3).shape == (0, 2)

def test_array_access():
    a = np.arange(6, dtype=float).reshape(3, 2)
    cplx = CubicalComplex(top_dimensional_cells=a)
    assert np.array_equal(cplx.top_dimensional_cells(), a)
    assert np.array_equal(
        cplx.vertices(), np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [2.0, 2.0, 3.0], [4.0, 4.0, 5.0]])
    )
    assert np.array_equal(
        cplx.all_cells(),
        np.array(
            [
                [0.0, 0.0, 0.0, 1.0, 1.0],
                [0.0, 0.0, 0.0, 1.0, 1.0],
                [0.0, 0.0, 0.0, 1.0, 1.0],
                [2.0, 2.0, 2.0, 3.0, 3.0],
                [2.0, 2.0, 2.0, 3.0, 3.0],
                [4.0, 4.0, 4.0, 5.0, 5.0],
                [4.0, 4.0, 4.0, 5.0, 5.0],
            ]
        ),
    )
    cplx.top_dimensional_cells()[0, 0] = 42
    assert cplx.all_cells()[1, 1] == 42

    a = np.array([[1.0, 2.0], [4.0, 3.0], [0.0, 5.0]])
    cplx = PeriodicCubicalComplex(top_dimensional_cells=a, periodic_dimensions=(True, False))
    assert np.array_equal(cplx.top_dimensional_cells(), a)
    assert np.array_equal(cplx.vertices(), np.array([[0.0, 0.0, 2.0], [1.0, 1.0, 2.0], [0.0, 0.0, 3.0]]))
    assert np.array_equal(
        cplx.all_cells(),
        np.array(
            [
                [0.0, 0.0, 0.0, 2.0, 2.0],
                [1.0, 1.0, 1.0, 2.0, 2.0],
                [1.0, 1.0, 1.0, 2.0, 2.0],
                [4.0, 4.0, 3.0, 3.0, 3.0],
                [0.0, 0.0, 0.0, 3.0, 3.0],
                [0.0, 0.0, 0.0, 5.0, 5.0]
            ]
        ),
    )
