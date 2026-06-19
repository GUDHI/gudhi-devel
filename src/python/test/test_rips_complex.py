""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - 2025/04 Hannah Schreiber: Add tests to verify possibility of tensor input
      - 2026/06 Joo-Heung Nahm and Vincent Rouvreau: Add tests for F-contiguous and non contiguous numpy arrays
      - YYYY/MM Author: Description of the modification
"""

from math import sqrt
import numpy as np

from gudhi import RipsComplex


def test_empty_rips():
    rips_complex = RipsComplex()


def _test_rips_from_points(point_list):
    rips_complex = RipsComplex(points=point_list, max_edge_length=42)

    simplex_tree = rips_complex.create_simplex_tree(max_dimension=1)

    assert simplex_tree._is_persistence_defined() == False

    assert simplex_tree.num_simplices() == 10
    assert simplex_tree.num_vertices() == 4

    assert list(simplex_tree.get_filtration()) == [
        ([0], 0.0),
        ([1], 0.0),
        ([2], 0.0),
        ([3], 0.0),
        ([0, 1], 1.0),
        ([0, 2], 1.0),
        ([1, 3], 1.0),
        ([2, 3], 1.0),
        ([1, 2], 1.4142135623730951),
        ([0, 3], 1.4142135623730951),
    ]

    assert simplex_tree.get_star([0]) == [
        ([0], 0.0),
        ([0, 1], 1.0),
        ([0, 2], 1.0),
        ([0, 3], 1.4142135623730951),
    ]
    assert simplex_tree.get_cofaces([0], 1) == [
        ([0, 1], 1.0),
        ([0, 2], 1.0),
        ([0, 3], 1.4142135623730951),
    ]


def test_rips_from_points():
    _test_rips_from_points([[0, 0], [1, 0], [0, 1], [1, 1]])


def test_rips_from_numpy_points():
    _test_rips_from_points(np.array([[0, 0], [1, 0], [0, 1], [1, 1]]))


def _test_filtered_rips_from_points(point_list):
    filtered_rips = RipsComplex(points=point_list, max_edge_length=1.0)

    simplex_tree = filtered_rips.create_simplex_tree(max_dimension=1)
    assert simplex_tree._is_persistence_defined() == False

    assert simplex_tree.num_simplices() == 8
    assert simplex_tree.num_vertices() == 4


def test_filtered_rips_from_points():
    _test_filtered_rips_from_points([[0, 0], [1, 0], [0, 1], [1, 1]])


def test_filtered_rips_from_numpy_points():
    _test_filtered_rips_from_points(np.array([[0, 0], [1, 0], [0, 1], [1, 1]]))


def _test_sparse_filtered_rips_from_points(point_list):
    filtered_rips = RipsComplex(points=point_list, max_edge_length=1.0, sparse=0.001)

    simplex_tree = filtered_rips.create_simplex_tree(max_dimension=1)

    assert simplex_tree._is_persistence_defined() == False

    assert simplex_tree.num_simplices() == 8
    assert simplex_tree.num_vertices() == 4


def test_sparse_filtered_rips_from_points():
    _test_sparse_filtered_rips_from_points([[0, 0], [1, 0], [0, 1], [1, 1]])


def test_sparse_filtered_rips_from_numpy_points():
    _test_sparse_filtered_rips_from_points(np.array([[0, 0], [1, 0], [0, 1], [1, 1]]))


def _test_rips_from_distance_matrix(distance_matrix):
    rips_complex = RipsComplex(distance_matrix=distance_matrix, max_edge_length=42)

    simplex_tree = rips_complex.create_simplex_tree(max_dimension=1)

    assert simplex_tree._is_persistence_defined() == False

    assert simplex_tree.num_simplices() == 10
    assert simplex_tree.num_vertices() == 4

    assert list(simplex_tree.get_filtration()) == [
        ([0], 0.0),
        ([1], 0.0),
        ([2], 0.0),
        ([3], 0.0),
        ([0, 1], 1.0),
        ([0, 2], 1.0),
        ([1, 3], 1.0),
        ([2, 3], 1.0),
        ([1, 2], 1.4142135623730951),
        ([0, 3], 1.4142135623730951),
    ]

    assert simplex_tree.get_star([0]) == [
        ([0], 0.0),
        ([0, 1], 1.0),
        ([0, 2], 1.0),
        ([0, 3], 1.4142135623730951),
    ]
    assert simplex_tree.get_cofaces([0], 1) == [
        ([0, 1], 1.0),
        ([0, 2], 1.0),
        ([0, 3], 1.4142135623730951),
    ]


def test_rips_from_distance_matrix():
    _test_rips_from_distance_matrix([[0], [1, 0], [1, sqrt(2), 0], [sqrt(2), 1, 1, 0]])


def test_rips_from_numpy_distance_matrix():
    _test_rips_from_distance_matrix(
        np.array([[0, 0, 0, 0], [1, 0, 0, 0], [1, sqrt(2), 0, 0], [sqrt(2), 1, 1, 0]])
    )


def _test_filtered_rips_from_distance_matrix(distance_matrix):
    filtered_rips = RipsComplex(distance_matrix=distance_matrix, max_edge_length=1.0)

    simplex_tree = filtered_rips.create_simplex_tree(max_dimension=1)

    assert simplex_tree._is_persistence_defined() == False

    assert simplex_tree.num_simplices() == 8
    assert simplex_tree.num_vertices() == 4


def test_filtered_rips_from_distance_matrix():
    _test_filtered_rips_from_distance_matrix(
        [[0], [1, 0], [1, sqrt(2), 0], [sqrt(2), 1, 1, 0]]
    )


def test_filtered_rips_from_numpy_distance_matrix():
    _test_filtered_rips_from_distance_matrix(
        np.array([[0, 0, 0, 0], [1, 0, 0, 0], [1, sqrt(2), 0, 0], [sqrt(2), 1, 1, 0]])
    )


def _test_sparse_with_multiplicity(points):
    rips = RipsComplex(points=points, sparse=0.01)
    simplex_tree = rips.create_simplex_tree(max_dimension=2)
    assert simplex_tree.num_simplices() == 7
    diag = simplex_tree.persistence()


def test_sparse_with_multiplicity():
    _test_sparse_with_multiplicity(
        [
            [3, 4],
            [0.1, 2],
            [0.1, 2],
            [0.1, 2],
            [0.1, 2],
            [0.1, 2],
            [0.1, 2],
            [0.1, 2],
            [0.1, 2],
            [0.1, 2],
            [0.1, 2],
            [3, 4.1],
        ]
    )


def test_sparse_with_numpy_multiplicity():
    _test_sparse_with_multiplicity(
        np.array(
            [
                [3, 4],
                [0.1, 2],
                [0.1, 2],
                [0.1, 2],
                [0.1, 2],
                [0.1, 2],
                [0.1, 2],
                [0.1, 2],
                [0.1, 2],
                [0.1, 2],
                [0.1, 2],
                [3, 4.1],
            ]
        )
    )


def test_sparse_with_numpy_transposed_multiplicity():
    points = np.array(
        [
            [3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 3],
            [4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4.1],
        ]
    )
    _test_sparse_with_multiplicity(points.transpose())


def test_tensors():
    try:
        import torch

        points = (torch.rand((5, 2)) * 2 - 1).requires_grad_()
        rips = RipsComplex(points=points)
    except ImportError:
        pass

    try:
        import tensorflow as tf

        points = tf.random.uniform(shape=[5, 2])
        rips = RipsComplex(points=points)
    except ImportError:
        pass


def test_contiguity_for_points():
    rng = np.random.default_rng(42)
    base = rng.standard_normal((10, 2))
    
    # Need to test with a small epsilon value for complete rips - or require a seed to be set
    for sparse in [None, 0.01]:
        st_c = RipsComplex(points=np.ascontiguousarray(base),
                           sparse=sparse).create_simplex_tree(max_dimension=2) # C-contiguous
        st_f = RipsComplex(points=np.asfortranarray(base),
                           sparse=sparse).create_simplex_tree(max_dimension=2) # F-contiguous (same data)
        
        assert st_c == st_f
        
        base_repeat = np.repeat(base, 2, axis=0)
        st_r = RipsComplex(points=base_repeat[::2],
                           sparse=sparse).create_simplex_tree(max_dimension=2) # no contiguity (but same data)
        
        assert st_c == st_r


def test_contiguity_for_distance_matrix():
    rng = np.random.default_rng(42)
    # Get a lower triangle 10x10 random positivevalues
    base = np.absolute(np.tril(rng.standard_normal((10, 10)), k=-1))
    
    # Need to test with a small epsilon value for complete rips - or require a seed to be set
    for sparse in [None, 0.01]:
        st_c = RipsComplex(distance_matrix=np.ascontiguousarray(base),
                           sparse=sparse).create_simplex_tree(max_dimension=2) # C-contiguous
        st_f = RipsComplex(distance_matrix=np.asfortranarray(base),
                           sparse=sparse).create_simplex_tree(max_dimension=2) # F-contiguous (same data)
        
        assert st_c == st_f
        
        base_repeat = np.repeat(base, 2, axis=0)
        st_r = RipsComplex(distance_matrix=base_repeat[::2],
                           sparse=sparse).create_simplex_tree(max_dimension=2) # no contiguity (but same data)
        
        assert st_c == st_r
