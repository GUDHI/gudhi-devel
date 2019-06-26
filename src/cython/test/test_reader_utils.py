import gudhi
import numpy as np

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2017 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2017 Inria"
__license__ = "MIT"


def test_non_existing_csv_file():
    # Try to open a non existing file
    matrix = gudhi.read_lower_triangular_matrix_from_csv_file(
        csv_file="pouetpouettralala.toubiloubabdou"
    )
    assert matrix == []


def test_full_square_distance_matrix_csv_file():
    # Create test file
    test_file = open("full_square_distance_matrix.csv", "w")
    test_file.write("0;1;2;3;\n1;0;4;5;\n2;4;0;6;\n3;5;6;0;")
    test_file.close()
    matrix = gudhi.read_lower_triangular_matrix_from_csv_file(
        csv_file="full_square_distance_matrix.csv"
    )
    assert matrix == [[], [1.0], [2.0, 4.0], [3.0, 5.0, 6.0]]


def test_lower_triangular_distance_matrix_csv_file():
    # Create test file
    test_file = open("lower_triangular_distance_matrix.csv", "w")
    test_file.write("\n1,\n2,3,\n4,5,6,\n7,8,9,10,")
    test_file.close()
    matrix = gudhi.read_lower_triangular_matrix_from_csv_file(
        csv_file="lower_triangular_distance_matrix.csv", separator=","
    )
    assert matrix == [[], [1.0], [2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0, 10.0]]


def test_non_existing_persistence_file():
    # Try to open a non existing file
    persistence = gudhi.read_persistence_intervals_grouped_by_dimension(
        persistence_file="pouetpouettralala.toubiloubabdou"
    )
    assert persistence == []
    persistence = gudhi.read_persistence_intervals_in_dimension(
        persistence_file="pouetpouettralala.toubiloubabdou", only_this_dim=1
    )
    np.testing.assert_array_equal(persistence, [])


def test_read_persistence_intervals_without_dimension():
    # Create test file
    test_file = open("persistence_intervals_without_dimension.pers", "w")
    test_file.write(
        "# Simple persistence diagram without dimension\n2.7 3.7\n9.6 14.\n34.2 34.974\n3. inf"
    )
    test_file.close()
    persistence = gudhi.read_persistence_intervals_in_dimension(
        persistence_file="persistence_intervals_without_dimension.pers"
    )
    np.testing.assert_array_equal(
        persistence, [(2.7, 3.7), (9.6, 14.0), (34.2, 34.974), (3.0, float("Inf"))]
    )
    persistence = gudhi.read_persistence_intervals_in_dimension(
        persistence_file="persistence_intervals_without_dimension.pers", only_this_dim=0
    )
    np.testing.assert_array_equal(persistence, [])
    persistence = gudhi.read_persistence_intervals_in_dimension(
        persistence_file="persistence_intervals_without_dimension.pers", only_this_dim=1
    )
    np.testing.assert_array_equal(persistence, [])
    persistence = gudhi.read_persistence_intervals_grouped_by_dimension(
        persistence_file="persistence_intervals_without_dimension.pers"
    )
    assert persistence == {
        -1: [(2.7, 3.7), (9.6, 14.0), (34.2, 34.974), (3.0, float("Inf"))]
    }


def test_read_persistence_intervals_with_dimension():
    # Create test file
    test_file = open("persistence_intervals_with_dimension.pers", "w")
    test_file.write(
        "# Simple persistence diagram with dimension\n0 2.7 3.7\n1 9.6 14.\n3 34.2 34.974\n1 3. inf"
    )
    test_file.close()
    persistence = gudhi.read_persistence_intervals_in_dimension(
        persistence_file="persistence_intervals_with_dimension.pers"
    )
    np.testing.assert_array_equal(
        persistence, [(2.7, 3.7), (9.6, 14.0), (34.2, 34.974), (3.0, float("Inf"))]
    )
    persistence = gudhi.read_persistence_intervals_in_dimension(
        persistence_file="persistence_intervals_with_dimension.pers", only_this_dim=0
    )
    np.testing.assert_array_equal(persistence, [(2.7, 3.7)])
    persistence = gudhi.read_persistence_intervals_in_dimension(
        persistence_file="persistence_intervals_with_dimension.pers", only_this_dim=1
    )
    np.testing.assert_array_equal(persistence, [(9.6, 14.0), (3.0, float("Inf"))])
    persistence = gudhi.read_persistence_intervals_in_dimension(
        persistence_file="persistence_intervals_with_dimension.pers", only_this_dim=2
    )
    np.testing.assert_array_equal(persistence, [])
    persistence = gudhi.read_persistence_intervals_in_dimension(
        persistence_file="persistence_intervals_with_dimension.pers", only_this_dim=3
    )
    np.testing.assert_array_equal(persistence, [(34.2, 34.974)])
    persistence = gudhi.read_persistence_intervals_grouped_by_dimension(
        persistence_file="persistence_intervals_with_dimension.pers"
    )
    assert persistence == {
        0: [(2.7, 3.7)],
        1: [(9.6, 14.0), (3.0, float("Inf"))],
        3: [(34.2, 34.974)],
    }
