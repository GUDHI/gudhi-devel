""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

import gudhi

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"


def test_write_off_file_for_tests():
    file = open("subsample.off", "w")
    file.write("nOFF\n")
    file.write("2 7 0 0\n")
    file.write("1.0 1.0\n")
    file.write("7.0 0.0\n")
    file.write("4.0 6.0\n")
    file.write("9.0 6.0\n")
    file.write("0.0 14.0\n")
    file.write("2.0 19.0\n")
    file.write("9.0 17.0\n")
    file.close()


def test_simple_choose_n_farthest_points_with_a_starting_point():
    point_set = [[0, 1], [0, 0], [1, 0], [1, 1]]
    i = 0
    for point in point_set:
        # The iteration starts with the given starting point
        sub_set = gudhi.choose_n_farthest_points(
            points=point_set, nb_points=1, starting_point=i
        )
        assert sub_set[0] == point_set[i]
        i = i + 1

    # The iteration finds then the farthest
    sub_set = gudhi.choose_n_farthest_points(
        points=point_set, nb_points=2, starting_point=1
    )
    assert sub_set[1] == point_set[3]
    sub_set = gudhi.choose_n_farthest_points(
        points=point_set, nb_points=2, starting_point=3
    )
    assert sub_set[1] == point_set[1]
    sub_set = gudhi.choose_n_farthest_points(
        points=point_set, nb_points=2, starting_point=0
    )
    assert sub_set[1] == point_set[2]
    sub_set = gudhi.choose_n_farthest_points(
        points=point_set, nb_points=2, starting_point=2
    )
    assert sub_set[1] == point_set[0]

    # Test the limits
    assert (
        gudhi.choose_n_farthest_points(points=[], nb_points=0, starting_point=0) == []
    )
    assert (
        gudhi.choose_n_farthest_points(points=[], nb_points=1, starting_point=0) == []
    )
    assert (
        gudhi.choose_n_farthest_points(points=[], nb_points=0, starting_point=1) == []
    )
    assert (
        gudhi.choose_n_farthest_points(points=[], nb_points=1, starting_point=1) == []
    )

    # From off file test
    for i in range(0, 7):
        assert (
            len(
                gudhi.choose_n_farthest_points(
                    off_file="subsample.off", nb_points=i, starting_point=i
                )
            )
            == i
        )


def test_simple_choose_n_farthest_points_randomed():
    point_set = [[0, 1], [0, 0], [1, 0], [1, 1]]
    # Test the limits
    assert gudhi.choose_n_farthest_points(points=[], nb_points=0) == []
    assert gudhi.choose_n_farthest_points(points=[], nb_points=1) == []
    assert gudhi.choose_n_farthest_points(points=point_set, nb_points=0) == []

    # Go furter than point set on purpose
    for iter in range(1, 10):
        sub_set = gudhi.choose_n_farthest_points(points=point_set, nb_points=iter)
        for sub in sub_set:
            found = False
            for point in point_set:
                if point == sub:
                    found = True
            # Check each sub set point is existing in the point set
            assert found == True

    # From off file test
    for i in range(0, 7):
        assert (
            len(gudhi.choose_n_farthest_points(off_file="subsample.off", nb_points=i))
            == i
        )


def test_simple_pick_n_random_points():
    point_set = [[0, 1], [0, 0], [1, 0], [1, 1]]
    # Test the limits
    assert gudhi.pick_n_random_points(points=[], nb_points=0) == []
    assert gudhi.pick_n_random_points(points=[], nb_points=1) == []
    assert gudhi.pick_n_random_points(points=point_set, nb_points=0) == []

    # Go furter than point set on purpose
    for iter in range(1, 10):
        sub_set = gudhi.pick_n_random_points(points=point_set, nb_points=iter)
        for sub in sub_set:
            found = False
            for point in point_set:
                if point == sub:
                    found = True
            # Check each sub set point is existing in the point set
            assert found == True

    # From off file test
    for i in range(0, 7):
        assert (
            len(gudhi.pick_n_random_points(off_file="subsample.off", nb_points=i)) == i
        )


def test_simple_sparsify_points():
    point_set = [[0, 1], [0, 0], [1, 0], [1, 1]]
    # Test the limits
    # assert gudhi.sparsify_point_set(points = [], min_squared_dist = 0.0) == []
    # assert gudhi.sparsify_point_set(points = [], min_squared_dist = 10.0) == []
    assert gudhi.sparsify_point_set(points=point_set, min_squared_dist=0.0) == point_set
    assert gudhi.sparsify_point_set(points=point_set, min_squared_dist=0.999) == point_set
    assert gudhi.sparsify_point_set(points=point_set, min_squared_dist=1.001) == [
        [0, 1],
        [1, 0],
    ]
    assert gudhi.sparsify_point_set(points=point_set, min_squared_dist=1.999) == [
        [0, 1],
        [1, 0],
    ]
    assert gudhi.sparsify_point_set(points=point_set, min_squared_dist=2.001) == [[0, 1]]

    assert (
        len(gudhi.sparsify_point_set(off_file="subsample.off", min_squared_dist=0.0))
        == 7
    )
    assert (
        len(gudhi.sparsify_point_set(off_file="subsample.off", min_squared_dist=30.0))
        == 5
    )
    assert (
        len(gudhi.sparsify_point_set(off_file="subsample.off", min_squared_dist=40.1))
        == 4
    )
    assert (
        len(gudhi.sparsify_point_set(off_file="subsample.off", min_squared_dist=89.9))
        == 3
    )
    assert (
        len(gudhi.sparsify_point_set(off_file="subsample.off", min_squared_dist=100.0))
        == 2
    )
    assert (
        len(gudhi.sparsify_point_set(off_file="subsample.off", min_squared_dist=324.9))
        == 2
    )
    assert (
        len(gudhi.sparsify_point_set(off_file="subsample.off", min_squared_dist=325.01))
        == 1
    )
