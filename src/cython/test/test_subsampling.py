import gudhi
import os

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


def test_write_off_file_for_tests():
    file = open("n_farthest.off", "w")
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
    point_set = [[0,1], [0,0], [1,0], [1,1]]
    i = 0
    for point in point_set:
        # The iteration starts with the given starting point
        sub_set = gudhi.choose_n_farthest_points(points = point_set, nb_points = 1, starting_point = i)
        assert sub_set[0] == point_set[i]
        i = i + 1

    # The iteration finds then the farthest
    sub_set = gudhi.choose_n_farthest_points(points = point_set, nb_points = 2, starting_point = 1)
    assert sub_set[1] == point_set[3]
    sub_set = gudhi.choose_n_farthest_points(points = point_set, nb_points = 2, starting_point = 3)
    assert sub_set[1] == point_set[1]
    sub_set = gudhi.choose_n_farthest_points(points = point_set, nb_points = 2, starting_point = 0)
    assert sub_set[1] == point_set[2]
    sub_set = gudhi.choose_n_farthest_points(points = point_set, nb_points = 2, starting_point = 2)
    assert sub_set[1] == point_set[0]

    # Test the limits
    assert gudhi.choose_n_farthest_points(points = [], nb_points = 0, starting_point = 0) == []
    assert gudhi.choose_n_farthest_points(points = [], nb_points = 1, starting_point = 0) == []
    assert gudhi.choose_n_farthest_points(points = [], nb_points = 0, starting_point = 1) == []
    assert gudhi.choose_n_farthest_points(points = [], nb_points = 1, starting_point = 1) == []

    print(os.getcwd())
    # From off file test
    for i in range (0, 7):
        assert len(gudhi.choose_n_farthest_points(off_file = 'n_farthest.off', nb_points = i, starting_point = i)) == i

def test_simple_choose_n_farthest_points_randomed():
    point_set = [[0,1], [0,0], [1,0], [1,1]]

    # Test the limits
    assert gudhi.choose_n_farthest_points(points = [], nb_points = 0) == []
    assert gudhi.choose_n_farthest_points(points = [], nb_points = 1) == []
    assert gudhi.choose_n_farthest_points(points = point_set, nb_points = 0) == []
    # Go furter than point set on purpose
    for iter in range(1,10):
        sub_set = gudhi.choose_n_farthest_points(points = point_set, nb_points = iter)
        for sub in sub_set:
            found = False
            for point in point_set:
                if point == sub:
                    found = True
            assert found == True

    print(os.getcwd())
    # From off file test
    for i in range (0, 7):
        assert len(gudhi.choose_n_farthest_points(off_file = 'n_farthest.off', nb_points = i)) == i
