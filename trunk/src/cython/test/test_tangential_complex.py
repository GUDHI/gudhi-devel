from gudhi import TangentialComplex, SimplexTree

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


def test_tangential():
    point_list = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
    tc = TangentialComplex(points=point_list)
    assert tc.__is_defined() == True
    assert tc.num_vertices() == 4

    st = tc.create_simplex_tree()
    assert st.__is_defined() == True
    assert st.__is_persistence_defined() == False

    assert st.num_simplices() == 6
    assert st.num_vertices() == 4
 
    assert st.get_filtration() == \
        [([0], 0.0), ([1], 0.0), ([2], 0.0), ([0, 2], 0.0), ([3], 0.0), ([1, 3], 0.0)]
    assert st.get_cofaces([0], 1) == [([0, 2], 0.0)]
 
    assert point_list[0] == tc.get_point(0)
    assert point_list[1] == tc.get_point(1)
    assert point_list[2] == tc.get_point(2)
    assert point_list[3] == tc.get_point(3)
    assert tc.get_point(4) == []
    assert tc.get_point(125) == []
