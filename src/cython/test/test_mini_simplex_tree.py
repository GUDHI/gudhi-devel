from gudhi import MiniSimplexTree

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


def test_mini():
    triangle012 = [0, 1, 2]
    edge03 = [0, 3]
    mini_st = MiniSimplexTree()
    assert mini_st.__is_defined() == True
    assert mini_st.__is_persistence_defined() == False
    assert mini_st.insert(triangle012) == True
    assert mini_st.insert(edge03) == True
    # FIXME: Remove this line
    mini_st.set_dimension(2)

    edge02 = [0, 2]
    assert mini_st.find(edge02) == True
    assert mini_st.get_coface_tree(edge02, 1) == \
            [([0, 1, 2], 0.0)]

    # remove_maximal_simplex test
    assert mini_st.get_coface_tree(triangle012, 1) == []
    mini_st.remove_maximal_simplex(triangle012)
    assert mini_st.find(edge02) == True
    assert mini_st.find(triangle012) == False
