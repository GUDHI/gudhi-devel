from gudhi import SimplexTree

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


def test_insertion():
    st = SimplexTree()
    assert st.__is_defined() == True
    assert st.__is_persistence_defined() == False

    # insert test
    assert st.insert([0, 1]) == True

    assert st.dimension() == 1

    assert st.insert([0, 1, 2], filtration=4.0) == True

    assert st.dimension() == 2

    assert st.num_simplices() == 7
    assert st.num_vertices() == 3

    # find test
    assert st.find([0, 1, 2]) == True
    assert st.find([0, 1]) == True
    assert st.find([0, 2]) == True
    assert st.find([0]) == True
    assert st.find([1]) == True
    assert st.find([2]) == True
    assert st.find([3])    == False
    assert st.find([0, 3]) == False
    assert st.find([1, 3]) == False
    assert st.find([2, 3]) == False

    # filtration test
    st.set_filtration(5.0)
    st.initialize_filtration()
    assert st.filtration([0, 1, 2]) == 4.0
    assert st.filtration([0, 2]) == 4.0
    assert st.filtration([1, 2]) == 4.0
    assert st.filtration([2]) == 4.0
    assert st.filtration([0, 1]) == 0.0
    assert st.filtration([0]) == 0.0
    assert st.filtration([1]) == 0.0

    # skeleton test
    assert st.get_skeleton(2) == \
        [([0, 1, 2], 4.0), ([0, 1], 0.0), ([0, 2], 4.0),
        ([0], 0.0), ([1, 2], 4.0), ([1], 0.0), ([2], 4.0)]
    assert st.get_skeleton(1) == \
        [([0, 1], 0.0), ([0, 2], 4.0), ([0], 0.0),
        ([1, 2], 4.0), ([1], 0.0), ([2], 4.0)]
    assert st.get_skeleton(0) == \
        [([0], 0.0), ([1], 0.0), ([2], 4.0)]

    # remove_maximal_simplex test
    assert st.get_cofaces([0, 1, 2], 1) == []
    st.remove_maximal_simplex([0, 1, 2])
    assert st.get_skeleton(2) == \
        [([0, 1], 0.0), ([0, 2], 4.0), ([0], 0.0),
        ([1, 2], 4.0), ([1], 0.0), ([2], 4.0)]
    assert st.find([0, 1, 2]) == False
    assert st.find([0, 1]) == True
    assert st.find([0, 2]) == True
    assert st.find([0]) == True
    assert st.find([1]) == True
    assert st.find([2]) == True

    st.initialize_filtration()
    assert st.persistence(persistence_dim_max = True) == [(1, (4.0, float('inf'))), (0, (0.0, float('inf')))]
    assert st.__is_persistence_defined() == True

    assert st.betti_numbers() == [1, 1]
    assert st.persistent_betti_numbers(-0.1, 10000.0) == [0, 0]
    assert st.persistent_betti_numbers(0.0, 10000.0) == [1, 0]
    assert st.persistent_betti_numbers(3.9, 10000.0) == [1, 0]
    assert st.persistent_betti_numbers(4.0, 10000.0) == [1, 1]
    assert st.persistent_betti_numbers(9999.0, 10000.0) == [1, 1]

def test_expansion():
    st = SimplexTree()
    assert st.__is_defined() == True
    assert st.__is_persistence_defined() == False

    # insert test
    assert st.insert([3, 2], 0.1) == True
    assert st.insert([2, 0], 0.2) == True
    assert st.insert([1, 0], 0.3) == True
    assert st.insert([3, 1], 0.4) == True
    assert st.insert([2, 1], 0.5) == True
    assert st.insert([6, 5], 0.6) == True
    assert st.insert([4, 2], 0.7) == True
    assert st.insert([3, 0], 0.8) == True
    assert st.insert([6, 4], 0.9) == True
    assert st.insert([6, 3], 1.0) == True

    assert st.num_vertices() == 7
    assert st.num_simplices() == 17
    assert st.get_filtration() == [([2], 0.1), ([3], 0.1), ([2, 3], 0.1),
    ([0], 0.2), ([0, 2], 0.2), ([1], 0.3), ([0, 1], 0.3), ([1, 3], 0.4),
    ([1, 2], 0.5), ([5], 0.6), ([6], 0.6), ([5, 6], 0.6), ([4], 0.7),
    ([2, 4], 0.7), ([0, 3], 0.8), ([4, 6], 0.9), ([3, 6], 1.0)]

    st.expansion(3)
    assert st.num_vertices() == 7
    assert st.num_simplices() == 22
    st.initialize_filtration()

    assert st.get_filtration() == [([2], 0.1), ([3], 0.1), ([2, 3], 0.1),
    ([0], 0.2), ([0, 2], 0.2), ([1], 0.3), ([0, 1], 0.3), ([1, 3], 0.4),
    ([1, 2], 0.5), ([0, 1, 2], 0.5), ([1, 2, 3], 0.5), ([5], 0.6), ([6], 0.6),
    ([5, 6], 0.6), ([4], 0.7), ([2, 4], 0.7), ([0, 3], 0.8), ([0, 1, 3], 0.8),
    ([0, 2, 3], 0.8), ([0, 1, 2, 3], 0.8), ([4, 6], 0.9), ([3, 6], 1.0)]
