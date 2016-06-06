import unittest

import gudhi

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


class TestSimplexTree(unittest.TestCase):

    def test_insertion(self):
        st = gudhi.SimplexTree()

        # insert test
        self.assertTrue(st.insert([0, 1]))
        self.assertTrue(st.insert([0, 1, 2], filtration=4.0))
        self.assertEqual(st.num_simplices(), 7)
        self.assertEqual(st.num_vertices(), 3)

        # find test
        self.assertTrue(st.find([0, 1, 2]))
        self.assertTrue(st.find([0, 1]))
        self.assertTrue(st.find([0, 2]))
        self.assertTrue(st.find([0]))
        self.assertTrue(st.find([1]))
        self.assertTrue(st.find([2]))
        self.assertFalse(st.find([3]))
        self.assertFalse(st.find([0, 3]))
        self.assertFalse(st.find([1, 3]))
        self.assertFalse(st.find([2, 3]))

        # filtration test
        st.set_filtration(5.0)
        st.initialize_filtration()
        self.assertEqual(st.get_filtration(), 5.0)
        self.assertEqual(st.filtration([0, 1, 2]), 4.0)
        self.assertEqual(st.filtration([0, 2]), 4.0)
        self.assertEqual(st.filtration([1, 2]), 4.0)
        self.assertEqual(st.filtration([2]), 4.0)
        self.assertEqual(st.filtration([0, 1]), 0.0)
        self.assertEqual(st.filtration([0]), 0.0)
        self.assertEqual(st.filtration([1]), 0.0)

        # skeleton_tree test
        self.assertEqual(st.get_skeleton_tree(2),
                         [([0, 1, 2], 4.0), ([0, 1], 0.0), ([0, 2], 4.0),
                          ([0], 0.0), ([1, 2], 4.0), ([1], 0.0), ([2], 4.0)])
        self.assertEqual(st.get_skeleton_tree(1),
                         [([0, 1], 0.0), ([0, 2], 4.0), ([0], 0.0),
                          ([1, 2], 4.0), ([1], 0.0), ([2], 4.0)])
        self.assertEqual(st.get_skeleton_tree(0),
                         [([0], 0.0), ([1], 0.0), ([2], 4.0)])

        # remove_maximal_simplex test
        self.assertEqual(st.get_coface_tree([0, 1, 2], 1), [])
        st.remove_maximal_simplex([0, 1, 2])
        self.assertEqual(st.get_skeleton_tree(2),
                         [([0, 1], 0.0), ([0, 2], 4.0), ([0], 0.0),
                          ([1, 2], 4.0), ([1], 0.0), ([2], 4.0)])
        self.assertFalse(st.find([0, 1, 2]))
        self.assertTrue(st.find([0, 1]))
        self.assertTrue(st.find([0, 2]))
        self.assertTrue(st.find([0]))
        self.assertTrue(st.find([1]))
        self.assertTrue(st.find([2]))

if __name__ == '__main__':
    unittest.main()
