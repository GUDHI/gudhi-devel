import unittest

import gudhi

class TestSimplexTree(unittest.TestCase):

  def test_insertion(self):
    st = gudhi.SimplexTree()

    # insert test
    self.assertTrue(st.insert([0,1]))
    self.assertTrue(st.insert([0,1,2], filtration=4.0))
    self.assertEqual(st.num_simplices(), 7)
    self.assertEqual(st.num_vertices(), 3)

    # find test
    self.assertTrue(st.find([0,1,2]))
    self.assertTrue(st.find([0,1]))
    self.assertTrue(st.find([0,2]))
    self.assertTrue(st.find([0]))
    self.assertTrue(st.find([1]))
    self.assertTrue(st.find([2]))
    self.assertFalse(st.find([3]))
    self.assertFalse(st.find([0,3]))
    self.assertFalse(st.find([1,3]))
    self.assertFalse(st.find([2,3]))

    # filtration test
    st.set_filtration(5.0)
    st.initialize_filtration()
    self.assertEqual(st.get_filtration(), 5.0)
    self.assertEqual(st.filtration([0,1,2]), 4.0)
    self.assertEqual(st.filtration([0,2]), 4.0)
    self.assertEqual(st.filtration([1,2]), 4.0)
    self.assertEqual(st.filtration([2]), 4.0)
    self.assertEqual(st.filtration([0,1]), 0.0)
    self.assertEqual(st.filtration([0]), 0.0)
    self.assertEqual(st.filtration([1]), 0.0)
      
    # skeleton_tree test
    self.assertEqual(st.get_skeleton_tree(2), [([0, 1, 2], 4.0), ([0, 1], 0.0), ([0, 2], 4.0), ([0], 0.0), ([1, 2], 4.0), ([1], 0.0), ([2], 4.0)])
    self.assertEqual(st.get_skeleton_tree(1), [([0, 1], 0.0), ([0, 2], 4.0), ([0], 0.0), ([1, 2], 4.0), ([1], 0.0), ([2], 4.0)])
    self.assertEqual(st.get_skeleton_tree(0), [([0], 0.0), ([1], 0.0), ([2], 4.0)])

    # remove_maximal_simplex test
    self.assertEqual(st.get_coface_tree([0,1,2], 1), [])
    st.remove_maximal_simplex([0,1,2])
    self.assertEqual(st.get_skeleton_tree(2), [([0, 1], 0.0), ([0, 2], 4.0), ([0], 0.0), ([1, 2], 4.0), ([1], 0.0), ([2], 4.0)])
    self.assertFalse(st.find([0,1,2]))
    self.assertTrue(st.find([0,1]))
    self.assertTrue(st.find([0,2]))
    self.assertTrue(st.find([0]))
    self.assertTrue(st.find([1]))
    self.assertTrue(st.find([2]))

if __name__ == '__main__':
    unittest.main()