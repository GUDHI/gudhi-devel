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

  def test_rips(self):
    rips_complex = gudhi.SimplexTree(points=[[0,0],[1,0],[0,1],[1,1]],max_dimension=1,max_edge_length=42)

    self.assertEqual(rips_complex.num_simplices(), 10)
    self.assertEqual(rips_complex.num_vertices(), 4)

    self.assertEqual(rips_complex.get_filtered_tree(), [([0], 0.0), ([1], 0.0), ([2], 0.0), ([3], 0.0), ([0, 1], 1.0), ([0, 2], 1.0), ([1, 3], 1.0), ([2, 3], 1.0), ([1, 2], 1.4142135623730951), ([0, 3], 1.4142135623730951)])
    self.assertEqual(rips_complex.get_star_tree([0]), [([0], 0.0), ([0, 1], 1.0), ([0, 2], 1.0), ([0, 3], 1.4142135623730951)])
    self.assertEqual(rips_complex.get_coface_tree([0], 1), [([0, 1], 1.0), ([0, 2], 1.0), ([0, 3], 1.4142135623730951)])

    filtered_rips = gudhi.SimplexTree(points=[[0,0],[1,0],[0,1],[1,1]],max_dimension=1,max_edge_length=1.0)
    self.assertEqual(filtered_rips.num_simplices(), 8)
    self.assertEqual(filtered_rips.num_vertices(), 4)

  def test_mini(self):
    triangle012 = [0,1,2]
    edge03 = [0,3]
    mini_st = gudhi.MiniSimplexTree()
    self.assertTrue(mini_st.insert(triangle012))
    self.assertTrue(mini_st.insert(edge03))
    # FIXME: Remove this line
    mini_st.set_dimension(2);

    edge02 = [0,2]
    self.assertTrue(mini_st.find(edge02))
    self.assertEqual(mini_st.get_coface_tree(edge02, 1), [([0, 1, 2], 0.0)])

    # remove_maximal_simplex test
    self.assertEqual(mini_st.get_coface_tree(triangle012, 1), [])
    mini_st.remove_maximal_simplex(triangle012)
    self.assertTrue(mini_st.find(edge02))
    self.assertFalse(mini_st.find(triangle012))

if __name__ == '__main__':
    unittest.main()