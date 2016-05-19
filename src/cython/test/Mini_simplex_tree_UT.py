import unittest

import gudhi

class TestMiniSimplexTree(unittest.TestCase):

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