import unittest

import gudhi

class TestRipsComplex(unittest.TestCase):

  def test_rips(self):
    point_list=[[0,0],[1,0],[0,1],[1,1]]
    rips_complex = gudhi.RipsComplex(points=point_list,max_dimension=1,max_edge_length=42)

    self.assertEqual(rips_complex.num_simplices(), 10)
    self.assertEqual(rips_complex.num_vertices(), 4)

    self.assertEqual(rips_complex.get_filtered_tree(), [([0], 0.0), ([1], 0.0), ([2], 0.0), ([3], 0.0), ([0, 1], 1.0), ([0, 2], 1.0), ([1, 3], 1.0), ([2, 3], 1.0), ([1, 2], 1.4142135623730951), ([0, 3], 1.4142135623730951)])
    self.assertEqual(rips_complex.get_star_tree([0]), [([0], 0.0), ([0, 1], 1.0), ([0, 2], 1.0), ([0, 3], 1.4142135623730951)])
    self.assertEqual(rips_complex.get_coface_tree([0], 1), [([0, 1], 1.0), ([0, 2], 1.0), ([0, 3], 1.4142135623730951)])

    filtered_rips = gudhi.RipsComplex(points=point_list,max_dimension=1,max_edge_length=1.0)

    self.assertEqual(filtered_rips.num_simplices(), 8)
    self.assertEqual(filtered_rips.num_vertices(), 4)

if __name__ == '__main__':
    unittest.main()