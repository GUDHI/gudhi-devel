import unittest

import gudhi


class TestAlphaComplex(unittest.TestCase):

    def test_infinite_alpha(self):
        point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
        alpha_complex = gudhi.AlphaComplex(points=point_list)

        self.assertEqual(alpha_complex.num_simplices(), 11)
        self.assertEqual(alpha_complex.num_vertices(), 4)

        self.assertEqual(alpha_complex.get_filtered_tree(),
                         [([0], 0.0), ([1], 0.0), ([2], 0.0), ([3], 0.0),
                          ([0, 1], 0.25), ([0, 2], 0.25), ([1, 3], 0.25),
                          ([2, 3], 0.25), ([1, 2], 0.5), ([0, 1, 2], 0.5),
                          ([1, 2, 3], 0.5)])
        self.assertEqual(alpha_complex.get_star_tree([0]),
                         [([0], 0.0), ([0, 1], 0.25), ([0, 1, 2], 0.5),
                         ([0, 2], 0.25)])
        self.assertEqual(alpha_complex.get_coface_tree([0], 1),
                         [([0, 1], 0.25), ([0, 2], 0.25)])

        self.assertEqual(point_list[0], alpha_complex.get_point(0))
        self.assertEqual(point_list[1], alpha_complex.get_point(1))
        self.assertEqual(point_list[2], alpha_complex.get_point(2))
        self.assertEqual(point_list[3], alpha_complex.get_point(3))
        self.assertEqual([], alpha_complex.get_point(4))
        self.assertEqual([], alpha_complex.get_point(125))

    def test_filtered_alpha(self):
        point_list = [[0, 0], [1, 0], [0, 1], [1, 1]]
        filtered_alpha = gudhi.AlphaComplex(points=point_list,
                                            max_alpha_square=0.25)

        self.assertEqual(filtered_alpha.num_simplices(), 8)
        self.assertEqual(filtered_alpha.num_vertices(), 4)

        self.assertEqual(point_list[0], filtered_alpha.get_point(0))
        self.assertEqual(point_list[1], filtered_alpha.get_point(1))
        self.assertEqual(point_list[2], filtered_alpha.get_point(2))
        self.assertEqual(point_list[3], filtered_alpha.get_point(3))
        self.assertEqual([], filtered_alpha.get_point(4))
        self.assertEqual([], filtered_alpha.get_point(125))

        self.assertEqual(filtered_alpha.get_filtered_tree(),
                         [([0], 0.0), ([1], 0.0), ([2], 0.0), ([3], 0.0),
                          ([0, 1], 0.25), ([0, 2], 0.25), ([1, 3], 0.25),
                          ([2, 3], 0.25)])
        self.assertEqual(filtered_alpha.get_star_tree([0]),
                         [([0], 0.0), ([0, 1], 0.25), ([0, 2], 0.25)])
        self.assertEqual(filtered_alpha.get_coface_tree([0], 1),
                         [([0, 1], 0.25), ([0, 2], 0.25)])

if __name__ == '__main__':
    unittest.main()
