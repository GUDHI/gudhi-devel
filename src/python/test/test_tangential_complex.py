""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - 2025/04 Hannah Schreiber: Add tests to verify possibility of tensor input
      - YYYY/MM Author: Description of the modification
"""

__author__ = "Vincent Rouvreau"
__maintainer__ = "Hannah Schreiber"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "GPL v3"


from gudhi import TangentialComplex, SimplexTree


def test_tangential():
    point_list = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
    tc = TangentialComplex(intrisic_dim=1, points=point_list)
    assert tc._is_defined() == True
    assert tc.num_vertices() == 4
    assert tc.num_simplices() == 0
    assert tc.num_inconsistent_simplices() == 0
    assert tc.num_inconsistent_stars() == 0

    tc.compute_tangential_complex()
    assert tc.num_vertices() == 4
    assert tc.num_simplices() == 4
    assert tc.num_inconsistent_simplices() == 0
    assert tc.num_inconsistent_stars() == 0

    st = tc.create_simplex_tree()
    assert st._is_defined() == True
    assert st._is_persistence_defined() == False

    assert st.num_simplices() == 6
    assert st.num_vertices() == 4

    assert list(st.get_filtration()) == [
        ([0], 0.0),
        ([1], 0.0),
        ([2], 0.0),
        ([0, 2], 0.0),
        ([3], 0.0),
        ([1, 3], 0.0),
    ]

    assert st.get_cofaces([0], 1) == [([0, 2], 0.0)]

    assert point_list[0] == tc.get_point(0)
    assert point_list[1] == tc.get_point(1)
    assert point_list[2] == tc.get_point(2)
    assert point_list[3] == tc.get_point(3)
    assert tc.get_point(4) == []
    assert tc.get_point(125) == []

def test_tensors():
    try:
        import torch

        point_list = (torch.rand((5, 2))).requires_grad_()
        cplex = TangentialComplex(intrisic_dim=1, points=point_list)
    except ImportError:
        pass

    try:
        import tensorflow as tf

        point_list = tf.random.uniform(shape=[5, 2])
        cplex = TangentialComplex(intrisic_dim=1, points=point_list)
    except ImportError:
        pass
