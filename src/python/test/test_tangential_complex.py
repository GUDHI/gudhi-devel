""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi import TangentialComplex, SimplexTree
import pytest

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"


def test_tangential():
    point_list = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
    tc = TangentialComplex(intrisic_dim=1, points=point_list)
    assert tc.__is_defined() == True
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
    assert st.__is_defined() == True
    assert st.__is_persistence_defined() == False

    assert st.num_simplices() == 6
    assert st.num_vertices() == 4

    filtration_generator = st.get_filtration()
    assert(next(filtration_generator) == ([0], 0.0))
    assert(next(filtration_generator) == ([1], 0.0))
    assert(next(filtration_generator) == ([2], 0.0))
    assert(next(filtration_generator) == ([0, 2], 0.0))
    assert(next(filtration_generator) == ([3], 0.0))
    assert(next(filtration_generator) == ([1, 3], 0.0))
    with pytest.raises(StopIteration):
        next(filtration_generator)

    assert st.get_cofaces([0], 1) == [([0, 2], 0.0)]

    assert point_list[0] == tc.get_point(0)
    assert point_list[1] == tc.get_point(1)
    assert point_list[2] == tc.get_point(2)
    assert point_list[3] == tc.get_point(3)
    assert tc.get_point(4) == []
    assert tc.get_point(125) == []
