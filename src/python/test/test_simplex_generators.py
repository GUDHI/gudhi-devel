""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Marc Glisse

    Copyright (C) 2020 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

import gudhi
import numpy as np


def test_flag_generators():
    pts = np.array([[0, 0], [0, 1.01], [1, 0], [1.02, 1.03], [100, 0], [100, 3.01], [103, 0], [103.02, 3.03]])
    r = gudhi.RipsComplex(pts, max_edge_length=4)
    st = r.create_simplex_tree(max_dimension=50)
    st.persistence()
    g = st.flag_persistence_generators()
    assert np.array_equal(g[0], [[2, 2, 0], [1, 1, 0], [3, 3, 1], [6, 6, 4], [5, 5, 4], [7, 7, 5]])
    assert len(g[1]) == 1
    assert np.array_equal(g[1][0], [[3, 2, 2, 1]])
    assert np.array_equal(g[2], [0, 4])
    assert len(g[3]) == 1
    assert np.array_equal(g[3][0], [[7, 6]])


def test_lower_star_generators():
    st = gudhi.SimplexTree()
    st.insert([0, 1, 2], -10)
    st.insert([0, 3], -10)
    st.insert([1, 3], -10)
    st.assign_filtration([2], -1)
    st.assign_filtration([3], 0)
    st.assign_filtration([0], 1)
    st.assign_filtration([1], 2)
    st.make_filtration_non_decreasing()
    st.persistence(min_persistence=-1)
    g = st.lower_star_persistence_generators(min_persistence=-1)
    assert len(g[0]) == 2
    assert np.array_equal(g[0][0], [[0, 0], [3, 0], [1, 1]])
    assert np.array_equal(g[0][1], [[1, 1]])
    assert len(g[1]) == 2
    assert np.array_equal(g[1][0], [2])
    assert np.array_equal(g[1][1], [1])


def test_empty():
    st = gudhi.SimplexTree()
    st.persistence()
    assert st.lower_star_persistence_generators() == ([], [])
    g = st.flag_persistence_generators()
    assert np.array_equal(g[0], np.empty((0, 3)))
    assert g[1] == []
    assert np.array_equal(g[2], [])
    assert g[3] == []
