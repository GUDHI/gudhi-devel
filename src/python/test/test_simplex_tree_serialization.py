""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2023 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""


import numpy as np
import pickle
import pytest

from gudhi import SimplexTree


def test_pickle_simplex_tree():
    st = SimplexTree.create_from_array(np.random.rand(10, 10))
    for dim in [1, 2, 3]:
        st.expansion(dim)
        with open("stree.pkl", "wb") as f:
            pickle.dump(st, f)
        with open("stree.pkl", "rb") as f:
            st_copy = pickle.load(f)
        assert st == st_copy


def test_simplex_tree_serialization_copy():
    st = SimplexTree()
    st.insert([1, 2, 3], 0.5)
    # compute persistence only on the original
    st.compute_persistence()

    st_copy = st.copy()
    st_data = pickle.dumps(st)
    # Modify the original before reloading it from its dump
    st.insert([4, 5], 2.0)
    st = pickle.loads(st_data)
    assert st == st_copy
