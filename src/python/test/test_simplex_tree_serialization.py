""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2023 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi import SimplexTree
import numpy as np
import pickle
import pytest

def test_pickle_simplex_tree():
    st = SimplexTree.create_from_array(np.random.rand(10, 10))
    for dim in [1, 2, 3]:
        st.expansion(dim)
        with open('stree.pkl','wb') as f:
            pickle.dump(st, f)
        with open('stree.pkl','rb') as f:
            st_copy = pickle.load(f)
        assert st == st_copy
