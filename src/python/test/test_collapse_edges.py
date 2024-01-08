""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Marc Glisse

    Copyright (C) 2023 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.flag_filtration.edge_collapse import reduce_graph as collapse_edges
import numpy as np
from scipy.sparse import coo_matrix
import pytest


def test_collapse_empty():
    x = coo_matrix(([], ([], [])), (3, 3))
    xo = collapse_edges(x)
    assert xo.shape == (3, 3)
    assert len(xo.data) == 0


def test_collapse():
    x = coo_matrix(([0.1, 0.2, 0.3, 0.4, 0.6], ([0, 1, 2, 3, 1], [1, 2, 3, 0, 3])), (6, 6))
    xo = collapse_edges(x)
    assert xo.shape == (6, 6)
    assert len(xo.data) == 5

    x = coo_matrix(([0.1, 0.2, 0.3, 0.4, -0.6], ([0, 1, 2, 3, 1], [1, 2, 3, 0, 3])), (4, 4))
    xo = collapse_edges(x)
    assert xo.shape == (4, 4)
    assert len(xo.data) == 3
    assert xo.dtype == np.dtype("float64")

    x = coo_matrix((np.array([0.1], dtype="float32"), ([0], [1])), (2, 2))
    xo = collapse_edges(x)
    assert xo.dtype == np.dtype("float32")
