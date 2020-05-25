""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Marc Glisse

    Copyright (C) 2020 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.clustering.tomato import Tomato
import numpy as np
import pytest


def test_tomato_something():
    a = [(1, 2), (1.1, 1.9), (0.9, 1.8), (10, 0), (10.1, 0.05), (10.2, -0.1), (5.4, 0)]
    t = Tomato(metric="euclidean", n_clusters=2, k=4, n_jobs=-1, eps=0.05)
    assert np.array_equal(t.fit_predict(a), [1,1,1,0,0,0,0]) # or with swapped 0 and 1

    t = Tomato(density_type='KDE', r=1, k=4)
    t.fit(a)
    assert np.array_equal(t.leaf_labels_, [1,1,1,0,0,0,0]) # or with swapped 0 and 1

    t = Tomato(graph_type='radius', r=4.7, k=4)
    t.fit(a)
    assert t.max_weight_per_cc_.size == 2
