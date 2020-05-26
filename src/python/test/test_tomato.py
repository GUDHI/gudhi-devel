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
import matplotlib.pyplot as plt

# Disable graphics for testing purposes
plt.show = lambda: None


def test_tomato_something():
    a = [(1, 2), (1.1, 1.9), (0.9, 1.8), (10, 0), (10.1, 0.05), (10.2, -0.1), (5.4, 0)]
    t = Tomato(metric="euclidean", n_clusters=2, k=4, n_jobs=-1, eps=0.05)
    assert np.array_equal(t.fit_predict(a), [1, 1, 1, 0, 0, 0, 0])  # or with swapped 0 and 1
    assert np.array_equal(t.children_, [[0, 1]])

    t = Tomato(density_type="KDE", r=1, k=4)
    t.fit(a)
    assert np.array_equal(t.leaf_labels_, [1, 1, 1, 0, 0, 0, 0])  # or with swapped 0 and 1
    assert t.n_clusters_ == 2
    t.merge_threshold_ = 10
    assert t.n_clusters_ == 1
    assert (t.labels_ == 0).all()

    t = Tomato(metric="euclidean", graph_type="radius", r=4.7, k=4)
    t.fit(a)
    assert t.max_weight_per_cc_.size == 2
    assert np.array_equal(t.neighbors_, [[0, 1, 2], [0, 1, 2], [0, 1, 2], [3, 4, 5, 6], [3, 4, 5], [3, 4, 5], [3, 6]])
    t.plot_diagram()

    t = Tomato(graph_type="radius", r=4.7, k=4, symmetrize_graph=True)
    t.fit(a)
    assert t.max_weight_per_cc_.size == 2
    assert [set(i) for i in t.neighbors_] == [{1, 2}, {0, 2}, {0, 1}, {4, 5, 6}, {3, 5}, {3, 4}, {3}]

    t = Tomato(n_clusters=2, k=4, symmetrize_graph=True)
    t.fit(a)
    assert [set(i) for i in t.neighbors_] == [
        {1, 2, 6},
        {0, 2, 6},
        {0, 1, 6},
        {4, 5, 6},
        {3, 5, 6},
        {3, 4, 6},
        {0, 1, 2, 3, 4, 5},
    ]
    t.plot_diagram()

    t = Tomato(k=6, metric="manhattan")
    t.fit(a)
    assert t.diagram_.size == 0
    assert t.max_weight_per_cc_.size == 1
    t.plot_diagram()
