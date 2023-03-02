""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2021 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.sklearn.cubical_persistence import CubicalPersistence
import numpy as np
from sklearn import datasets

CUBICAL_PERSISTENCE_H0_IMG0 = np.array([[0.0, 6.0], [0.0, 8.0], [0.0, np.inf]])


def test_simple_constructor_from_top_cells_list():
    digits = datasets.load_digits().images[:10]
    cp = CubicalPersistence(homology_dimensions=0, n_jobs=-2)

    diags = cp.fit_transform(digits)
    assert len(diags) == 10
    np.testing.assert_array_equal(diags[0], CUBICAL_PERSISTENCE_H0_IMG0)

    cp = CubicalPersistence(homology_dimensions=[0, 1], n_jobs=-1)
    diagsH0H1 = cp.fit_transform(digits)
    assert len(diagsH0H1) == 10
    for idx in range(10):
        np.testing.assert_array_equal(diags[idx], diagsH0H1[idx][0])


def test_simple_constructor_from_flattened_cells():
    cells = datasets.load_digits().images[0]
    # Not squared (extended) flatten cells
    flat_cells = np.hstack((cells, np.zeros((cells.shape[0], 2)))).flatten()

    cp = CubicalPersistence(homology_dimensions=0, newshape=[-1, 8, 10])
    diags = cp.fit_transform([flat_cells])

    np.testing.assert_array_equal(diags[0], CUBICAL_PERSISTENCE_H0_IMG0)

    # Not squared (extended) non-flatten cells
    cells = np.hstack((cells, np.zeros((cells.shape[0], 2))))

    # The aim of this second part of the test is to resize even if not mandatory
    cp = CubicalPersistence(homology_dimensions=0, newshape=[-1, 8, 10])
    diags = cp.fit_transform([cells])

    np.testing.assert_array_equal(diags[0], CUBICAL_PERSISTENCE_H0_IMG0)


def test_1d():
    a = np.array([2, 4, 3, 5])
    r = np.array([[3, 4], [2, np.inf]])
    np.testing.assert_array_equal(r, CubicalPersistence(0).fit_transform([a])[0])
    np.testing.assert_array_equal(r, CubicalPersistence(0).fit_transform([a.astype("float32")[::-1]])[0])
    np.testing.assert_array_equal(r, CubicalPersistence([0]).fit_transform([a.astype("float64")])[0][0])
