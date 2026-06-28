""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2025 Inria

    Modification(s):
      - 2026/06 Vincent Rouvreau: Add tests for F-contiguous and non contiguous numpy arrays
      - YYYY/MM Author: Description of the modification
"""


import numpy as np
from numpy.testing import assert_almost_equal

from gudhi._pers_cub_low_dim_ext import _persistence_on_rectangle_from_top_cells
from gudhi._pers_cub_low_dim_ext import _persistence_on_a_line
from gudhi import CubicalComplex as Cubical


def test_basic_persistence_on_a_line():
    data = np.array([0.0, 1.5, 0.7, 2.8, 3.1, -1.0, 2.0])
    persline = _persistence_on_a_line(data)
    # sort by line, just in case
    persline = persline[persline[:, 0].argsort()]

    cc = Cubical(top_dimensional_cells=data)
    cc.compute_persistence()
    pers0 = cc.persistence_intervals_in_dimension(0)
    # sort by line, just in case
    pers0 = pers0[pers0[:, 0].argsort()]

    assert_almost_equal(pers0, persline)


def test_basic_persistence_on_rectangle_from_top_cells():
    # cf. sklearn.datasets.load_digits().images[0]
    img = np.array(
        [
            [0.0, 0.0, 5.0, 13.0, 9.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 13.0, 15.0, 10.0, 15.0, 5.0, 0.0],
            [0.0, 3.0, 15.0, 2.0, 0.0, 11.0, 8.0, 0.0],
            [0.0, 4.0, 12.0, 0.0, 0.0, 8.0, 8.0, 0.0],
            [0.0, 5.0, 8.0, 0.0, 0.0, 9.0, 8.0, 0.0],
            [0.0, 4.0, 11.0, 0.0, 1.0, 12.0, 7.0, 0.0],
            [0.0, 2.0, 14.0, 5.0, 10.0, 12.0, 0.0, 0.0],
            [0.0, 0.0, 6.0, 13.0, 10.0, 0.0, 0.0, 0.0],
        ]
    )
    cc = Cubical(top_dimensional_cells=img)
    cc.compute_persistence()
    pers0 = cc.persistence_intervals_in_dimension(0)
    pers1 = cc.persistence_intervals_in_dimension(1)

    persrec = _persistence_on_rectangle_from_top_cells(img, 0.0)
    # Check persistence are equal, but maybe not sorted, maybe some duplicates
    assert_almost_equal(np.unique(pers0, axis=0), np.unique(persrec[0], axis=0))
    assert_almost_equal(np.unique(pers1, axis=0), np.unique(persrec[1], axis=0))

def test_contiguity_for_persistence_on_a_line():
    rng = np.random.default_rng(42)
    # Get a 10x10 random values
    base = rng.standard_normal((10))
    persline_c = _persistence_on_a_line(np.ascontiguousarray(base))
    persline_f = _persistence_on_a_line(np.asfortranarray(base))
    assert_almost_equal(persline_c, persline_f)

    base_repeat = np.repeat(base, 2, axis=0)
    persline_r = _persistence_on_a_line(base_repeat[::2])
    assert_almost_equal(persline_c, persline_r)

def test_contiguity_for_persistence_on_rectangle():
    rng = np.random.default_rng(42)
    # Get a 10x10 random values
    base = rng.standard_normal((10, 10))
    persrec_c = _persistence_on_rectangle_from_top_cells(np.ascontiguousarray(base), 0.)
    persrec_f = _persistence_on_rectangle_from_top_cells(np.asfortranarray(base), 0.)
    for idx in range(2):
        assert_almost_equal(persrec_c[idx], persrec_f[idx])

    base_repeat = np.repeat(base, 2, axis=0)
    persrec_r = _persistence_on_rectangle_from_top_cells(base_repeat[::2], 0.)
    for idx in range(2):
        assert_almost_equal(persrec_c[idx], persrec_r[idx])
