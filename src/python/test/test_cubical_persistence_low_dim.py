""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       ???

    Copyright (C) 20?? Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__author__ = "???"
__maintainer__ = ""
__copyright__ = "Copyright (C) 20?? Inria"
__license__ = "MIT"


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
    # sort by line, just in case
    pers0 = pers0[pers0[:, 0].argsort()]
    pers1 = pers1[pers1[:, 0].argsort()]

    persrec = _persistence_on_rectangle_from_top_cells(img, 0.0)
    # sort by line, just in case
    persrec0 = persrec[0][persrec[0][:, 0].argsort()]
    persrec1 = persrec[1][persrec[1][:, 0].argsort()]
    assert_almost_equal(pers0, persrec0)
    assert_almost_equal(pers1, persrec1)
