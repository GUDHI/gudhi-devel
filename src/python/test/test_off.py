""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Marc Glisse

    Copyright (C) 2022 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

import gudhi as gd
import numpy as np
import pytest


def test_off_rw():
    for dim in range(2, 6):
        X = np.random.rand(123, dim)
        gd.write_points_to_off_file("rand.off", X)
        Y = gd.read_points_from_off_file("rand.off")
        assert Y == pytest.approx(X)
