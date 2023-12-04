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
from tempfile import NamedTemporaryFile

def test_off_rw():
    for dim in range(2, 6):
        X = np.random.rand(123, dim)
        gd.write_points_to_off_file("rand.off", X)
        Y = gd.read_points_from_off_file("rand.off")
        assert Y == pytest.approx(X)

def test_human_off():
    pts = gd.read_points_from_off_file("human.off")
    # Should not try to read faces
    assert pts.shape == (4706, 3)

def test_invalid_off_file():
    name = NamedTemporaryFile().name
    with open(name, "w") as f:
        f.write("nOFL\n1000 200 2\n1.2 1.3 1.4")
    with pytest.raises(ValueError):
        gd.read_points_from_off_file(name)
    # Try to open a non-existing file - a new temp file name should not exist
    with pytest.raises(FileNotFoundError):
        gd.read_points_from_off_file(NamedTemporaryFile().name)
