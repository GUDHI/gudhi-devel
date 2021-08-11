""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Hind Montassif

    Copyright (C) 2021 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.datasets.generators import points
from gudhi.datasets.generators import _points

import pytest

def test_sphere():
    assert _points.sphere(n_samples = 10, ambient_dim = 2, radius = 1., sample = 'random').shape == (10, 2)

    with pytest.raises(ValueError):
        _points.sphere(n_samples = 10, ambient_dim = 2, radius = 1., sample = 'other')

def test_torus():
    assert _points.torus(n_samples = 64, dim = 3, sample = 'random').shape == (64, 6)
    assert _points.torus(n_samples = 64, dim = 3, sample = 'grid').shape == (64, 6)
    
    assert _points.torus(n_samples = 10, dim = 4, sample = 'random').shape == (10, 8)
    assert _points.torus(n_samples = 10, dim = 4, sample = 'grid').shape == (1, 8)

    with pytest.raises(ValueError):
        _points.torus(n_samples = 10, dim = 4, sample = 'other')

def test_torus_full_python():
    assert points.torus(n_samples = 64, dim = 3, sample = 'random').shape == (64, 6)
    assert points.torus(n_samples = 64, dim = 3, sample = 'grid').shape == (64, 6)
    
    assert points.torus(n_samples = 10, dim = 4, sample = 'random').shape == (10, 8)
    assert points.torus(n_samples = 10, dim = 4, sample = 'grid').shape == (1, 8)

    with pytest.raises(ValueError):
        points.torus(n_samples = 10, dim = 4, sample = 'other')
