""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Hind Montassif

    Copyright (C) 2021 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.datasets.generators import points

import pytest

def test_sphere():
    assert points.sphere(n_samples = 10, ambient_dim = 2, radius = 1., sample = 'random').shape == (10, 2)

    with pytest.raises(ValueError):
        points.sphere(n_samples = 10, ambient_dim = 2, radius = 1., sample = 'other')

def _basic_torus(impl):
    assert impl(n_samples = 64, dim = 3, sample = 'random').shape == (64, 6)
    assert impl(n_samples = 64, dim = 3, sample = 'grid').shape == (64, 6)

    assert impl(n_samples = 10, dim = 4, sample = 'random').shape == (10, 8)

    # Here 1**dim < n_samples < 2**dim, the output shape is therefore (1, 2*dim) = (1, 8), where shape[0] is rounded down to the closest perfect 'dim'th power
    assert impl(n_samples = 10, dim = 4, sample = 'grid').shape == (1, 8)

    with pytest.raises(ValueError):
        impl(n_samples = 10, dim = 4, sample = 'other')

def test_torus():
    for torus_impl in [points.torus, points.ctorus]:
        _basic_torus(torus_impl)
    # Check that the two versions (torus and ctorus) generate the same output
    assert points.ctorus(n_samples = 64, dim = 3, sample = 'random').all() == points.torus(n_samples = 64, dim = 3, sample = 'random').all()
    assert points.ctorus(n_samples = 64, dim = 3, sample = 'grid').all() == points.torus(n_samples = 64, dim = 3, sample = 'grid').all()
    assert points.ctorus(n_samples = 10, dim = 3, sample = 'grid').all() == points.torus(n_samples = 10, dim = 3, sample = 'grid').all()
