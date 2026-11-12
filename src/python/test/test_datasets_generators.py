""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Hind Montassif

    Copyright (C) 2021 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

import pytest
from numpy.testing import assert_almost_equal

from gudhi.datasets.generators import points
from gudhi.random import set_seed


def test_sphere():
    assert points.sphere(n_samples=10, ambient_dim=2, radius=1.0, sample="random").shape == (
        10,
        2,
    )

    with pytest.raises(ValueError):
        points.sphere(n_samples=10, ambient_dim=2, radius=1.0, sample="other")


def _basic_torus(impl):
    assert impl(n_samples=64, dim=3, sample="random").shape == (64, 6)
    assert impl(n_samples=64, dim=3, sample="grid").shape == (64, 6)

    assert impl(n_samples=10, dim=4, sample="random").shape == (10, 8)

    # Here 1**dim < n_samples < 2**dim, the output shape is therefore (1, 2*dim) = (1, 8), where shape[0] is rounded down to the closest perfect 'dim'th power
    assert impl(n_samples=10, dim=4, sample="grid").shape == (1, 8)

    with pytest.raises(ValueError):
        impl(n_samples=10, dim=4, sample="other")


def test_torus():
    for torus_impl in [points.torus, points.ctorus]:
        _basic_torus(torus_impl)
    # Check that the two versions (torus and ctorus) generate the same output
    assert (
        points.ctorus(n_samples=64, dim=3, sample="random").all()
        == points.torus(n_samples=64, dim=3, sample="random").all()
    )
    assert (
        points.ctorus(n_samples=64, dim=3, sample="grid").all()
        == points.torus(n_samples=64, dim=3, sample="grid").all()
    )
    assert (
        points.ctorus(n_samples=10, dim=3, sample="grid").all()
        == points.torus(n_samples=10, dim=3, sample="grid").all()
    )


def test_reproductibility_sphere():
    set_seed(42)
    rand_pts1 = points.sphere(n_samples=20, ambient_dim=4)
    set_seed(42)
    rand_pts2 = points.sphere(n_samples=20, ambient_dim=4)
    assert_almost_equal(rand_pts1, rand_pts2)


def test_reproductibility_ctorus():
    for sample in ["grid", "random"]:
        set_seed(42)
        rand_pts1 = points.ctorus(n_samples=20, dim=4, sample=sample)
        set_seed(42)
        rand_pts2 = points.ctorus(n_samples=20, dim=4, sample=sample)
        assert_almost_equal(rand_pts1, rand_pts2)


def test_reproductibility_torus_embedded_r3():
    for sample in ["grid", "random"]:
        set_seed(42)
        rand_pts1 = points.torus_embedded_r3(n_samples=20, major_radius=4., minor_radius=1., sample=sample)
        set_seed(42)
        rand_pts2 = points.torus_embedded_r3(n_samples=20, major_radius=4., minor_radius=1., sample=sample)
        assert_almost_equal(rand_pts1, rand_pts2)


def test_torus_embedded_r3():
    random_torus = points.torus_embedded_r3(n_samples=64, sample="random")
    grid_torus = points.torus_embedded_r3(n_samples=64, sample="grid")
    assert random_torus.shape == (64, 3)
    assert grid_torus.shape   == (64, 3)

    for torus in [random_torus, grid_torus]:
        for point in torus:
            # Test the point is inside a sphere of radius 1.5 (major radius = 1. + minor radius = 0.5)
            assert (point[0]*point[0] + point[1]*point[1] + point[2]*point[2]) <= 2.25
    
    with pytest.raises(ValueError):
        points.torus_embedded_r3(n_samples=10, sample="other")
