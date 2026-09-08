"""This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
Author(s):       Henrique Ennes & Raphaël Tinarrage

Copyright (C) 2016 Inria & Institute of Science and Technology Austria
"""

import numpy as np
from gudhi import liedetect
from gudhi.datasets.linear_orbits import (
    get_canonical_pushforward_algebra,
    sample_from_lie_algebra,
    sample_from_lie_group,
    sample_from_lie_group_rep,
)
from gudhi.liedetect.utils import sample_orbit_from_algebra_su2


def test_sampling():
    group = "torus"
    ambient_dim = 4
    nb_points = 100
    frequency_max = 2
    group_dim = 1
    method = "evenly_spaced"
    span_ambient_space = True
    rep_type = ((1, 2),)
    algebra = [np.array([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 2], [0, 0, -2, 0]])]
    x = np.array([0.5, 0.5, 0.5, 0.5])

    pts_algebra = sample_from_lie_algebra(
        group, rep_type, algebra, nb_points, x=x, method=method
    )

    pts_rep = sample_from_lie_group_rep(group, rep_type, nb_points, method=method)

    assert np.allclose(pts_algebra, pts_rep, atol=1e-6)

    pts_group, _ = sample_from_lie_group(
        group,
        ambient_dim,
        nb_points,
        frequency_max=frequency_max,
        group_dim=group_dim,
        method=method,
        span_ambient_space=span_ambient_space,
    )

    assert pts_group.shape[1] == ambient_dim

    orbit_fitter = liedetect.OrbitFitter(pts_rep)
    nb_neighbors = 10
    orbit_dim = 1

    _ = orbit_fitter.lie_pca(
        nb_neighbors=nb_neighbors, orbit_dim=orbit_dim, method="PCA", verbose=False
    )

    frequency_max = 4

    method = "abelian"
    la, _ = orbit_fitter.closest_algebra(
        group="torus",
        group_dim=1,
        frequency_max=frequency_max,
        method=method,
        verbose=False,
    )

    assert liedetect.are_representations_equivalent("torus", rep_type, la)


def test_reproducible():
    # Seed on SU(2), sample orbit from lie group
    orbit_1, rep_type_1 = sample_from_lie_group(
        group="SU(2)",
        ambient_dim=7,
        nb_points=20,
        seed=123,
    )
    orbit_2, rep_type_2 = sample_from_lie_group(
        group="SU(2)",
        ambient_dim=7,
        nb_points=20,
        seed=123,
    )

    assert rep_type_1 == rep_type_2
    np.testing.assert_allclose(orbit_1, orbit_2)

    # Seed on SU(2), "random_uniform" sample from lie algebra
    rep_type = (3,)
    algebra = get_canonical_pushforward_algebra(
        group="SU(2)",
        rep_type=rep_type,
    )
    x = np.ones(3)
    x /= np.linalg.norm(x)

    orbit_1 = sample_orbit_from_algebra_su2(
        rep_type=rep_type,
        algebra=algebra,
        x=x,
        nb_points=20,
        method="random_uniform",
        seed=123,
    )
    orbit_2 = sample_orbit_from_algebra_su2(
        rep_type=rep_type,
        algebra=algebra,
        x=x,
        nb_points=20,
        method="random_uniform",
        seed=123,
    )

    np.testing.assert_allclose(orbit_1, orbit_2)


def test_optional_base_point_x():
    group = "SO(3)"
    rep_type = (3,)
    ambient_dim = 3
    nb_points = 10

    x_default = np.ones(ambient_dim) / np.sqrt(ambient_dim)
    algebra = get_canonical_pushforward_algebra(group=group, rep_type=rep_type)

    # sample_from_lie_algebra
    orbit_default = sample_from_lie_algebra(
        group=group,
        rep_type=rep_type,
        algebra=algebra,
        nb_points=nb_points,
    )
    orbit_explicit = sample_from_lie_algebra(
        group=group,
        rep_type=rep_type,
        algebra=algebra,
        nb_points=nb_points,
        x=x_default,
    )
    np.testing.assert_allclose(orbit_default, orbit_explicit)

    # sample_from_lie_group_rep
    orbit_default = sample_from_lie_group_rep(
        group=group,
        rep_type=rep_type,
        nb_points=nb_points,
    )
    orbit_explicit = sample_from_lie_group_rep(
        group=group,
        rep_type=rep_type,
        nb_points=nb_points,
        x=x_default,
    )
    np.testing.assert_allclose(orbit_default, orbit_explicit)

    # sample_from_lie_group
    orbit_default, rep_default = sample_from_lie_group(
        group=group,
        ambient_dim=ambient_dim,
        nb_points=nb_points,
        seed=42,
    )
    orbit_explicit, rep_explicit = sample_from_lie_group(
        group=group,
        ambient_dim=ambient_dim,
        nb_points=nb_points,
        x=x_default,
        seed=42,
    )

    assert rep_default == rep_explicit
    np.testing.assert_allclose(orbit_default, orbit_explicit)

    # Test custom base point
    x = np.array([1.0, 0.0, 0.0])

    orbit = sample_from_lie_group_rep(
        group="SO(3)",
        rep_type=(3,),
        nb_points=10,
        x=x,
    )

    # Orthogonal representations preserve the norm of x.
    np.testing.assert_allclose(
        np.linalg.norm(orbit, axis=1),
        np.linalg.norm(x),
    )
