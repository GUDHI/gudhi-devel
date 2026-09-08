"""This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
Author(s):       Henrique Ennes & Raphaël Tinarrage

Copyright (C) 2016 Inria & Institute of Science and Technology Austria
"""

import pytest
from gudhi import liedetect
from gudhi.datasets import linear_orbits
from gudhi.liedetect import OrbitFitter

"""
Fixes data set.
"""

rep_type = ((3),)


@pytest.fixture
def fix_pts():
    """
    Builds orbit of weights 3 in R^3.
    """
    group = "SO(3)"
    seed = 42
    nb_points = 10000

    return linear_orbits.sample_from_lie_group_rep(
        group=group,
        rep_type=rep_type,
        nb_points=nb_points,
        seed=seed,
    )


"""
Test functions
"""


def test_lie_pca(fix_pts):
    orbit_fitter = OrbitFitter(fix_pts)
    nb_neighbors = 100
    orbit_dim = 3

    _ = orbit_fitter.lie_pca(
        nb_neighbors=nb_neighbors, orbit_dim=orbit_dim, verbose=False
    )
    vals = orbit_fitter.print_lie_pca_eigenvalues(return_vals=True)
    assert abs(vals[0]) < 1e-2
    assert abs(vals[1]) < 1e-2
    assert abs(vals[2]) < 1e-2

    for i in vals[3:]:
        assert abs(i) > 1e-2


def test_project_lie_algebra(fix_pts):
    orbit_fitter = OrbitFitter(fix_pts)
    nb_neighbors = 100
    orbit_dim = 3

    _ = orbit_fitter.lie_pca(
        nb_neighbors=nb_neighbors, orbit_dim=orbit_dim, method="PCA", verbose=False
    )

    method = "bottom_lie_pca"
    la, _ = orbit_fitter.closest_algebra(
        group="SO(3)", group_dim=1, method=method, verbose=False
    )
    assert liedetect.are_representations_equivalent("SO(3)", la, (3,))


def test_sample_orbit(fix_pts):
    orbit_fitter = OrbitFitter(fix_pts)
    nb_neighbors = 100
    orbit_dim = 3

    _ = orbit_fitter.lie_pca(
        nb_neighbors=nb_neighbors, orbit_dim=orbit_dim, method="PCA", verbose=False
    )
    method = "bottom_lie_pca"
    la, _ = orbit_fitter.closest_algebra(group="SO(3)", method=method, verbose=False)

    _ = orbit_fitter.sample_orbit(nb_points=10000)

    assert min(orbit_fitter.hausdorff_distances_) < 1e-1
