""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2023 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.datasets.generators import points
from gudhi.sklearn.rips_persistence import RipsPersistence
import numpy as np
import random
import pytest


def test_rips_persistence_of_points_on_a_circle():
    # Let's test with 5 point clouds
    NB_PC = 5
    point_clouds = [points.sphere(n_samples=random.randint(100, 150), ambient_dim=2) for _ in range(NB_PC)]

    rips = RipsPersistence(homology_dimensions=[0, 1], n_jobs=-2)
    diags = rips.fit_transform(point_clouds)
    # list of one array as an input means list of one array as an output
    assert len(diags) == NB_PC
    for idx in range(NB_PC):
        # H0 + H1
        assert len(diags[idx]) == 2
        # H0
        assert len(diags[idx][0]) == len(point_clouds[idx])
        # H1
        assert len(diags[idx][1]) == 1


def test_h1_only_rips_persistence_of_points_on_a_circle():
    rips = RipsPersistence(homology_dimensions=1, n_jobs=-2)
    diags = rips.fit_transform([points.sphere(n_samples=150, ambient_dim=2)])[0]
    assert len(diags) == 1
    assert 0. < diags[0][0] < 0.6
    assert 1. < diags[0][1] < 2.


def test_invalid_input_type():
    rips = RipsPersistence(homology_dimensions=0, input_type='whatsoever')
    with pytest.raises(ValueError):
        rips.fit_transform([points.sphere(n_samples=10, ambient_dim=2)])


def test_distance_matrix_rips_persistence_of_points_on_a_circle():
    try:
        from scipy.spatial.distance import cdist
        pts = points.sphere(n_samples=150, ambient_dim=2)
        distance_matrix = cdist(pts, pts)
        rips = RipsPersistence(homology_dimensions=1, input_type='lower distance matrix')
        diags = rips.fit_transform([distance_matrix])[0]
        assert len(diags) == 1
        assert 0. < diags[0][0] < 0.6
        assert 1. < diags[0][1] < 2.
    except ValueError:
        pass


def test_set_output():
    try:
        import pandas
        NB_PC = 5
        point_clouds = [points.sphere(n_samples=random.randint(100, 150), ambient_dim=2) for _ in range(NB_PC)]

        rips = RipsPersistence(homology_dimensions=[0, 2], n_jobs=-2)
        diags_pandas = rips.set_output(transform="pandas").fit_transform(point_clouds)
        assert 'H0' == diags_pandas.columns[0]
        assert 'H2' == diags_pandas.columns[1]
        assert len(diags_pandas.index) == NB_PC
    except ImportError:
        print("Missing pandas, skipping set_output test")


def test_big():
    # A circle + many isolated points
    n=1000000
    X=np.zeros((n,2))
    X[:,0]=np.arange(n)*100
    X[:24]=points.torus(24,1,'grid')
    # Ripser cannot handle it, have to fall back to SimplexTree
    # Computing the full distance matrix would require too much memory -> kd-tree
    RipsPersistence(range(25), threshold=10).fit_transform([X])
