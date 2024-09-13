""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2023 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.datasets.generators import points
from gudhi.sklearn.rips_persistence import RipsPersistence
from gudhi import RipsComplex, SimplexTree
from gudhi._ripser import _lower, _full, _sparse, _lower_to_coo, _lower_cone_radius
from scipy.sparse import coo_matrix
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import cKDTree
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
    assert 0.0 < diags[0][0] < 0.6
    assert 1.0 < diags[0][1] < 2.0


def test_invalid_input_type():
    rips = RipsPersistence(homology_dimensions=0, input_type="whatsoever")
    with pytest.raises(ValueError):
        rips.fit_transform([points.sphere(n_samples=10, ambient_dim=2)])


def test_distance_matrix_rips_persistence_of_points_on_a_circle():
    try:
        from scipy.spatial.distance import cdist

        pts = points.sphere(n_samples=150, ambient_dim=2)
        distance_matrix = cdist(pts, pts)
        rips = RipsPersistence(homology_dimensions=1, input_type="lower distance matrix")
        diags = rips.fit_transform([distance_matrix])[0]
        assert len(diags) == 1
        assert 0.0 < diags[0][0] < 0.6
        assert 1.0 < diags[0][1] < 2.0
    except ValueError:
        pass


def test_set_output():
    try:
        import pandas

        NB_PC = 5
        point_clouds = [points.sphere(n_samples=random.randint(100, 150), ambient_dim=2) for _ in range(NB_PC)]

        rips = RipsPersistence(homology_dimensions=[0, 2], n_jobs=-2)
        diags_pandas = rips.set_output(transform="pandas").fit_transform(point_clouds)
        assert "H0" == diags_pandas.columns[0]
        assert "H2" == diags_pandas.columns[1]
        assert len(diags_pandas.index) == NB_PC
    except ImportError:
        print("Missing pandas, skipping set_output test")


def test_big():
    # A circle + many isolated points
    n = 1000000
    X = np.zeros((n, 2))
    X[:, 0] = np.arange(n) * 100
    X[:24] = points.torus(24, 1, "grid")
    # Ripser cannot handle it, have to fall back to SimplexTree
    # Computing the full distance matrix would require too much memory -> kd-tree
    RipsPersistence(range(25), threshold=10).fit_transform([X])


def test_ripser_interfaces():
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    random_prime = primes[random.randint(0, 9)]
    print(f"random prime = {random_prime}")

    point_cloud = points.sphere(n_samples=random.randint(100, 150), ambient_dim=2)
    print(f"nb points = {len(point_cloud)}")
    inp = squareform(pdist(point_cloud))

    # Check cone radius
    assert inp.max(-1).min() < 2.0
    assert _lower_cone_radius(inp) < 2.0

    ## As there is no easy way to force the use of SimplexTree, let's build it
    stree = RipsComplex(distance_matrix=inp).create_simplex_tree(max_dimension=2)
    stree.compute_persistence(homology_coeff_field=random_prime, persistence_dim_max=True)
    dgm0 = stree.persistence_intervals_in_dimension(0)
    dgm1 = stree.persistence_intervals_in_dimension(1)

    dgm = _full(inp, max_dimension=2, max_edge_length=float("inf"), homology_coeff_field=random_prime)
    np.testing.assert_almost_equal(dgm0, dgm[0])
    np.testing.assert_almost_equal(dgm1, dgm[1])

    dgm = _lower(inp, max_dimension=2, max_edge_length=float("inf"), homology_coeff_field=random_prime)
    np.testing.assert_almost_equal(dgm0, dgm[0])
    np.testing.assert_almost_equal(dgm1, dgm[1])

    # From a coo matrix
    n = len(point_cloud)
    tree = cKDTree(point_cloud)
    pairs = tree.query_pairs(r=float("inf"), output_type="ndarray")
    data = np.ravel(np.linalg.norm(np.diff(point_cloud[pairs], axis=1), axis=-1))
    inp = coo_matrix((data, (pairs[:, 0], pairs[:, 1])), shape=(n,) * 2)
    ## As there is no easy way to force the use of SimplexTree, let's build it
    stree = SimplexTree()
    stree.insert_batch(np.arange(n).reshape(1, -1), np.zeros(n))
    stree.insert_edges_from_coo_matrix(inp)
    stree.expansion(2)
    stree.compute_persistence(homology_coeff_field=random_prime, persistence_dim_max=True)
    dgm0 = stree.persistence_intervals_in_dimension(0)
    dgm1 = stree.persistence_intervals_in_dimension(1)

    dgm = _sparse(
        inp.row,
        inp.col,
        inp.data,
        inp.shape[0],
        max_dimension=2,
        max_edge_length=float("inf"),
        homology_coeff_field=random_prime,
    )
    np.testing.assert_almost_equal(dgm0, dgm[0])
    np.testing.assert_almost_equal(dgm1, dgm[1])
