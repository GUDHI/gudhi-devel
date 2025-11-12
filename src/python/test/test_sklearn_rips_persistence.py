""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2023 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""


from scipy.sparse import coo_matrix
from scipy.spatial.distance import cdist
import numpy as np
import random
import pytest

from gudhi.datasets.generators import points
from gudhi.sklearn import RipsPersistence
from gudhi import RipsComplex, SimplexTree
from gudhi._ripser_ext import _lower, _full, _sparse, _lower_to_coo, _lower_cone_radius
from gudhi import bottleneck_distance


def test_rips_persistence_of_points_on_a_circle():
    # Let's test with 5 point clouds
    NB_PC = 5
    point_clouds = [
        points.sphere(n_samples=random.randint(100, 150), ambient_dim=2) for _ in range(NB_PC)
    ]

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
        point_clouds = [
            points.sphere(n_samples=random.randint(100, 150), ambient_dim=2)
            for _ in range(NB_PC)
        ]

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


def cmp_rips(point_cloud):
    primes = [2, 3, 11, 17, 29]  # Small list so 2 is often selected
    field = random.choice(primes)
    print(f"random prime = {field}")

    print(f"nb points = {len(point_cloud)}, dim = {point_cloud.shape[1]}")
    dists = cdist(point_cloud, point_cloud)

    # Check cone radius
    cr = dists.max(-1).min()
    assert cr == _lower_cone_radius(dists)
    assert cr < 2.0

    ## Compute with the SimplexTree
    stree = RipsComplex(distance_matrix=dists).create_simplex_tree(max_dimension=2)
    stree.compute_persistence(homology_coeff_field=field, persistence_dim_max=True)
    dgm0 = stree.persistence_intervals_in_dimension(0)
    dgm1 = stree.persistence_intervals_in_dimension(1)

    # Compute with Ripser
    dgm = _full(
        dists, max_dimension=1, max_edge_length=float("inf"), homology_coeff_field=field
    )
    # The order of the intervals may differ, so we cannot compare the arrays with np.testing.assert_almost_equal
    assert bottleneck_distance(dgm0, dgm[0]) < 1e-8
    assert bottleneck_distance(dgm1, dgm[1]) < 1e-8

    dgm = _lower(
        dists, max_dimension=1, max_edge_length=float("inf"), homology_coeff_field=field
    )
    assert bottleneck_distance(dgm0, dgm[0]) < 1e-8
    assert bottleneck_distance(dgm1, dgm[1]) < 1e-8

    # Convert to coo matrix
    n = len(point_cloud)
    dists_copy = np.array(dists, copy=True)
    dists_copy[np.triu_indices_from(dists_copy)] = 0  # Keep only the lower entries
    dists_sparse = coo_matrix(dists_copy)
    s1 = _lower_to_coo(dists, float("inf"))
    s2 = _lower_to_coo(dists_copy, float("inf"))
    s1 = coo_matrix((s1[2], (s1[0], s1[1])), shape=(n, n))
    s2 = coo_matrix((s2[2], (s2[0], s2[1])), shape=(n, n))
    assert np.array_equal(s1.toarray(), s2.toarray())
    assert np.array_equal(dists_sparse.toarray(), s1.toarray())
    ## Compute with the SimplexTree
    stree = SimplexTree()
    stree.insert_batch(np.arange(n).reshape(1, -1), np.zeros(n))
    stree.insert_edges_from_coo_matrix(dists_sparse)
    stree.expansion(2)
    stree.compute_persistence(homology_coeff_field=field, persistence_dim_max=True)
    sp_dgm0 = stree.persistence_intervals_in_dimension(0)
    sp_dgm1 = stree.persistence_intervals_in_dimension(1)

    # Compute with Ripser
    sp_dgm = _sparse(
        dists_sparse.row,
        dists_sparse.col,
        dists_sparse.data,
        dists_sparse.shape[0],
        max_dimension=1,
        max_edge_length=float("inf"),
        homology_coeff_field=field,
    )
    assert bottleneck_distance(sp_dgm0, sp_dgm[0]) < 1e-8
    assert bottleneck_distance(sp_dgm1, sp_dgm[1]) < 1e-8
    assert bottleneck_distance(sp_dgm0, dgm0) < 1e-8
    assert bottleneck_distance(sp_dgm1, dgm1) < 1e-8


def test_ripser_interfaces():
    cmp_rips(points.sphere(n_samples=random.randint(100, 150), ambient_dim=2))
    cmp_rips(np.random.rand(random.randint(100, 150), random.randint(2, 3)))
