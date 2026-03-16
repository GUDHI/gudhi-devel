""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2025 Inria

    Modification(s):
      - 2026/03 Vincent Rouvreau: Rewrite H1 tests and only check lifetimes > 1e-11
      - YYYY/MM Author: Description of the modification
"""


import numpy as np
from numpy.testing import assert_almost_equal
import pytest
import math

from gudhi.datasets.generators import points
from gudhi.sklearn import CechPersistence, WeightedCechPersistence


def test_cech_persistence_of_points_on_a_circle():
    # Let's test with 5 point clouds with N points (100 <= N <= 150) on 5 circles
    NB_PC = 5
    n_samples_array = np.random.randint(100, high=150, size=NB_PC)
    pts = [points.sphere(n_samples=n_samples, ambient_dim=2) for n_samples in n_samples_array]

    for precision in ["safe", "fast", "exact"]:
        diags = CechPersistence(homology_dimensions=[0, 1], precision=precision, n_jobs=-2).fit_transform(pts)
        # list of one array as an input means list of one array as an output
        assert len(diags) == NB_PC
        for idx in range(NB_PC):
            # H0 + H1
            assert len(diags[idx]) == 2
            # H0
            assert len(diags[idx][0]) == len(pts[idx])
            # H1
            h1 = diags[idx][1]
            lifetimes = np.abs(h1[:, 1] - h1[:, 0])
            # Degenerate case with points on a circle, some 1e-12 H1 lifetimes can appear
            epsilon = 1e-11
            assert np.sum(lifetimes > epsilon) == 1

        # Same test but for the weighted version
        wpts = [np.append(pts[idx], np.zeros((len(pts[idx]), 1)), axis=1) for idx in range(NB_PC)]
        diags = WeightedCechPersistence(homology_dimensions=[0, 1], precision=precision, n_jobs=-2).fit_transform(wpts)
        assert len(diags) == NB_PC
        for idx in range(NB_PC):
            # H0 + H1
            assert len(diags[idx]) == 2
            # H0
            assert len(diags[idx][0]) == len(wpts[idx])
            # H1 - can have some [1., 1.] - Maybe points on a circle is not the best test case
            h1 = diags[idx][1]
            lifetimes = np.abs(h1[:, 1] - h1[:, 0])
            # Degenerate case with points on a circle, some 1e-12 H1 lifetimes can appear
            epsilon = 1e-11
            assert np.sum(lifetimes > epsilon) == 1


def test_h1_only_cech_persistence_of_points_on_a_circle():
    pts = [points.sphere(n_samples=150, ambient_dim=2)]
    diags = CechPersistence(homology_dimensions=1, n_jobs=1).fit_transform(pts)[0]
    assert_almost_equal(np.array([[0.1, 1.]]), diags, decimal=1)

    # Same test but for the weighted version
    wpts = [np.append(pts[0], np.zeros((len(pts[0]), 1)), axis=1)]
    diags = WeightedCechPersistence(homology_dimensions=1, n_jobs=1).fit_transform(wpts)[0]
    for diag in diags:
        # H1 - can have some [1., 1.] - Maybe points on a circle is not the best test case
        if not math.isclose(diag[0], 1.):
            assert_almost_equal(np.array([0.1, 1.]), diag, decimal=1)

def test_set_output():
    try:
        import pandas

        NB_PC = 5
        n_samples_array = np.random.randint(100, high=150, size=NB_PC)
        pts = [points.sphere(n_samples=n_samples, ambient_dim=2) for n_samples in n_samples_array]

        cech = CechPersistence(homology_dimensions=[0, 2], n_jobs=-2)
        diags_pandas = cech.set_output(transform="pandas").fit_transform(pts)
        assert "H0" == diags_pandas.columns[0]
        assert "H2" == diags_pandas.columns[1]
        assert len(diags_pandas.index) == NB_PC

        # Same test but for the weighted version
        wpts = [np.append(pts[idx], np.zeros((len(pts[idx]), 1)), axis=1) for idx in range(NB_PC)]

        wcech = WeightedCechPersistence(homology_dimensions=[0, 2], n_jobs=-2)
        diags_pandas = wcech.set_output(transform="pandas").fit_transform(wpts)
        assert "H0" == diags_pandas.columns[0]
        assert "H2" == diags_pandas.columns[1]
        assert len(diags_pandas.index) == NB_PC

    except ImportError:
        print("Missing pandas, skipping set_output test")


def test_cech_persistence_constructor_exception():
    for CechConstructor in [CechPersistence, WeightedCechPersistence]:
        with pytest.raises(ValueError):
            cech = CechConstructor(homology_dimensions=[[0], [2]])
        with pytest.raises(ValueError):
            cech = CechConstructor(homology_dimensions=0, precision="some unvalid precision")
