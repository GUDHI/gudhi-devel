""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2021 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.sklearn.cubical_persistence import CubicalPersistence
import gudhi
import numpy as np
from sklearn import datasets

CUBICAL_PERSISTENCE_H0_IMG0 = np.array([[0.0, 6.0], [0.0, 8.0], [0.0, np.inf]])


def test_simple_constructor_from_top_cells_list():
    digits = datasets.load_digits().images[:10]
    cp = CubicalPersistence(homology_dimensions=0, n_jobs=-2)

    diags = cp.fit_transform(digits)
    assert len(diags) == 10
    np.testing.assert_array_equal(diags[0], CUBICAL_PERSISTENCE_H0_IMG0)

    cp = CubicalPersistence(homology_dimensions=[0, 1], n_jobs=-1)
    diagsH0H1 = cp.fit_transform(digits)
    assert len(diagsH0H1) == 10
    for idx in range(10):
        np.testing.assert_array_equal(diags[idx], diagsH0H1[idx][0])


def test_simple_constructor_from_top_cells_list():
    digits = datasets.load_digits().images[:10]
    cp = CubicalPersistence(homology_dimensions=0, input_type='vertices', n_jobs=-2)
    diags = cp.fit_transform(digits)
    assert len(diags) == 10
    np.testing.assert_array_equal(diags[2], [[8, 13], [0, 15], [0, np.inf]])


def test_1d():
    a = np.array([2, 4, 3, 5])
    r = np.array([[3, 4], [2, np.inf]])
    ri = CubicalPersistence(0).fit_transform([a])[0]
    rf = CubicalPersistence(0).fit_transform([a.astype("float32")[::-1]])[0]
    rd = CubicalPersistence([0]).fit_transform([a.astype("float64")])[0][0]
    assert ri.dtype == np.dtype("float64")
    assert rf.dtype == np.dtype("float32")
    assert rd.dtype == np.dtype("float64")
    np.testing.assert_array_equal(r, ri)
    np.testing.assert_array_equal(r, rf)
    np.testing.assert_array_equal(r, rd)

def cmp(a):
    dim_list = list(range(len(a.shape)))
    d1 = CubicalPersistence(homology_dimensions=dim_list, input_type='top_dimensional_cells').fit_transform([a])[0]
    p = gudhi.CubicalComplex(top_dimensional_cells=a).persistence()
    # Ideally we would compare multisets, but length + set should be good enough.
    assert sum(len(d) for d in d1) == len(p)
    d1 = [{(a, b) for [a, b] in d} for d in d1]
    d2 = [set() for _ in dim_list]
    for (dim, i) in p:
        d2[dim].add(i)
    assert d1 == d2

def test_compare_top_cells():
    cmp(np.random.rand(20,15,10))
    cmp(np.random.rand(20,10))
    cmp(np.random.rand(20,2))
    cmp(np.random.rand(2,20))
    cmp(np.random.rand(20,1))
    cmp(np.random.rand(1,20))
    cmp(np.random.rand(20))
    a = np.random.rand(20,10)
    cmp(np.array(a, order='F'))
    cmp(a[0:18])
    cmp(a[::2,::-1])
