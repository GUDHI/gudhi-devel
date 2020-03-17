from gudhi.barycenter import lagrangian_barycenter
import numpy as np

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Theo Lacombe

    Copyright (C) 2019 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__author__ = "Theo Lacombe"
__copyright__ = "Copyright (C) 2019 Inria"
__license__ = "MIT"


def test_lagrangian_barycenter():
     
    dg1 = np.array([[0.2, 0.5]])
    dg2 = np.array([[0.2, 0.7]])
    dg3 = np.array([[0.3, 0.6], [0.7, 0.8], [0.2, 0.3]])
    dg4 = np.array([])
    dg5 = np.array([])
    dg6 = np.array([])
    res = np.array([[0.27916667, 0.55416667], [0.7375,  0.7625], [0.2375, 0.2625]])

    dg7 = np.array([[0.1, 0.15], [0.1, 0.7], [0.2, 0.22], [0.55, 0.84], [0.11, 0.91], [0.61, 0.75], [0.33, 0.46], [0.12, 0.41], [0.32, 0.48]])
    dg8 = np.array([[0., 4.], [4, 8]])
   
    # error crit.
    eps = 1e-7


    assert np.linalg.norm(lagrangian_barycenter(pdiagset=[dg1, dg2, dg3, dg4],init=3, verbose=False) - res) < eps 
    assert np.array_equal(lagrangian_barycenter(pdiagset=[dg4, dg5, dg6], verbose=False), np.empty(shape=(0,2)))
    assert np.linalg.norm(lagrangian_barycenter(pdiagset=[dg7], verbose=False) - dg7) < eps
    Y, log = lagrangian_barycenter(pdiagset=[dg4, dg8], verbose=True)
    assert np.linalg.norm(Y - np.array([[1,3], [5, 7]])) < eps
    assert np.abs(log["energy"] - 2) < eps
    assert np.array_equal(log["groupings"][0] , np.array([[0, -1], [1, -1]]))
    assert np.array_equal(log["groupings"][1] , np.array([[0, 0], [1, 1]]))
    assert np.linalg.norm(lagrangian_barycenter(pdiagset=[dg8, dg4], init=np.array([[0.2, 0.6], [0.5, 0.7]]), verbose=False) - np.array([[1, 3], [5, 7]])) < eps
    assert lagrangian_barycenter(pdiagset = []) is None

