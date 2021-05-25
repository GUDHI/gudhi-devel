""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2021 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.sklearn.cubical_persistence import CubicalPersistence
import numpy as np

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2021 Inria"
__license__ = "MIT"

def test_simple_constructor_from_top_cells():
    cp = CubicalPersistence(persistence_dim = 0)

    # The first "0" from sklearn.datasets.load_digits()
    bmp = np.array([[ 0.,  0.,  5., 13.,  9.,  1.,  0.,  0.],
                    [ 0.,  0., 13., 15., 10., 15.,  5.,  0.],
                    [ 0.,  3., 15.,  2.,  0., 11.,  8.,  0.],
                    [ 0.,  4., 12.,  0.,  0.,  8.,  8.,  0.],
                    [ 0.,  5.,  8.,  0.,  0.,  9.,  8.,  0.],
                    [ 0.,  4., 11.,  0.,  1., 12.,  7.,  0.],
                    [ 0.,  2., 14.,  5., 10., 12.,  0.,  0.],
                    [ 0.,  0.,  6., 13., 10.,  0.,  0.,  0.]])

    assert cp.fit_transform(bmp) == np.array([[0., 6.], [0., 8.]])

# from gudhi.representations import PersistenceImage
# PersistenceImage(bandwidth=50, weight=lambda x: x[1]**2, im_range=[0,256,0,256], resolution=[20, 20])
# PI.fit_transform([diag])