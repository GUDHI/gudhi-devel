""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2021 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.sklearn.cubical_persistence import CubicalPersistence
import numpy as np
from sklearn import datasets

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2021 Inria"
__license__ = "MIT"

def test_simple_constructor_from_top_cells():
    cells = datasets.load_digits().images[0]
    cp = CubicalPersistence(persistence_dim = 0)
    np.testing.assert_array_equal(cp._CubicalPersistence__transform(cells),
                                  np.array([[0., 6.], [0., 8.]]))

def test_simple_constructor_from_top_cells_list():
    digits = datasets.load_digits().images[:10]
    cp = CubicalPersistence(persistence_dim = 0, n_jobs=-2)

    diags = cp.fit_transform(digits)
    assert len(diags) == 10
    np.testing.assert_array_equal(diags[0],
                                  np.array([[0., 6.], [0., 8.]]))

# from gudhi.representations import PersistenceImage
# pi = PersistenceImage(bandwidth=50, weight=lambda x: x[1]**2, im_range=[0,256,0,256], resolution=[20, 20])
# pi.fit_transform(diags)