""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2023 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi.datasets.generators import points
from gudhi.sklearn.rips_persistence import RipsPersistence
import random

def test_rips_persistence_of_points_on_a_circle():
    nb_points = random.randint(100,150)
    rips = RipsPersistence(homology_dimensions=[0, 1], n_jobs=-2)
    diags = rips.fit_transform([points.sphere(n_samples = nb_points, ambient_dim = 2)])
    # list of one array as an input means list of one array as an output
    assert len(diags) == 1
    # H0 + H1
    assert len(diags[0]) == 2
    # H0
    assert len(diags[0][0]) == nb_points
    # H1
    assert len(diags[0][1]) == 1

def test_h1_only_rips_persistence_of_points_on_a_circle():
    rips = RipsPersistence(homology_dimensions=1, n_jobs=-2)
    diags = rips.fit_transform([points.sphere(n_samples = 100, ambient_dim = 2)])[0]
    assert len(diags) == 1
    assert 0. < diags[0][0] < 0.5
    assert 1. < diags[0][1] < 2.

