# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Hind Montassif
#
# Copyright (C) 2021 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

import numpy as np
import itertools

def _generate_random_points(n_samples, dim):

    # Generate random angles of size n_samples*dim
    alpha = 2*np.pi*np.random.rand(n_samples*dim)

    # Based on angles, construct points of size n_samples*dim on a circle and reshape the result in a n_samples*2*dim array
    array_points = np.column_stack([np.cos(alpha), np.sin(alpha)]).reshape(-1, 2*dim)
    
    return array_points

def _generate_grid_points(n_samples, dim):
    
    n_samples_grid = int(n_samples**(1./dim))
    alpha = np.linspace(0, 2*np.pi, n_samples_grid, endpoint=False)
    
    array_points_inter = np.column_stack([np.cos(alpha), np.sin(alpha)])
    array_points = np.array(list(itertools.product(array_points_inter, repeat=dim))).reshape(-1, 2*dim)
    
    return array_points

def torus(n_samples, dim, sample='random'):
    if sample == 'random':
        print("Sample is random")
        return _generate_random_points(n_samples, dim)
    elif sample == 'grid':
        print("Sample is grid")
        return _generate_grid_points(n_samples, dim)
    else:
        raise Exception("Sample type '{}' is not supported".format(sample))
        return
