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

def _generate_random_points_on_torus(n_samples, dim):

    # Generate random angles of size n_samples*dim
    alpha = 2*np.pi*np.random.rand(n_samples*dim)

    # Based on angles, construct points of size n_samples*dim on a circle and reshape the result in a n_samples*2*dim array
    array_points = np.column_stack([np.cos(alpha), np.sin(alpha)]).reshape(-1, 2*dim)
    
    return array_points

def _generate_grid_points_on_torus(n_samples, dim):
    
    # Generate points on a dim-torus as a grid
    n_samples_grid = int(n_samples**(1./dim))
    alpha = np.linspace(0, 2*np.pi, n_samples_grid, endpoint=False)
    
    array_points_inter = np.column_stack([np.cos(alpha), np.sin(alpha)])
    array_points = np.array(list(itertools.product(array_points_inter, repeat=dim))).reshape(-1, 2*dim)
    
    return array_points

def torus(n_samples, dim, sample='random'):
    """ 
    Generate points on a dim-torus in R^2dim either randomly or on a grid
    
    :param n_samples: The number of points to be generated.
    :param dim: The dimension of the torus on which points would be generated in R^2*dim.
    :param sample: The sample type of the generated points. Can be 'random' or 'grid'.
    :returns: numpy array containing the generated points on a torus.
        The shape of returned numpy array is :
        if sample is 'random' : (n_samples, 2*dim).
        if sample is 'grid' : ([n_samples**(1./dim)]**dim, 2*dim).
    """
    if sample == 'random':
        # Generate points randomly
        print("Sample is random")
        return _generate_random_points_on_torus(n_samples, dim)
    elif sample == 'grid':
        # Generate points on a grid
        print("Sample is grid")
        return _generate_grid_points_on_torus(n_samples, dim)
    else:
        raise ValueError("Sample type '{}' is not supported".format(sample))
        return