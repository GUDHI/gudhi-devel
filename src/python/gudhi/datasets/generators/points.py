# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Hind Montassif
#
# Copyright (C) 2021 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

import numpy as np

from ._points import ctorus
from ._points import sphere

def _generate_random_points_on_torus(n_samples, dim):

    # Generate random angles of size n_samples*dim
    alpha = 2*np.pi*np.random.rand(n_samples*dim)

    # Based on angles, construct points of size n_samples*dim on a circle and reshape the result in a
    # n_samples*2*dim array
    array_points = np.column_stack([np.cos(alpha), np.sin(alpha)]).reshape(-1, 2*dim)

    return array_points

def _generate_grid_points_on_torus(n_samples, dim):

    # Generate points on a dim-torus as a grid
    n_samples_grid = int((n_samples+.5)**(1./dim)) # add .5 to avoid rounding down with numerical approximations
    alpha = np.linspace(0, 2*np.pi, n_samples_grid, endpoint=False)

    array_points = np.column_stack([np.cos(alpha), np.sin(alpha)])
    array_points_idx = np.empty([n_samples_grid]*dim + [dim], dtype=int)
    for i, x in enumerate(np.ix_(*([np.arange(n_samples_grid)]*dim))):
        array_points_idx[...,i] = x
    return array_points[array_points_idx].reshape(-1, 2*dim)

def torus(n_samples, dim, sample='random'):
    """
    Generate points on a flat dim-torus in R^2dim either randomly or on a grid

    :param n_samples: The number of points to be generated.
    :param dim: The dimension of the torus on which points would be generated in R^2*dim.
    :param sample: The sample type of the generated points. Can be 'random' or 'grid'.
    :returns: numpy array containing the generated points on a torus.

    The shape of returned numpy array is:

    If sample is 'random': (n_samples, 2*dim).

    If sample is 'grid': (⌊n_samples**(1./dim)⌋**dim, 2*dim), where shape[0] is rounded down to the closest perfect
    'dim'th power.
    """
    if sample == 'random':
        # Generate points randomly
        return _generate_random_points_on_torus(n_samples, dim)
    elif sample == 'grid':
        # Generate points on a grid
        return _generate_grid_points_on_torus(n_samples, dim)
    else:
        raise ValueError(f"Sample type '{sample}' is not supported")
