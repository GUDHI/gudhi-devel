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

def generate_random_points(num_points, dim):

    # Generate random angles of size num_points*dim
    alpha = 2*np.pi*np.random.rand(num_points*dim)

    # Based on angles, construct points of size num_points*dim on a circle and reshape the result in a num_points*2*dim array    
    array_points = np.column_stack([np.cos(alpha), np.sin(alpha)]).reshape(-1, 2*dim)
    
    return array_points

def generate_grid_points(num_points, dim):
    
    num_points_grid = int(num_points**(1./dim))
    alpha = np.linspace(0, 2*np.pi, num_points_grid, endpoint=False)
    
    array_points_inter = np.column_stack([np.cos(alpha), np.sin(alpha)])
    array_points = np.array(list(itertools.product(array_points_inter, repeat=dim))).reshape(-1, 2*dim)
    
    return array_points

def generate_points(num_points, dim, sample='random'):
    if sample == 'random':
        print("Sample is random")
        generate_random_points(num_points, dim)
    elif sample == 'grid':
        print("Sample is grid")
        generate_grid_points(num_points, dim)
    else:
        print("Sample type '{}' is not supported".format(sample))
        return
