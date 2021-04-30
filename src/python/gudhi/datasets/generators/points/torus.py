# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Hind Montassif
#
# Copyright (C) 2021 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

import numpy as np
import math


def generate_random_points(num_points, dim):

    # Generate random angles of size num_points*dim
    alpha = 2*math.pi*np.random.rand(num_points*dim)

    # Based on angles, construct points of size num_points*dim on a circle and reshape the result in a num_points*2*dim array
    array_points = np.asarray([[np.cos(a), np.sin(a)] for a in alpha]).ravel().reshape(num_points, 2*dim)
    
    return array_points


def generate_grid_points(num_points, dim):
    
    num_points_grid = (int(num_points**(1./dim)))**dim
    
    alpha = 2*math.pi*np.random.rand(num_points_grid*dim)
    
    array_points = np.asarray([[np.cos(a), np.sin(a)] for a in alpha]).ravel().reshape(num_points_grid, 2*dim)
    
    return array_points

def generate_points(num_points, dim, sample='random'):
    if sample == 'random':
        print("Sample is random")
        npoints = num_points
    elif sample == 'grid':
        print("Sample is grid")
        npoints = (int(num_points**(1./dim)))**dim
    else:
        print("Sample type '{}' is not supported".format(sample))
        return
    
    # Generate random angles of size num_points*dim
    alpha = 2*math.pi*np.random.rand(npoints*dim)
    
    # Based on angles, construct points of size num_points*dim on a circle and reshape the result in a num_points*2*dim array
    array_points = np.asarray([[np.cos(a), np.sin(a)] for a in alpha]).ravel().reshape(npoints, 2*dim)
    
    return array_points
