# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Hind Montassif
#
# Copyright (C) 2021 Inria
#
# Modification(s):
#   - 2026/01 Vincent Rouvreau: Add gudhi.random.GudhiBitGenerator support to set the seed
#   - YYYY/MM Author: Description of the modification

__license__ = "MIT"


import numpy as np
import os
from pathlib import Path

from numpy.random import Generator
from gudhi.random import GudhiBitGenerator

# Specific for Windows - To be done before import _ctorus and _sphere
if os.name == 'nt':
    # Add random.dll that is located in "../.."
    dll_dir = Path(__file__).parent / "../.."
    os.add_dll_directory(str(dll_dir.resolve()))
    
from ._points_ext import _ctorus, _sphere

def _generate_random_points_on_torus(n_samples: int, dim: int, rng: GudhiBitGenerator = GudhiBitGenerator()):

    # Generate random angles of size n_samples*dim
    alpha = 2 * np.pi * Generator(rng).random(n_samples * dim)

    # Based on angles, construct points of size n_samples*dim on a circle and reshape the result in a
    # n_samples*2*dim array
    array_points = np.column_stack([np.cos(alpha), np.sin(alpha)]).reshape(-1, 2 * dim)

    return array_points


def _generate_grid_points_on_torus(n_samples: int, dim: int):

    # Generate points on a dim-torus as a grid
    n_samples_grid = int(
        (n_samples + 0.5) ** (1.0 / dim)
    )  # add .5 to avoid rounding down with numerical approximations
    alpha = np.linspace(0, 2 * np.pi, n_samples_grid, endpoint=False)

    array_points = np.column_stack([np.cos(alpha), np.sin(alpha)])
    array_points_idx = np.empty([n_samples_grid] * dim + [dim], dtype=int)
    for i, x in enumerate(np.ix_(*([np.arange(n_samples_grid)] * dim))):
        array_points_idx[..., i] = x
    return array_points[array_points_idx].reshape(-1, 2 * dim)


def torus(n_samples: int, dim: int, sample: str = "random", rng: GudhiBitGenerator = GudhiBitGenerator()):
    """
    Generate points on a flat dim-torus in R^2dim either randomly or on a grid

    :param n_samples: The number of points to be generated.
    :param dim: The dimension of the torus on which points would be generated in R^2*dim.
    :param sample: The sample type of the generated points. Can be 'random' or 'grid'.
    :param rng: Default is `GudhiBitGenerator()`, but it can be set with a seed with `gudhi.random.GudhiBitGenerator(seed)`.
        It has no effect when sample is 'grid'.
    :type rng: gudhi.random.GudhiBitGenerator
    :returns: numpy array containing the generated points on a torus.

    The shape of returned numpy array is:

    If sample is 'random': (n_samples, 2*dim).

    If sample is 'grid': (⌊n_samples**(1./dim)⌋**dim, 2*dim), where shape[0] is rounded down to the closest perfect
    'dim'th power.
    """
    if sample == "random":
        # Generate points randomly
        return _generate_random_points_on_torus(n_samples, dim)
    elif sample == "grid":
        # Generate points on a grid
        return _generate_grid_points_on_torus(n_samples, dim)
    else:
        raise ValueError(f"Sample type '{sample}' is not supported")

def sphere(n_samples: int, ambient_dim: int, radius: float = 1., sample: str = "random",
           rng: GudhiBitGenerator = GudhiBitGenerator()):
    """
    Generate random i.i.d. points uniformly on a (d-1)-sphere in R^d
    
    :param n_samples: The number of points to be generated.
    :type n_samples: integer
    :param ambient_dim: The ambient dimension d.
    :type ambient_dim: integer
    :param radius: The radius. Default value is `1.`.
    :type radius: float
    :param sample: The sample type. Default and only available value is `"random"`.
    :type sample: string
    :param rng: Default is `GudhiBitGenerator()`, but it can be set with a seed with `gudhi.random.GudhiBitGenerator(seed)`.
    :type rng: gudhi.random.GudhiBitGenerator
    :returns: the generated points on a sphere.
    """
    return _sphere(n_samples, ambient_dim, radius, sample, rng.rng)

def ctorus(n_samples: int, dim: int, sample: str = "random", rng: GudhiBitGenerator = GudhiBitGenerator()):
    """
    Generate random i.i.d. points on a d-torus in R^2d or as a grid
    
    :param n_samples: The number of points to be generated.
    :type n_samples: integer
    :param dim: The dimension of the torus on which points would be generated in R^2*dim.
    :type dim: integer
    :param sample: The sample type. Available values are: `"random"` and `"grid"`. Default value is `"random"`.
    :type sample: string
    :param rng: Default is `GudhiBitGenerator()`, but it can be set with a seed with `gudhi.random.GudhiBitGenerator(seed)`.
        It has no effect when sample is 'grid'.
    :type rng: gudhi.random.GudhiBitGenerator
    :returns: the generated points on a torus.
    
    The shape of returned numpy array is:
        
    If sample is 'random': (n_samples, 2*dim).
    
    If sample is 'grid': (⌊n_samples**(1./dim)⌋**dim, 2*dim), where shape[0] is rounded down to the closest perfect 'dim'th power.
    """
    return _ctorus(n_samples, dim, sample, rng.rng)
