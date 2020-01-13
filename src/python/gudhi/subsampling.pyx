# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
import os

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "GPL v3"

cdef extern from "Subsampling_interface.h" namespace "Gudhi::subsampling":
    vector[vector[double]] subsampling_n_farthest_points(vector[vector[double]] points, unsigned nb_points)
    vector[vector[double]] subsampling_n_farthest_points(vector[vector[double]] points, unsigned nb_points, unsigned starting_point)
    vector[vector[double]] subsampling_n_farthest_points_from_file(string off_file, unsigned nb_points)
    vector[vector[double]] subsampling_n_farthest_points_from_file(string off_file, unsigned nb_points, unsigned starting_point)
    vector[vector[double]] subsampling_n_random_points(vector[vector[double]] points, unsigned nb_points)
    vector[vector[double]] subsampling_n_random_points_from_file(string off_file, unsigned nb_points)
    vector[vector[double]] subsampling_sparsify_points(vector[vector[double]] points, double min_squared_dist)
    vector[vector[double]] subsampling_sparsify_points_from_file(string off_file, double min_squared_dist)

def choose_n_farthest_points(points=None, off_file='', nb_points=0, starting_point = ''):
    """Subsample by a greedy strategy of iteratively adding the farthest point
    from the current chosen point set to the subsampling.
    The iteration starts with the landmark `starting point`.

    :param points: The input point set.
    :type points: vector[vector[double]].

    Or

    :param off_file: An OFF file style name.
    :type off_file: string

    :param nb_points: Number of points of the subsample.
    :type nb_points: unsigned.
    :param starting_point: The iteration starts with the landmark `starting \
    point`,which is the index of the point to start with. If not set, this \
    index is chosen randomly.
    :type starting_point: unsigned.
    :returns:  The subsample point set.
    :rtype: vector[vector[double]]
    """
    if off_file:
        if os.path.isfile(off_file):
            if starting_point == '':
                return subsampling_n_farthest_points_from_file(str.encode(off_file),
                                                               nb_points)
            else:
                return subsampling_n_farthest_points_from_file(str.encode(off_file),
                                                               nb_points,
                                                               starting_point)
        else:
            print("file " + off_file + " not found.")
    else:
        if points is None:
            # Empty points
            points=[]
        if starting_point == '':
            return subsampling_n_farthest_points(points, nb_points)
        else:
            return subsampling_n_farthest_points(points, nb_points,
                                                 starting_point)

def pick_n_random_points(points=None, off_file='', nb_points=0):
    """Subsample a point set by picking random vertices.

    :param points: The input point set.
    :type points: vector[vector[double]].

    Or

    :param off_file: An OFF file style name.
    :type off_file: string

    :param nb_points: Number of points of the subsample.
    :type nb_points: unsigned.
    :returns:  The subsample point set.
    :rtype: vector[vector[double]]
    """
    if off_file:
        if os.path.isfile(off_file):
            return subsampling_n_random_points_from_file(str.encode(off_file),
                nb_points)
        else:
            print("file " + off_file + " not found.")
    else:
        if points is None:
            # Empty points
            points=[]
        return subsampling_n_random_points(points, nb_points)

def sparsify_point_set(points=None, off_file='', min_squared_dist=0.0):
    """Outputs a subset of the input points so that the squared distance
    between any two points is greater than or equal to min_squared_dist.

    :param points: The input point set.
    :type points: vector[vector[double]].

    Or

    :param off_file: An OFF file style name.
    :type off_file: string

    :param min_squared_dist: Minimum squared distance separating the output \
    points.
    :type min_squared_dist: float.
    :returns:  The subsample point set.
    :rtype: vector[vector[double]]
    """
    if off_file:
        if os.path.isfile(off_file):
            return subsampling_sparsify_points_from_file(str.encode(off_file),
                                                         min_squared_dist)
        else:
            print("file " + off_file + " not found.")
    else:
        if points is None:
            # Empty points
            points=[]
        return subsampling_sparsify_points(points, min_squared_dist)
