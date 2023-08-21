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
__license__ = "MIT (GPL v3 for sparsify_point_set)"

cdef extern from "Subsampling_interface.h" namespace "Gudhi::subsampling":
    vector[vector[double]] subsampling_n_farthest_points(bool, vector[vector[double]] points, size_t nb_points)
    vector[vector[double]] subsampling_n_farthest_points(bool, vector[vector[double]] points, size_t nb_points, size_t starting_point)
    vector[vector[double]] subsampling_n_farthest_points_from_file(bool, string off_file, size_t nb_points)
    vector[vector[double]] subsampling_n_farthest_points_from_file(bool, string off_file, size_t nb_points, size_t starting_point)
    vector[vector[double]] subsampling_n_random_points(vector[vector[double]] points, unsigned nb_points)
    vector[vector[double]] subsampling_n_random_points_from_file(string off_file, unsigned nb_points)
    vector[vector[double]] subsampling_sparsify_points(vector[vector[double]] points, double min_squared_dist)
    vector[vector[double]] subsampling_sparsify_points_from_file(string off_file, double min_squared_dist)

cdef extern from "Subsampling_interface.h":
    cdef int _GUDHI_SUBSAMPLING_USE_CGAL

GUDHI_SUBSAMPLING_USE_CGAL = _GUDHI_SUBSAMPLING_USE_CGAL

def choose_n_farthest_points(points=None, off_file='', nb_points=-<size_t>1, starting_point=None, fast=True):
    """Subsample by a greedy strategy of iteratively adding the farthest point
    from the current chosen point set to the subsampling.
    The iteration starts with the landmark `starting point`.

    :param points: The input point set.
    :type points: Iterable[Iterable[float]]

    Or

    :param off_file: An OFF file style name.
    :type off_file: string

    And in both cases

    :param nb_points: Number of points of the subsample (the subsample may be
        smaller if there are fewer than nb_points distinct input points). Default: all of them.
    :type nb_points: int
    :param starting_point: The iteration starts with the landmark `starting
        point`, which is the index of the point to start with. If not set, this
        index is chosen randomly.
    :type starting_point: int
    :param fast: If True (default), use an implementation that is efficient when the doubling dimension
        and spread are small, but slow otherwise. If False, use the standard quadratic algorithm.
    :type fast: bool
    :returns:  The subsample point set, in the order they were selected by the greedy strategy.
    :rtype: List[List[float]]
    """
    if off_file:
        if os.path.isfile(off_file):
            if starting_point is None:
                return subsampling_n_farthest_points_from_file(fast, off_file.encode('utf-8'),
                                                               nb_points)
            else:
                return subsampling_n_farthest_points_from_file(fast, off_file.encode('utf-8'),
                                                               nb_points, starting_point)
        else:
            print("file " + off_file + " not found.")
    else:
        if points is None:
            # Empty points
            points=[]
        if starting_point is None:
            return subsampling_n_farthest_points(fast, points, nb_points)
        else:
            return subsampling_n_farthest_points(fast, points, nb_points, starting_point)

def pick_n_random_points(points=None, off_file='', nb_points=0):
    """Subsample a point set by picking random vertices.

    :param points: The input point set.
    :type points: Iterable[Iterable[float]]

    Or

    :param off_file: An OFF file style name.
    :type off_file: string

    And in both cases

    :param nb_points: Number of points of the subsample.
    :type nb_points: int
    :returns:  The subsample point set.
    :rtype: List[List[float]]
    """
    if off_file:
        if os.path.isfile(off_file):
            return subsampling_n_random_points_from_file(off_file.encode('utf-8'),
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
    between any two points is greater than min_squared_dist.

    :param points: The input point set.
    :type points: Iterable[Iterable[float]]

    Or

    :param off_file: An OFF file style name.
    :type off_file: string

    And in both cases

    :param min_squared_dist: Minimum squared distance separating the output \
    points.
    :type min_squared_dist: float
    :returns:  The subsample point set.
    :rtype: List[List[float]]
    """
    if not GUDHI_SUBSAMPLING_USE_CGAL:
        raise NotImplementedError("subsampling sparsify_point_set is only available with CGAL >= 4.11 and Eigen3")
    
    if off_file:
        if os.path.isfile(off_file):
            return subsampling_sparsify_points_from_file(off_file.encode('utf-8'),
                                                         min_squared_dist)
        else:
            print("file " + off_file + " not found.")
    else:
        if points is None:
            # Empty points
            points=[]
        return subsampling_sparsify_points(points, min_squared_dist)
