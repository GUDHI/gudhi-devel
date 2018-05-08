from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
import os

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Vincent Rouvreau

   Copyright (C) 2016 Inria

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

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
    point`,which is the index of the poit to start with. If not set, this \
    index is choosen randomly.
    :type starting_point: unsigned.
    :returns:  The subsample point set.
    :rtype: vector[vector[double]]
    """
    if off_file is not '':
        if os.path.isfile(off_file):
            if starting_point is '':
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
        if starting_point is '':
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
    if off_file is not '':
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
    """Subsample a point set by picking random vertices.

    :param points: The input point set.
    :type points: vector[vector[double]].

    Or

    :param off_file: An OFF file style name.
    :type off_file: string

    :param min_squared_dist: Number of points of the subsample.
    :type min_squared_dist: unsigned.
    :returns:  The subsample point set.
    :rtype: vector[vector[double]]
    """
    if off_file is not '':
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
