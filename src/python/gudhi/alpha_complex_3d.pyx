# This file is part of the Gudhi Library - https://gudhi.inria.fr/ -
# which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full
# license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2021 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from __future__ import print_function
from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp cimport bool
from libc.stdint cimport intptr_t
import errno
import os
import warnings

from gudhi.simplex_tree cimport *
from gudhi.simplex_tree import SimplexTree
from gudhi import read_points_from_off_file

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2021 Inria"
__license__ = "GPL v3"

cdef extern from "Alpha_complex_interface_3d.h" namespace "Gudhi":
    cdef cppclass Alpha_complex_interface_3d "Gudhi::alpha_complex::Alpha_complex_interface_3d":
        Alpha_complex_interface_3d(vector[vector[double]] points, vector[double] weights, bool fast_version, bool exact_version) nogil except +
        vector[double] get_point(int vertex) nogil except +
        void create_simplex_tree(Simplex_tree_interface_full_featured* simplex_tree, double max_alpha_square) nogil except +

# AlphaComplex3D python interface
cdef class AlphaComplex3D:
    """AlphaComplex3D is a simplicial complex constructed from the finite cells of a Delaunay Triangulation.

    The filtration value of each simplex is computed as the square of the circumradius of the simplex if the
    circumsphere is empty (the simplex is then said to be Gabriel), and as the minimum of the filtration values of the
    codimension 1 cofaces that make it not Gabriel otherwise.

    All simplices that have a filtration value strictly greater than a given alpha squared value are not inserted into
    the complex.

    .. note::

        When AlphaComplex3D is constructed with an infinite value of alpha, the complex is a Delaunay complex.

    .. warning::

        Contrary to the dD version, with the 3d version, the vertices in the output simplex tree are not guaranteed to
        match the order of the input points. One can use :func:`~gudhi.AlphaComplex3D.get_point` to get the initial
        point back.
    """

    cdef Alpha_complex_interface_3d * this_ptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, points=[], weights=[], precision='safe'):
        """AlphaComplex3D constructor.

        :param points: A list of points in 3d.
        :type points: Iterable[Iterable[float]]

        :param weights: A list of weights. If set, the number of weights must correspond to the number of points.
        :type weights: Iterable[float]

        :param precision: Alpha complex precision can be 'fast', 'safe' or 'exact'. Default is 'safe'.
        :type precision: string

        :raises ValueError: If the given points are not in 3d.
        :raises ValueError: In case of inconsistency between the number of points and weights.
        """

    # The real cython constructor
    def __cinit__(self, points = [], weights=[], precision = 'safe'):
        assert precision in ['fast', 'safe', 'exact'], "Alpha complex precision can only be 'fast', 'safe' or 'exact'"
        cdef bool fast = precision == 'fast'
        cdef bool exact = precision == 'exact'

        if len(points) > 0:
            if len(points[0]) != 3:
                raise ValueError("AlphaComplex3D only accepts 3d points as an input")

        # weights are set but is inconsistent with the number of points
        if len(weights) != 0 and len(weights) != len(points):
            raise ValueError("Inconsistency between the number of points and weights")

        # need to copy the points to use them without the gil
        cdef vector[vector[double]] pts
        cdef vector[double] wgts
        pts = points
        wgts = weights
        with nogil:
            self.this_ptr = new Alpha_complex_interface_3d(pts, wgts, fast, exact)

    def __dealloc__(self):
        if self.this_ptr != NULL:
            del self.this_ptr

    def __is_defined(self):
        """Returns true if AlphaComplex3D pointer is not NULL.
         """
        return self.this_ptr != NULL

    def get_point(self, vertex):
        """This function returns the point corresponding to a given vertex from the :class:`~gudhi.SimplexTree`.

        :param vertex: The vertex.
        :type vertex: int
        :rtype: list of float
        :returns: the point.
        """
        return self.this_ptr.get_point(vertex)

    def create_simplex_tree(self, max_alpha_square = float('inf')):
        """
        :param max_alpha_square: The maximum alpha square threshold the simplices shall not exceed. Default is set to
            infinity, and there is very little point using anything else since it does not save time.
        :type max_alpha_square: float
        :returns: A simplex tree created from the Delaunay Triangulation.
        :rtype: SimplexTree

        :raises ValueError: If the points given at construction time are on a plane.
        """
        stree = SimplexTree()
        cdef double mas = max_alpha_square
        cdef intptr_t stree_int_ptr=stree.thisptr
        with nogil:
            self.this_ptr.create_simplex_tree(<Simplex_tree_interface_full_featured*>stree_int_ptr,
                                              mas)
        return stree
