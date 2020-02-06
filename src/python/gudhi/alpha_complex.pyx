# This file is part of the Gudhi Library - https://gudhi.inria.fr/ -
# which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full
# license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
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
import os

from gudhi.simplex_tree cimport *
from gudhi.simplex_tree import SimplexTree

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "GPL v3"

cdef extern from "Alpha_complex_interface.h" namespace "Gudhi":
    cdef cppclass Alpha_complex_interface "Gudhi::alpha_complex::Alpha_complex_interface":
        Alpha_complex_interface(vector[vector[double]] points) except +
        # bool from_file is a workaround for cython to find the correct signature
        Alpha_complex_interface(string off_file, bool from_file) except +
        vector[double] get_point(int vertex) except +
        void create_simplex_tree(Simplex_tree_interface_full_featured* simplex_tree, double max_alpha_square) except +

# AlphaComplex python interface
cdef class AlphaComplex:
    """AlphaComplex is a simplicial complex constructed from the finite cells
    of a Delaunay Triangulation.

    The filtration value of each simplex is computed as the square of the
    circumradius of the simplex if the circumsphere is empty (the simplex is
    then said to be Gabriel), and as the minimum of the filtration values of
    the codimension 1 cofaces that make it not Gabriel otherwise.

    All simplices that have a filtration value strictly greater than a given
    alpha squared value are not inserted into the complex.

    .. note::

        When Alpha_complex is constructed with an infinite value of alpha, the
        complex is a Delaunay complex.

    """

    cdef Alpha_complex_interface * thisptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, points=None, off_file=''):
        """AlphaComplex constructor.

        :param points: A list of points in d-Dimension.
        :type points: list of list of double

        Or

        :param off_file: An OFF file style name.
        :type off_file: string
        """

    # The real cython constructor
    def __cinit__(self, points = None, off_file = ''):
        if off_file:
            if os.path.isfile(off_file):
                self.thisptr = new Alpha_complex_interface(
                    off_file.encode('utf-8'), True)
            else:
                print("file " + off_file + " not found.")
        else:
            if points is None:
                # Empty Alpha construction
                points=[]
            self.thisptr = new Alpha_complex_interface(points)
                

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr

    def __is_defined(self):
        """Returns true if AlphaComplex pointer is not NULL.
         """
        return self.thisptr != NULL

    def get_point(self, vertex):
        """This function returns the point corresponding to a given vertex.

        :param vertex: The vertex.
        :type vertex: int
        :rtype: list of float
        :returns: the point.
        """
        return self.thisptr.get_point(vertex)

    def create_simplex_tree(self, max_alpha_square = float('inf')):
        """
        :param max_alpha_square: The maximum alpha square threshold the
            simplices shall not exceed. Default is set to infinity, and
            there is very little point using anything else since it does
            not save time.
        :type max_alpha_square: float
        :returns: A simplex tree created from the Delaunay Triangulation.
        :rtype: SimplexTree
        """
        stree = SimplexTree()
        cdef intptr_t stree_int_ptr=stree.thisptr
        self.thisptr.create_simplex_tree(<Simplex_tree_interface_full_featured*>stree_int_ptr, max_alpha_square)
        return stree
