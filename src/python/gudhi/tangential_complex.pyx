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

cdef extern from "Tangential_complex_interface.h" namespace "Gudhi":
    cdef cppclass Tangential_complex_interface "Gudhi::tangential_complex::Tangential_complex_interface":
        Tangential_complex_interface(int intrisic_dim, vector[vector[double]] points)
        # bool from_file is a workaround for cython to find the correct signature
        Tangential_complex_interface(int intrisic_dim, string off_file, bool from_file)
        void compute_tangential_complex() except +
        vector[double] get_point(unsigned vertex)
        unsigned number_of_vertices()
        unsigned number_of_simplices()
        unsigned number_of_inconsistent_simplices()
        unsigned number_of_inconsistent_stars()
        void create_simplex_tree(Simplex_tree_interface_full_featured* simplex_tree)
        void fix_inconsistencies_using_perturbation(double max_perturb, double time_limit)
        void set_max_squared_edge_length(double max_squared_edge_length)

# TangentialComplex python interface
cdef class TangentialComplex:
    """The class Tangential_complex represents a tangential complex. After the
    computation of the complex, an optional post-processing called perturbation
    can be run to attempt to remove inconsistencies.
    """

    cdef Tangential_complex_interface * thisptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, intrisic_dim, points=None, off_file=''):
        """TangentialComplex constructor.

        :param intrisic_dim: Intrinsic dimension of the manifold.
        :type intrisic_dim: integer

        :param points: A list of points in d-Dimension.
        :type points: list of list of double

        Or

        :param off_file: An OFF file style name.
        :type off_file: string
        """

    # The real cython constructor
    def __cinit__(self, intrisic_dim, points=None, off_file=''):
        if off_file:
            if os.path.isfile(off_file):
                self.thisptr = new Tangential_complex_interface(intrisic_dim, off_file.encode('utf-8'), True)
            else:
                print("file " + off_file + " not found.")
        else:
            if points is None:
                # Empty tangential construction
                points=[]
            self.thisptr = new Tangential_complex_interface(intrisic_dim, points)
                

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr

    def __is_defined(self):
        """Returns true if TangentialComplex pointer is not NULL.
         """
        return self.thisptr != NULL

    def compute_tangential_complex(self):
        """This function computes the tangential complex.

        Raises:
            ValueError: In debug mode, if the computed star dimension is too
                low. Try to set a bigger maximal edge length value with
                :meth:`set_max_squared_edge_length`
                if this happens.
        """
        self.thisptr.compute_tangential_complex()

    def get_point(self, vertex):
        """This function returns the point corresponding to a given vertex.

        :param vertex: The vertex.
        :type vertex: int.
        :returns:  The point.
        :rtype: list of float
        """
        cdef vector[double] point = self.thisptr.get_point(vertex)
        return point

    def num_vertices(self):
        """
        :returns:  The number of vertices.
        :rtype: unsigned
        """
        return self.thisptr.number_of_vertices()

    def num_simplices(self):
        """
        :returns:  Total number of simplices in stars (including duplicates that appear in several stars).
        :rtype: unsigned
        """
        return self.thisptr.number_of_simplices()

    def num_inconsistent_simplices(self):
        """
        :returns:  The number of inconsistent simplices.
        :rtype: unsigned
        """
        return self.thisptr.number_of_inconsistent_simplices()

    def num_inconsistent_stars(self):
        """
        :returns:  The number of stars containing at least one inconsistent simplex.
        :rtype: unsigned
        """
        return self.thisptr.number_of_inconsistent_stars()

    def create_simplex_tree(self):
        """Exports the complex into a simplex tree.

        :returns: A simplex tree created from the complex.
        :rtype: SimplexTree
        """
        stree = SimplexTree()
        cdef intptr_t stree_int_ptr=stree.thisptr
        self.thisptr.create_simplex_tree(<Simplex_tree_interface_full_featured*>stree_int_ptr)
        return stree

    def fix_inconsistencies_using_perturbation(self, max_perturb, time_limit=-1.0):
        """Attempts to fix inconsistencies by perturbing the point positions.

        :param max_perturb: Maximum length of the translations used by the
            perturbation.
        :type max_perturb: double
        :param time_limit: Time limit in seconds. If -1, no time limit is set.
        :type time_limit: double
        """
        self.thisptr.fix_inconsistencies_using_perturbation(max_perturb,
                                                            time_limit)

    def set_max_squared_edge_length(self, max_squared_edge_length):
        """Sets the maximal possible squared edge length for the edges in the
        triangulations.

        :param max_squared_edge_length: Maximal possible squared edge length.
        :type max_squared_edge_length: double

        If the maximal edge length value is too low
        :meth:`compute_tangential_complex`
        will throw an exception in debug mode.
        """
        self.thisptr.set_max_squared_edge_length(max_squared_edge_length)
