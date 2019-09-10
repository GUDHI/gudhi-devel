from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libc.stdint cimport intptr_t

from gudhi.simplex_tree cimport *
from gudhi.simplex_tree import SimplexTree

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "GPL v3"

cdef extern from "Euclidean_strong_witness_complex_interface.h" namespace "Gudhi":
    cdef cppclass Euclidean_strong_witness_complex_interface "Gudhi::witness_complex::Euclidean_strong_witness_complex_interface":
        Euclidean_strong_witness_complex_interface(vector[vector[double]] landmarks, vector[vector[double]] witnesses)
        void create_simplex_tree(Simplex_tree_interface_full_featured* simplex_tree, double max_alpha_square)
        void create_simplex_tree(Simplex_tree_interface_full_featured* simplex_tree, double max_alpha_square,
            unsigned limit_dimension)
        vector[double] get_point(unsigned vertex)

# EuclideanStrongWitnessComplex python interface
cdef class EuclideanStrongWitnessComplex:
    """Constructs strong witness complex for given sets of witnesses and
    landmarks in Euclidean space.
    """

    cdef Euclidean_strong_witness_complex_interface * thisptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, landmarks=None, witnesses=None):
        """WitnessComplex constructor.

        :param landmarks: A list of landmarks (in the point cloud).
        :type landmarks: list of list of double

        :param witnesses: The point cloud.
        :type witnesses: list of list of double
        """

    # The real cython constructor
    def __cinit__(self, landmarks=None, witnesses=None):
        if landmarks is not None and witnesses is not None:
            self.thisptr = new Euclidean_strong_witness_complex_interface(landmarks, witnesses)

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr

    def __is_defined(self):
        """Returns true if WitnessComplex pointer is not NULL.
         """
        return self.thisptr != NULL

    def create_simplex_tree(self, max_alpha_square, limit_dimension = -1):
        """
        :param max_alpha_square: The maximum alpha square threshold the
            simplices shall not exceed. Default is set to infinity.
        :type max_alpha_square: float
        :returns: A simplex tree created from the Delaunay Triangulation.
        :rtype: SimplexTree
        """
        stree = SimplexTree()
        cdef intptr_t stree_int_ptr=stree.thisptr
        cdef Simplex_tree_interface_full_featured* stree_ptr = <Simplex_tree_interface_full_featured*>stree_int_ptr
        if limit_dimension is not -1:
            self.thisptr.create_simplex_tree(stree_ptr, max_alpha_square, limit_dimension)
        else:
            self.thisptr.create_simplex_tree(stree_ptr, max_alpha_square)
        return stree

    def get_point(self, vertex):
        """This function returns the point corresponding to a given vertex.

        :param vertex: The vertex.
        :type vertex: int.
        :returns:  The point.
        :rtype: list of float
        """
        cdef vector[double] point = self.thisptr.get_point(vertex)
        return point

