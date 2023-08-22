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

from gudhi.simplex_tree cimport *
from gudhi.simplex_tree import SimplexTree

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

cdef extern from "Rips_complex_interface.h" namespace "Gudhi":
    cdef cppclass Rips_complex_interface "Gudhi::rips_complex::Rips_complex_interface":
        Rips_complex_interface() nogil
        void init_points(vector[vector[double]] values, double threshold) nogil
        void init_matrix(vector[vector[double]] values, double threshold) nogil
        void init_points_sparse(vector[vector[double]] values, double threshold, double sparse) nogil
        void init_matrix_sparse(vector[vector[double]] values, double threshold, double sparse) nogil
        void create_simplex_tree(Simplex_tree_interface_full_featured* simplex_tree, int dim_max) nogil except +

# RipsComplex python interface
cdef class RipsComplex:
    """The data structure is a one skeleton graph, or Rips graph, containing edges when the edge length is less or
    equal to a given threshold. Edge length is computed from a user given point cloud with a given distance function,
    or a distance matrix.
    """

    cdef Rips_complex_interface thisref

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, *, points=None, distance_matrix=None, max_edge_length=float('inf'), sparse=None):
        """RipsComplex constructor.

        :param points: A list of points in d-Dimension.
        :type points: List[List[float]]

        Or

        :param distance_matrix: A distance matrix (full square or lower triangular).
        :type distance_matrix: List[List[float]]

        And in both cases

        :param max_edge_length: Maximal edge length. All edges of the graph strictly greater than `threshold` are not
            inserted in the graph.
        :type max_edge_length: float
        :param sparse: If this is not None, it switches to building a sparse Rips and represents the approximation
            parameter epsilon.
        :type sparse: float
        """

    # The real cython constructor
    def __cinit__(self, *, points=None, distance_matrix=None, max_edge_length=float('inf'), sparse=None):
        if sparse is not None:
          if distance_matrix is not None:
              self.thisref.init_matrix_sparse(distance_matrix, max_edge_length, sparse)
          else:
              if points is None:
                  # Empty Rips construction
                  points=[]
              self.thisref.init_points_sparse(points, max_edge_length, sparse)
        else:
          if distance_matrix is not None:
              self.thisref.init_matrix(distance_matrix, max_edge_length)
          else:
              if points is None:
                  # Empty Rips construction
                  points=[]
              self.thisref.init_points(points, max_edge_length)


    def create_simplex_tree(self, max_dimension=1):
        """
        :param max_dimension: graph expansion for Rips until this given maximal dimension.
        :type max_dimension: int
        :returns: A simplex tree encoding the Vietoris–Rips filtration.
        :rtype: SimplexTree
        """
        stree = SimplexTree()
        cdef intptr_t stree_int_ptr=stree.thisptr
        cdef int maxdim = max_dimension
        with nogil:
            self.thisref.create_simplex_tree(<Simplex_tree_interface_full_featured*>stree_int_ptr, maxdim)
        return stree
