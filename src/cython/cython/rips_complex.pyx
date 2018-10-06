from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair
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

cdef extern from "Rips_complex_interface.h" namespace "Gudhi":
    cdef cppclass Rips_complex_interface "Gudhi::rips_complex::Rips_complex_interface":
        Rips_complex_interface()
        void init_points(vector[vector[double]] values, double threshold)
        void init_matrix(vector[vector[double]] values, double threshold)
        void init_points_sparse(vector[vector[double]] values, double threshold, double sparse)
        void init_matrix_sparse(vector[vector[double]] values, double threshold, double sparse)
        void create_simplex_tree(Simplex_tree_interface_full_featured* simplex_tree, int dim_max)

# RipsComplex python interface
cdef class RipsComplex:
    """The data structure is a one skeleton graph, or Rips graph, containing
    edges when the edge length is less or equal to a given threshold. Edge
    length is computed from a user given point cloud with a given distance
    function, or a distance matrix.
    """

    cdef Rips_complex_interface thisref

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, points=None, distance_matrix=None, max_edge_length=float('inf'), sparse=None):
        """RipsComplex constructor.

        :param max_edge_length: Rips value.
        :type max_edge_length: float

        :param points: A list of points in d-Dimension.
        :type points: list of list of double

        Or

        :param distance_matrix: A distance matrix (full square or lower triangular).
        :type points: list of list of double

        And in both cases

        :param sparse: If this is not None, it switches to building a sparse Rips and represents the approximation parameter epsilon.
        :type sparse: float
        """

    # The real cython constructor
    def __cinit__(self, points=None, distance_matrix=None, max_edge_length=float('inf'), sparse=None):
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
        :param max_dimension: graph expansion for rips until this given maximal
            dimension.
        :type max_dimension: int
        :returns: A simplex tree created from the Delaunay Triangulation.
        :rtype: SimplexTree
        """
        simplex_tree = SimplexTree()
        self.thisref.create_simplex_tree(simplex_tree.thisptr, max_dimension)
        return simplex_tree
