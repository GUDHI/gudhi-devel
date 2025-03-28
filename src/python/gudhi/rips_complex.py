# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from gudhi import _rips_complex_ext as t
from gudhi.simplex_tree import SimplexTree

__author__ = "Thibaud Kloczko"
__copyright__ = "Copyright (C) 2025 Inria"
__license__ = "MIT"

# cdef extern from "Rips_complex_interface.h" namespace "Gudhi":
#     cdef cppclass Rips_complex_interface "Gudhi::rips_complex::Rips_complex_interface":
#         Rips_complex_interface() nogil
#         void init_points(vector[vector[double]] values, double threshold) nogil
#         void init_matrix(vector[vector[double]] values, double threshold) nogil
#         void init_points_sparse(vector[vector[double]] values, double threshold, double sparse) nogil
#         void init_matrix_sparse(vector[vector[double]] values, double threshold, double sparse) nogil
#         void create_simplex_tree(Simplex_tree_python_interface* simplex_tree, int dim_max) nogil except +

# RipsComplex python interface
class RipsComplex(t.Rips_complex_interfcae):
    """The data structure is a one skeleton graph, or Rips graph, containing edges when the edge length is less or
    equal to a given threshold. Edge length is computed from a user given point cloud with a given distance function,
    or a distance matrix.
    """

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
        super().__init__()
        if sparse is not None:
          if distance_matrix is not None:
              super().init_matrix_sparse(distance_matrix, max_edge_length, sparse)
          else:
              if points is None:
                  # Empty Rips construction
                  points=[]
              super().init_points_sparse(points, max_edge_length, sparse)
        else:
          if distance_matrix is not None:
              super().init_matrix(distance_matrix, max_edge_length)
          else:
              if points is None:
                  # Empty Rips construction
                  points=[]
              super().init_points(points, max_edge_length)super()


    def create_simplex_tree(self, max_dimension=1):
        """
        :param max_dimension: graph expansion for Rips until this given maximal dimension.
        :type max_dimension: int
        :returns: A simplex tree encoding the Vietoris–Rips filtration.
        :rtype: SimplexTree
        """
        stree = SimplexTree()
        super().create_simplex_tree(stree, maxdim)
        return stree
