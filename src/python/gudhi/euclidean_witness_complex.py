# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - 2025/03 Vincent Rouvreau: Use nanobind instead of Cython for python bindings.
#   - YYYY/MM Author: Description of the modification

from gudhi import _euclidean_witness_complex_ext as t
from gudhi.simplex_tree import SimplexTree

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "GPL v3"

# cdef extern from "Euclidean_witness_complex_interface.h" namespace "Gudhi":
#     cdef cppclass Euclidean_witness_complex_interface "Gudhi::witness_complex::Euclidean_witness_complex_interface":
#         Euclidean_witness_complex_interface(vector[vector[double]] landmarks, vector[vector[double]] witnesses)
#         void create_simplex_tree(Simplex_tree_python_interface* simplex_tree, double max_alpha_square) except +
#         void create_simplex_tree(Simplex_tree_python_interface* simplex_tree, double max_alpha_square,
#             unsigned limit_dimension) except +
#         vector[double] get_point(unsigned vertex)

# EuclideanWitnessComplex python interface
class EuclideanWitnessComplex(t.Euclidean_witness_complex_interface):
    """Constructs (weak) witness complex for given sets of witnesses and
    landmarks in Euclidean space.
    """
    def __init__(self, landmarks, witnesses):
        """WitnessComplex constructor.

        :param landmarks: A list of landmarks (in the point cloud).
        :type landmarks: list of list of double

        :param witnesses: The point cloud.
        :type witnesses: list of list of double
        """
        super().__init__(landmarks, witnesses)

    def create_simplex_tree(self, max_alpha_square, limit_dimension = -1):
        """
        :param max_alpha_square: The maximum alpha square threshold the
            simplices shall not exceed. Default is set to infinity.
        :type max_alpha_square: float
        :returns: A simplex tree created from the Delaunay Triangulation.
        :rtype: SimplexTree
        """
        stree = SimplexTree()
        if limit_dimension != -1:
            super().create_simplex_tree(stree, max_alpha_square, limit_dimension)
        else:
            super().create_simplex_tree(stree, max_alpha_square)
        return stree

    def get_point(self, vertex):
        """This function returns the point corresponding to a given vertex.

        :param vertex: The vertex.
        :type vertex: int.
        :returns:  The point.
        :rtype: list of float
        """
        return super().get_point(vertex)
