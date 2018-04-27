from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair

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

cdef extern from "Strong_witness_complex_interface.h" namespace "Gudhi":
    cdef cppclass Strong_witness_complex_interface "Gudhi::witness_complex::Strong_witness_complex_interface":
        Strong_witness_complex_interface(vector[vector[pair[size_t, double]]] nearest_landmark_table)
        void create_simplex_tree(Simplex_tree_interface_full_featured* simplex_tree, double max_alpha_square)
        void create_simplex_tree(Simplex_tree_interface_full_featured* simplex_tree, double max_alpha_square,
            unsigned limit_dimension)

# StrongWitnessComplex python interface
cdef class StrongWitnessComplex:
    """Constructs (strong) witness complex for a given table of nearest
    landmarks with respect to witnesses.
    """

    cdef Strong_witness_complex_interface * thisptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, nearest_landmark_table=None):
        """StrongWitnessComplex constructor.

        :param nearest_landmark_table: A list of nearest landmark.
        :type nearest_landmark_table: list of list of pair of unsigned and double
        """

    # The real cython constructor
    def __cinit__(self, nearest_landmark_table=None):
        if nearest_landmark_table is not None:
            self.thisptr = new Strong_witness_complex_interface(nearest_landmark_table)

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr

    def __is_defined(self):
        """Returns true if StrongWitnessComplex pointer is not NULL.
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
        simplex_tree = SimplexTree()
        if limit_dimension is not -1:
            self.thisptr.create_simplex_tree(simplex_tree.thisptr, max_alpha_square, limit_dimension)
        else:
            self.thisptr.create_simplex_tree(simplex_tree.thisptr, max_alpha_square)
        return simplex_tree
