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

cdef extern from "Witness_complex_interface.h" namespace "Gudhi":
    cdef cppclass Witness_complex_interface "Gudhi::witness_complex::Witness_complex_interface":
        Witness_complex_interface(vector[vector[pair[size_t, double]]] nearest_landmark_table)
        void create_simplex_tree(Simplex_tree_interface_full_featured* simplex_tree, double max_alpha_square)
        void create_simplex_tree(Simplex_tree_interface_full_featured* simplex_tree, double max_alpha_square,
            unsigned limit_dimension)

# WitnessComplex python interface
cdef class WitnessComplex:
    """Constructs (weak) witness complex for a given table of nearest landmarks
    with respect to witnesses.
    """

    cdef Witness_complex_interface * thisptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, nearest_landmark_table=None):
        """WitnessComplex constructor.

        :param nearest_landmark_table: A list of lists of nearest landmarks and their distances.
            `nearest_landmark_table[w][k]==(l,d)` means that l is the k-th nearest landmark to
            witness w, and d is the (squared) distance between l and w.
        :type nearest_landmark_table: list of list of pair of int and float
        """

    # The real cython constructor
    def __cinit__(self, nearest_landmark_table=None):
        if nearest_landmark_table is not None:
            self.thisptr = new Witness_complex_interface(nearest_landmark_table)

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr

    def __is_defined(self):
        """Returns true if WitnessComplex pointer is not NULL.
         """
        return self.thisptr != NULL

    def create_simplex_tree(self, max_alpha_square = float('inf'), limit_dimension = -1):
        """
        :param max_alpha_square: The maximum relaxation parameter.
            Default is set to infinity.
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
