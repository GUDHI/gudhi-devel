from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Vincent Rouvreau

   Copyright (C) 2016 INRIA

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
__copyright__ = "Copyright (C) 2016 INRIA"
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
        simplex_tree = SimplexTree()
        if limit_dimension is not -1:
            self.thisptr.create_simplex_tree(simplex_tree.thisptr, max_alpha_square, limit_dimension)
        else:
            self.thisptr.create_simplex_tree(simplex_tree.thisptr, max_alpha_square)
        return simplex_tree

    def get_point(self, vertex):
        """This function returns the point corresponding to a given vertex.

        :param vertex: The vertex.
        :type vertex: int.
        :returns:  The point.
        :rtype: list of float
        """
        cdef vector[double] point = self.thisptr.get_point(vertex)
        return point

