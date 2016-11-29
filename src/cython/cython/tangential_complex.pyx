from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp cimport bool
from cython.operator cimport dereference as deref
import os

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

cdef extern from "Tangential_complex_interface.h" namespace "Gudhi":
    cdef cppclass Tangential_complex_interface "Gudhi::tangential_complex::Tangential_complex_interface":
        Tangential_complex_interface(vector[vector[double]] points)
        # bool from_file is a workaround for cython to find the correct signature
        Tangential_complex_interface(string off_file, bool from_file)
        vector[double] get_point(unsigned vertex)
        unsigned number_of_vertices()
        unsigned number_of_simplices()
        unsigned number_of_inconsistent_simplices()
        unsigned number_of_inconsistent_stars()
        void create_simplex_tree(Simplex_tree_interface_full_featured* simplex_tree)

# TangentialComplex python interface
cdef class TangentialComplex:
    """TangentialComplex is a simplicial complex constructed from the finite cells
    of a Delaunay Triangulation.

    The filtration value of each simplex is computed as the square of the
    circumradius of the simplex if the circumsphere is empty (the simplex is
    then said to be Gabriel), and as the minimum of the filtration values of
    the codimension 1 cofaces that make it not Gabriel otherwise.

    All simplices that have a filtration value strictly greater than a given
    alpha squared value are not inserted into the complex.

    .. note::

        When Tangential_complex is constructed with an infinite value of alpha, the
        complex is a Delaunay complex.

    """

    cdef Tangential_complex_interface * thisptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, points=None, off_file=''):
        """TangentialComplex constructor.

        :param points: A list of points in d-Dimension.
        :type points: list of list of double

        Or

        :param off_file: An OFF file style name.
        :type off_file: string
        """

    # The real cython constructor
    def __cinit__(self, points=[], off_file=''):
        if off_file is not '':
            if os.path.isfile(off_file):
                self.thisptr = new Tangential_complex_interface(off_file, True)
            else:
                print("file " + off_file + " not found.")
        else:
            self.thisptr = new Tangential_complex_interface(points)
                

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr

    def __is_defined(self):
        """Returns true if TangentialComplex pointer is not NULL.
         """
        return self.thisptr != NULL

    def get_point(self, vertex):
        """This function returns the point corresponding to a given vertex.

        :param vertex: The vertex.
        :type vertex: int.
        :returns:  list of float -- the point.
        """
        cdef vector[double] point = self.thisptr.get_point(vertex)
        return point

    def num_vertices(self):
        """This function returns the number of vertices.

        :returns:  unsigned -- the number of vertices.
        """
        return self.thisptr.number_of_vertices()

    def num_simplices(self):
        """This function returns the number of simplices.

        :returns:  unsigned -- the number of simplices.
        """
        return self.thisptr.number_of_simplices()

    def num_inconsistent_simplices(self):
        """This function returns the number of inconsistent simplices.

        :returns:  unsigned -- the number of inconsistent simplices.
        """
        return self.thisptr.number_of_inconsistent_simplices()

    def num_inconsistent_stars(self):
        """This function returns the number of inconsistent stars.

        :returns:  unsigned -- the number of inconsistent stars.
        """
        return self.thisptr.number_of_inconsistent_stars()

    def create_simplex_tree(self):
        """This function creates the given simplex tree from the Delaunay
        Triangulation.

        :param simplex_tree: The simplex tree to create (must be empty)
        :type simplex_tree: SimplexTree
        :returns: A simplex tree created from the Delaunay Triangulation.
        :rtype: SimplexTree
        """
        simplex_tree = SimplexTree()
        self.thisptr.create_simplex_tree(simplex_tree.thisptr)
        return simplex_tree
