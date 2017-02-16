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

cdef extern from "Rips_complex_interface.h" namespace "Gudhi":
    cdef cppclass Rips_complex_interface "Gudhi::rips_complex::Rips_complex_interface":
        Rips_complex_interface(vector[vector[double]] values, double threshold, bool euclidean)
        # bool from_file is a workaround for cython to find the correct signature
        Rips_complex_interface(string file_name, double threshold, bool euclidean, bool from_file)
        void create_simplex_tree(Simplex_tree_interface_full_featured* simplex_tree, int dim_max)

# RipsComplex python interface
cdef class RipsComplex:
    """The data structure is a one skeleton graph, or Rips graph, containing
    edges when the edge length is less or equal to a given threshold. Edge
    length is computed from a user given point cloud with a given distance
    function, or a distance matrix.
    """

    cdef Rips_complex_interface * thisptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, points=None, off_file='', distance_matrix=None, csv_file='', max_edge_length=float('inf')):
        """RipsComplex constructor.

        :param max_edge_length: Rips value.
        :type max_edge_length: int

        :param points: A list of points in d-Dimension.
        :type points: list of list of double

        Or

        :param off_file: An OFF file style name.
        :type off_file: string

        Or

        :param distance_matrix: A distance matrix (full square or lower
        triangular).
        :type points: list of list of double

        Or

        :param csv_file: A csv file style name containing a full square or a
        lower triangular distance matrix.
        :type csv_file: string
        """

    # The real cython constructor
    def __cinit__(self, points=None, off_file='', distance_matrix=None, csv_file='', max_edge_length=float('inf')):
        if off_file is not '':
            if os.path.isfile(off_file):
                self.thisptr = new Rips_complex_interface(off_file,
                                                          max_edge_length,
                                                          True,
                                                          True)
            else:
                print("file " + off_file + " not found.")
        elif csv_file is not '':
            if os.path.isfile(csv_file):
                self.thisptr = new Rips_complex_interface(csv_file,
                                                          max_edge_length,
                                                          False,
                                                          True)
            else:
                print("file " + csv_file + " not found.")
        elif distance_matrix is not None:
            self.thisptr = new Rips_complex_interface(distance_matrix, max_edge_length, False)
        else:
            if points is None:
                # Empty Rips construction
                points=[]
            self.thisptr = new Rips_complex_interface(points, max_edge_length, True)


    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr

    def __is_defined(self):
        """Returns true if RipsComplex pointer is not NULL.
         """
        return self.thisptr != NULL

    def create_simplex_tree(self, max_dimension=1):
        """
        :param max_dimension: graph expansion for rips until this given maximal
            dimension.
        :type max_dimension: int
        :returns: A simplex tree created from the Delaunay Triangulation.
        :rtype: SimplexTree
        """
        simplex_tree = SimplexTree()
        self.thisptr.create_simplex_tree(simplex_tree.thisptr, max_dimension)
        return simplex_tree
