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

   Copyright (C) 2018 Inria

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
__copyright__ = "Copyright (C) 2018 Inria"
__license__ = "GPL v3"

cdef extern from "Nerve_gic_interface.h" namespace "Gudhi":
    cdef cppclass Nerve_gic_interface "Gudhi::cover_complex::Nerve_gic_interface":
        Nerve_gic_interface()
        double compute_confidence_level_from_distance(double distance)
        double compute_distance_from_confidence_level(double alpha)
        void compute_distribution(int N)
        double compute_p_value()
        void compute_PD()
        void create_simplex_tree(Simplex_tree_interface_full_featured* simplex_tree)
        bool read_point_cloud(string off_file_name)
        double set_automatic_resolution()
        void set_color_from_coordinate(int k)
        void set_color_from_file(string color_file_name)
        void set_color_from_vector(vector[double] color)
        void set_cover_from_file(string cover_file_name)
        void set_cover_from_function()
        void set_cover_from_Euclidean_Voronoi(int m)
        void set_function_from_coordinate(int k)
        void set_function_from_file(string func_file_name)
        void set_function_from_range(vector[double] function)
        void set_gain(double g)

# CoverComplex python interface
cdef class CoverComplex:
    """Cover complex data structure.

    The data structure is a simplicial complex, representing a Graph Induced
    simplicial Complex (GIC) or a Nerve, and whose simplices are computed with
    a cover C of a point cloud P, which often comes from the preimages of
    intervals covering the image of a function f defined on P. These intervals
    are parameterized by their resolution (either their length or their number)
    and their gain (percentage of overlap). To compute a GIC, one also needs a
    graph G built on top of P, whose cliques with vertices belonging to
    different elements of C correspond to the simplices of the GIC.
    """

    cdef Nerve_gic_interface * thisptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self):
        """CoverComplex constructor.
        """

    # The real cython constructor
    def __cinit__(self):
        self.thisptr = new Nerve_gic_interface()

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr

    def __is_defined(self):
        """Returns true if CoverComplex pointer is not NULL.
         """
        return self.thisptr != NULL

    def compute_confidence_level_from_distance(self, distance):
        """Computes the confidence level of a specific bottleneck distance
        threshold.

        :param distance: Bottleneck distance.
        :type distance: double
        :rtype: double
        :returns: Confidence level.
        """
        return self.thisptr.compute_confidence_level_from_distance(distance)

    def compute_distance_from_confidence_level(self, alpha):
        """Computes the bottleneck distance threshold corresponding to a
        specific confidence level.

        :param alpha: Confidence level.
        :type alpha: double
        :rtype: double
        :returns: Bottleneck distance.
        """
        return self.thisptr.compute_distance_from_confidence_level(alpha)

    def compute_distribution(self, N=100):
        """Computes bootstrapped distances distribution.

        :param N: Loop number (default value is 100).
        :type alpha: int
        """
        self.thisptr.compute_distribution(N)

    def compute_p_value(self):
        """Computes the p-value, i.e. the opposite of the confidence level of
        the largest bottleneck distance preserving the points in the
        persistence diagram of the output simplicial complex.

        :rtype: double
        :returns: p-value.
        """
        return self.thisptr.compute_p_value()

    def compute_PD(self):
        """Computes the extended persistence diagram of the complex.
        """
        self.thisptr.compute_PD()

    def create_simplex_tree(self):
        """
        :returns: A simplex tree created from the Cover complex.
        :rtype: SimplexTree
        """
        simplex_tree = SimplexTree()
        self.thisptr.create_simplex_tree(simplex_tree.thisptr)
        return simplex_tree

    def read_point_cloud(self, off_file):
        """Reads and stores the input point cloud.

        :param off_file: Name of the input .OFF or .nOFF file.
        :type off_file: string
        :rtype: bool
        :returns: Read file status.
        """
        if os.path.isfile(off_file):
            return self.thisptr.read_point_cloud(str.encode(off_file))
        else:
            print("file " + off_file + " not found.")
            return False

    def set_automatic_resolution(self):
        """Computes the optimal length of intervals (i.e. the smallest interval
        length avoiding discretization artifacts—see [12]) for a functional
        cover.

        :rtype: double
        :returns: reso interval length used to compute the cover.
        """
        return self.thisptr.set_automatic_resolution()

    def set_color_from_coordinate(self, k=0):
        """Computes the function used to color the nodes of the simplicial
        complex from the k-th coordinate.

        :param k: Coordinate to use (start at 0). Default value is 0.
        :type k: int
        """
        return self.thisptr.set_color_from_coordinate(k)

    def set_color_from_file(self, color_file_name):
        """Computes the function used to color the nodes of the simplicial
        complex from a file containing the function values.

        :param color_file_name: Name of the input color file.
        :type color_file_name: string
        """
        if os.path.isfile(color_file_name):
            self.thisptr.set_color_from_file(str.encode(color_file_name))
        else:
            print("file " + color_file_name + " not found.")

    def set_color_from_vector(self, color):
        """Computes the function used to color the nodes of the simplicial
        complex from a vector stored in memory.

        :param color: Input vector of values.
        :type color: vector[double]
        """
        self.thisptr.set_color_from_vector(color)

    def set_cover_from_file(self, cover_file_name):
        """Creates the cover C from a file containing the cover elements of
        each point (the order has to be the same as in the input file!).

        :param cover_file_name: Name of the input cover file.
        :type cover_file_name: string
        """
        if os.path.isfile(cover_file_name):
            self.thisptr.set_cover_from_file(str.encode(cover_file_name))
        else:
            print("file " + cover_file_name + " not found.")

    def set_cover_from_function(self):
        """Creates a cover C from the preimages of the function f.
        """
        self.thisptr.set_cover_from_function()

    def set_cover_from_Voronoi(self, m=100):
        """Creates the cover C from the Voronoï cells of a subsampling of the
        point cloud.

        :param m: Number of points in the subsample. Default value is 100.
        :type m: int
        """
        self.thisptr.set_cover_from_Euclidean_Voronoi(m)

    def set_function_from_coordinate(self, k):
        """Creates the function f from the k-th coordinate of the point cloud.

        :param k: coordinate to use (start at 0).
        :type k: int
        """
        self.thisptr.set_function_from_coordinate(k)

    def set_function_from_file(self, func_file_name):
        """Creates the function f from a file containing the function values.

        :param func_file_name: name of the input function file.
        :type func_file_name: string
        """
        if os.path.isfile(func_file_name):
            self.thisptr.set_function_from_file(str.encode(func_file_name))
        else:
            print("file " + func_file_name + " not found.")

    def set_function_from_range(self, function):
        """Creates the function f from a vector stored in memory.

        :param function: Input vector of values.
        :type function: vector[double]
        """
        self.thisptr.set_function_from_range(function)

    def set_gain(self, g):
        """Sets a gain from a value stored in memory.

        :param g: gain (default value is 0.3).
        :type g: double
        """
        self.thisptr.set_gain(g)
