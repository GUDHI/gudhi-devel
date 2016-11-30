from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp cimport bool

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

cdef extern from "Simplex_tree_interface.h" namespace "Gudhi":
    cdef cppclass Simplex_tree_options_full_featured:
        pass

    cdef cppclass Simplex_tree_interface_full_featured "Gudhi::Simplex_tree_interface<Gudhi::Simplex_tree_options_full_featured>":
        Simplex_tree()
        double filtration()
        double simplex_filtration(vector[int] simplex)
        void set_filtration(double filtration)
        void initialize_filtration()
        int num_vertices()
        int num_simplices()
        void set_dimension(int dimension)
        int dimension()
        bint find_simplex(vector[int] simplex)
        bint insert_simplex_and_subfaces(vector[int] simplex,
                                         double filtration)
        vector[pair[vector[int], double]] get_filtered_tree()
        vector[pair[vector[int], double]] get_skeleton_tree(int dimension)
        vector[pair[vector[int], double]] get_star_tree(vector[int] simplex)
        vector[pair[vector[int], double]] get_coface_tree(vector[int] simplex,
                                                          int dimension)
        void remove_maximal_simplex(vector[int] simplex)

cdef extern from "Persistent_cohomology_interface.h" namespace "Gudhi":
    cdef cppclass Simplex_tree_persistence_interface "Gudhi::Persistent_cohomology_interface<Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_full_featured>>":
        Simplex_tree_persistence_interface(Simplex_tree_interface_full_featured * st, bool persistence_dim_max)
        vector[pair[int, pair[double, double]]] get_persistence(int homology_coeff_field, double min_persistence)
        vector[int] betti_numbers()
        vector[int] persistent_betti_numbers(double from_value, double to_value)

# SimplexTree python interface
cdef class SimplexTree:
    """The simplex tree is an efficient and flexible data structure for
    representing general (filtered) simplicial complexes. The data structure
    is described in Jean-Daniel Boissonnat and Clément Maria. The Simplex
    Tree: An Efficient Data Structure for General Simplicial Complexes.
    Algorithmica, pages 1–22, 2014.

    This class is a filtered, with keys, and non contiguous vertices version
    of the simplex tree.
    """
    cdef Simplex_tree_interface_full_featured * thisptr

    cdef Simplex_tree_persistence_interface * pcohptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self):
        """SimplexTree constructor.
        """

    # The real cython constructor
    def __cinit__(self):
        self.thisptr = new Simplex_tree_interface_full_featured()

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr
        if self.pcohptr != NULL:
            del self.pcohptr

    def __is_defined(self):
        """Returns true if SimplexTree pointer is not NULL.
         """
        return self.thisptr != NULL

    def __is_persistence_defined(self):
        """Returns true if Persistence pointer is not NULL.
         """
        return self.pcohptr != NULL

    def get_filtration(self):
        """This function returns the main simplicial complex filtration value.

        :returns:  The simplicial complex filtration value.
        :rtype:  float
        """
        return self.thisptr.filtration()

    def filtration(self, simplex):
        """This function returns the simplicial complex filtration value for a
        given N-simplex.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int.
        :returns:  The simplicial complex filtration value.
        :rtype:  float
        """
        return self.thisptr.simplex_filtration(simplex)

    def set_filtration(self, filtration):
        """This function sets the main simplicial complex filtration value.

        :param filtration: The filtration value.
        :type filtration: float.
        """
        self.thisptr.set_filtration(<double> filtration)

    def initialize_filtration(self):
        """This function initializes and sorts the simplicial complex
        filtration vector.

        .. note::

            This function must be launched before persistence, betti_numbers,
            persistent_betti_numbers or get_filtered_tree after inserting or
            removing simplices.
        """
        self.thisptr.initialize_filtration()

    def num_vertices(self):
        """This function returns the number of vertices of the simplicial
        complex.

        :returns:  The simplicial complex number of vertices.
        :rtype:  int
        """
        return self.thisptr.num_vertices()

    def num_simplices(self):
        """This function returns the number of simplices of the simplicial
        complex.

        :returns:  the simplicial complex number of simplices.
        :rtype:  int
        """
        return self.thisptr.num_simplices()

    def dimension(self):
        """This function returns the dimension of the simplicial complex.

        :returns:  the simplicial complex dimension.
        :rtype:  int
        """
        return self.thisptr.dimension()

    def set_dimension(self, dimension):
        """This function sets the dimension of the simplicial complex.

        :param dimension: The new dimension value.
        :type dimension: int.
        """
        self.thisptr.set_dimension(<int>dimension)

    def find(self, simplex):
        """This function returns if the N-simplex was found in the simplicial
        complex or not.

        :param simplex: The N-simplex to find, represented by a list of vertex.
        :type simplex: list of int.
        :returns:  true if the simplex was found, false otherwise.
        :rtype:  bool
        """
        cdef vector[int] complex
        for i in simplex:
            complex.push_back(i)
        return self.thisptr.find_simplex(complex)

    def insert(self, simplex, filtration=0.0):
        """This function inserts the given N-simplex with the given filtration
        value (default value is '0.0').

        :param simplex: The N-simplex to insert, represented by a list of
            vertex.
        :type simplex: list of int.
        :param filtration: The filtration value of the simplex.
        :type filtration: float.
        :returns:  true if the simplex was found, false otherwise.
        :rtype:  bool
        """
        cdef vector[int] complex
        for i in simplex:
            complex.push_back(i)
        return self.thisptr.insert_simplex_and_subfaces(complex,
                                                        <double>filtration)

    def get_filtered_tree(self):
        """This function returns the tree sorted by increasing filtration
        values.

        :returns:  The tree sorted by increasing filtration values.
        :rtype:  list of tuples(simplex, filtration)
        """
        cdef vector[pair[vector[int], double]] coface_tree \
            = self.thisptr.get_filtered_tree()
        ct = []
        for filtered_complex in coface_tree:
            v = []
            for vertex in filtered_complex.first:
                v.append(vertex)
            ct.append((v, filtered_complex.second))
        return ct

    def get_skeleton_tree(self, dimension):
        """This function returns the tree skeleton of a maximum given
        dimension.

        :param dimension: The skeleton dimension value.
        :type dimension: int.
        :returns:  The skeleton tree of a maximum dimension.
        :rtype:  list of tuples(simplex, filtration)
        """
        cdef vector[pair[vector[int], double]] coface_tree \
            = self.thisptr.get_skeleton_tree(<int>dimension)
        ct = []
        for filtered_complex in coface_tree:
            v = []
            for vertex in filtered_complex.first:
                v.append(vertex)
            ct.append((v, filtered_complex.second))
        return ct

    def get_star_tree(self, simplex):
        """This function returns the star tree of a given N-simplex.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int.
        :returns:  The star tree of a simplex.
        :rtype:  list of tuples(simplex, filtration)
        """
        cdef vector[int] complex
        for i in simplex:
            complex.push_back(i)
        cdef vector[pair[vector[int], double]] coface_tree \
            = self.thisptr.get_star_tree(complex)
        ct = []
        for filtered_complex in coface_tree:
            v = []
            for vertex in filtered_complex.first:
                v.append(vertex)
            ct.append((v, filtered_complex.second))
        return ct

    def get_coface_tree(self, simplex, codimension):
        """This function returns the coface tree of a given N-simplex with a
        given codimension.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int.
        :param codimension: The codimension. If codimension = 0, all cofaces
            are returned (equivalent of get_star_tree function)
        :type codimension: int.
        :returns:  The coface tree of a simplex.
        :rtype:  list of tuples(simplex, filtration)
        """
        cdef vector[int] complex
        for i in simplex:
            complex.push_back(i)
        cdef vector[pair[vector[int], double]] coface_tree \
            = self.thisptr.get_coface_tree(complex, <int>codimension)
        ct = []
        for filtered_complex in coface_tree:
            v = []
            for vertex in filtered_complex.first:
                v.append(vertex)
            ct.append((v, filtered_complex.second))
        return ct

    def remove_maximal_simplex(self, simplex):
        """This function removes a given maximal N-simplex from the simplicial
        complex.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int.
        """
        self.thisptr.remove_maximal_simplex(simplex)

    def persistence(self, homology_coeff_field=11, min_persistence=0, persistence_dim_max = False):
        """This function returns the persistence of the simplicial complex.

        :param homology_coeff_field: The homology coefficient field. Must be a
            prime number
        :type homology_coeff_field: int.
        :param min_persistence: The minimum persistence value to take into
            account (strictly greater than min_persistence). Default value is
            0.0.
            Sets min_persistence to -1.0 to see all values.
        :type min_persistence: float.
        :returns: The persistence of the simplicial complex.
        :rtype:  list of pairs(dimension, pair(birth, death))
        """
        if self.pcohptr != NULL:
            del self.pcohptr
        self.pcohptr = new Simplex_tree_persistence_interface(self.thisptr, persistence_dim_max)
        cdef vector[pair[int, pair[double, double]]] persistence_result
        if self.pcohptr != NULL:
            persistence_result = self.pcohptr.get_persistence(homology_coeff_field, min_persistence)
        return persistence_result

    def betti_numbers(self):
        """This function returns the Betti numbers of the simplicial complex.

        :returns: The Betti numbers ([B0, B1, ..., Bn]).
        :rtype:  list of int

        :note: betti_numbers function requires persistence function to be
            launched first.
        """
        cdef vector[int] bn_result
        if self.pcohptr != NULL:
            bn_result = self.pcohptr.betti_numbers()
        else:
            print("betti_numbers function requires persistence function"
                  " to be launched first.")
        return bn_result

    def persistent_betti_numbers(self, from_value, to_value):
        """This function returns the persistent Betti numbers of the
        simplicial complex.

        :param from_value: The persistence birth limit to be added in the
            numbers (persistent birth <= from_value).
        :type from_value: float.
        :param to_value: The persistence death limit to be added in the
            numbers (persistent death > to_value).
        :type to_value: float.

        :returns: The persistent Betti numbers ([B0, B1, ..., Bn]).
        :rtype:  list of int

        :note: persistent_betti_numbers function requires persistence
            function to be launched first.
        """
        cdef vector[int] pbn_result
        if self.pcohptr != NULL:
            pbn_result = self.pcohptr.persistent_betti_numbers(<double>from_value, <double>to_value)
        else:
            print("persistent_betti_numbers function requires persistence function"
                  " to be launched first.")
        return pbn_result
