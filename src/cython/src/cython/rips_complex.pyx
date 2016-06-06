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

cdef extern from "Simplex_tree_interface.h" namespace "Gudhi":
    cdef cppclass Simplex_tree_options_full_featured:
        pass

    cdef cppclass Rips_complex_interface "Gudhi::Simplex_tree_interface<Gudhi::Simplex_tree_options_full_featured>":
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
        void graph_expansion(vector[vector[double]] points, int max_dimension,
                             double max_edge_length)

cdef extern from "Persistent_cohomology_interface.h" namespace "Gudhi":
    cdef cppclass Rips_complex_persistence_interface "Gudhi::Persistent_cohomology_interface<Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_full_featured>>":
        Rips_complex_persistence_interface(Rips_complex_interface * st)
        vector[pair[int, pair[double, double]]] get_persistence(int homology_coeff_field, double min_persistence)
        vector[int] betti_numbers()
        vector[int] persistent_betti_numbers(double from_value, double to_value)

# RipsComplex python interface
cdef class RipsComplex:
    cdef Rips_complex_interface * thisptr

    cdef Rips_complex_persistence_interface * pcohptr

    def __cinit__(self, points=None, max_dimension=3,
                  max_edge_length=float('inf')):
        self.thisptr = new Rips_complex_interface()
        # Constructor from graph expansion
        if points is not None:
            self.thisptr.graph_expansion(points, max_dimension,
                                         max_edge_length)

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr
        if self.pcohptr != NULL:
            del self.pcohptr

    def get_filtration(self):
        return self.thisptr.filtration()

    def filtration(self, simplex):
        return self.thisptr.simplex_filtration(simplex)

    def set_filtration(self, filtration):
        self.thisptr.set_filtration(<double> filtration)

    def initialize_filtration(self):
        self.thisptr.initialize_filtration()

    def num_vertices(self):
        return self.thisptr.num_vertices()

    def num_simplices(self):
        return self.thisptr.num_simplices()

    def dimension(self):
        return self.thisptr.dimension()

    def set_dimension(self, dim):
        self.thisptr.set_dimension(<int>dim)

    def find(self, simplex):
        cdef vector[int] complex
        for i in simplex:
            complex.push_back(i)
        return self.thisptr.find_simplex(complex)

    def insert(self, simplex, filtration=0.0):
        cdef vector[int] complex
        for i in simplex:
            complex.push_back(i)
        return self.thisptr.insert_simplex_and_subfaces(complex,
                                                        <double>filtration)

    def get_filtered_tree(self):
        cdef vector[pair[vector[int], double]] coface_tree \
            = self.thisptr.get_filtered_tree()
        ct = []
        for filtered_complex in coface_tree:
            v = []
            for vertex in filtered_complex.first:
                v.append(vertex)
            ct.append((v, filtered_complex.second))
        return ct

    def get_skeleton_tree(self, dim):
        cdef vector[pair[vector[int], double]] coface_tree \
            = self.thisptr.get_skeleton_tree(<int>dim)
        ct = []
        for filtered_complex in coface_tree:
            v = []
            for vertex in filtered_complex.first:
                v.append(vertex)
            ct.append((v, filtered_complex.second))
        return ct

    def get_star_tree(self, simplex):
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

    def get_coface_tree(self, simplex, dim):
        cdef vector[int] complex
        for i in simplex:
            complex.push_back(i)
        cdef vector[pair[vector[int], double]] coface_tree \
            = self.thisptr.get_coface_tree(complex, <int>dim)
        ct = []
        for filtered_complex in coface_tree:
            v = []
            for vertex in filtered_complex.first:
                v.append(vertex)
            ct.append((v, filtered_complex.second))
        return ct

    def remove_maximal_simplex(self, simplex):
        self.thisptr.remove_maximal_simplex(simplex)

    def persistence(self, homology_coeff_field=11, min_persistence=0):
        if self.pcohptr != NULL:
            del self.pcohptr
        self.pcohptr = new Rips_complex_persistence_interface(self.thisptr)
        cdef vector[pair[int, pair[double, double]]] persistence_result
        if self.pcohptr != NULL:
            persistence_result = self.pcohptr.get_persistence(homology_coeff_field, min_persistence)
        return persistence_result

    def betti_numbers(self):
        cdef vector[int] bn_result
        if self.pcohptr != NULL:
            bn_result = self.pcohptr.betti_numbers()
        else:
            print("betti_numbers function requires persistence function"
                  " to be launched first.")
        return bn_result

    def persistent_betti_numbers(self, from_value, to_value):
        cdef vector[int] pbn_result
        if self.pcohptr != NULL:
            pbn_result = self.pcohptr.persistent_betti_numbers(<double>from_value, <double>to_value)
        else:
            print("persistent_betti_numbers function requires persistence function"
                  " to be launched first.")
        return pbn_result
