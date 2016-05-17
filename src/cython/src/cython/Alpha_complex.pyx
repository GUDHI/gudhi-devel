"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Vincent Rouvreau

   Copyright (C) 2016  INRIA Saclay (France)

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>."""

__author__      = "Vincent Rouvreau"
__copyright__   = "Copyright (C) 2016  INRIA Saclay (France)"
__license__     = "GPL v3"

from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair

cdef extern from "Alpha_complex_interface.h" namespace "Gudhi":
    cdef cppclass Alpha_complex_interface "Gudhi::alphacomplex::Alpha_complex_interface":
        Alpha_complex_interface(vector[vector[double]] points,double max_alpha_square)
        double filtration()
        double simplex_filtration(vector[int] simplex)
        void set_filtration(double filtration)
        void initialize_filtration()
        int num_vertices()
        int num_simplices()
        void set_dimension(int dimension)
        int dimension()
        bint find_simplex(vector[int] simplex)
        bint insert_simplex_and_subfaces(vector[int] simplex, double filtration)
        vector[pair[vector[int], double]] get_filtered_tree()
        vector[pair[vector[int], double]] get_skeleton_tree(int dimension)
        vector[pair[vector[int], double]] get_star_tree(vector[int] simplex)
        vector[pair[vector[int], double]] get_coface_tree(vector[int] simplex, int dimension)
        void remove_maximal_simplex(vector[int] simplex)
        vector[double] get_point(int vertex)

# AlphaComplex python interface
cdef class AlphaComplex:
    cdef Alpha_complex_interface *thisptr
    def __cinit__(self, points=None, max_alpha_square=float('inf')):
        if points is not None:
            self.thisptr = new Alpha_complex_interface(points, max_alpha_square)
    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr
    def get_filtration(self):
        return self.thisptr.filtration()
    def filtration(self, simplex):
        return self.thisptr.simplex_filtration(simplex)
    def set_filtration(self, filtration):
        self.thisptr.set_filtration(<double>filtration)
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
    def insert(self, simplex, filtration = 0.0):
        cdef vector[int] complex
        for i in simplex:
          complex.push_back(i)
        return self.thisptr.insert_simplex_and_subfaces(complex, <double>filtration)
    def get_filtered_tree(self):
        cdef vector[pair[vector[int], double]] coface_tree = self.thisptr.get_filtered_tree()
        ct = []
        for filtered_complex in coface_tree:
            v = []
            for vertex in filtered_complex.first:
                v.append(vertex)
            ct.append((v,filtered_complex.second))
        return ct
    def get_skeleton_tree(self, dim):
        cdef vector[pair[vector[int], double]] coface_tree = self.thisptr.get_skeleton_tree(<int>dim)
        ct = []
        for filtered_complex in coface_tree:
            v = []
            for vertex in filtered_complex.first:
                v.append(vertex)
            ct.append((v,filtered_complex.second))
        return ct
    def get_star_tree(self, simplex):
        cdef vector[int] complex
        for i in simplex:
          complex.push_back(i)
        cdef vector[pair[vector[int], double]] coface_tree = self.thisptr.get_star_tree(complex)
        ct = []
        for filtered_complex in coface_tree:
            v = []
            for vertex in filtered_complex.first:
                v.append(vertex)
            ct.append((v,filtered_complex.second))
        return ct
    def get_coface_tree(self, simplex, dim):
        cdef vector[int] complex
        for i in simplex:
          complex.push_back(i)
        cdef vector[pair[vector[int], double]] coface_tree = self.thisptr.get_coface_tree(complex, <int>dim)
        ct = []
        for filtered_complex in coface_tree:
            v = []
            for vertex in filtered_complex.first:
                v.append(vertex)
            ct.append((v,filtered_complex.second))
        return ct
    def remove_maximal_simplex(self, simplex):
        self.thisptr.remove_maximal_simplex(simplex)
    def get_point(self, vertex):
        cdef vector[double] point = self.thisptr.get_point(vertex)
        return point
