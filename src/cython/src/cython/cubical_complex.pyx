from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair

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
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016  INRIA Saclay (France)"
__license__ = "GPL v3"

cdef extern from "Cubical_complex_interface.h" namespace "Gudhi":
    cdef cppclass Bitmap_cubical_complex_base_interface "Gudhi::Cubical_complex::Cubical_complex_interface<>":
          Bitmap_cubical_complex_base_interface(vector[unsigned] dimensions, vector[double] top_dimensional_cells)

cdef extern from "Persistent_cohomology_interface.h" namespace "Gudhi":
    cdef cppclass Cubical_complex_persistence_interface "Gudhi::Persistent_cohomology_interface<Gudhi::Cubical_complex::Cubical_complex_interface<>>":
        Cubical_complex_persistence_interface(Bitmap_cubical_complex_base_interface * st)
        vector[pair[int, pair[double, double]]] get_persistence(int homology_coeff_field, double min_persistence)
        vector[int] betti_numbers()
        vector[int] persistent_betti_numbers(double from_value, double to_value)


# CubicalComplex python interface
cdef class CubicalComplex:
    cdef Bitmap_cubical_complex_base_interface * thisptr

    cdef Cubical_complex_persistence_interface * pcohptr

    def __cinit__(self, dimensions=None, top_dimensional_cells=None):
        if dimensions is not None:
            if top_dimensional_cells is not None:
                self.thisptr = new Bitmap_cubical_complex_base_interface(dimensions, top_dimensional_cells)

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr
        if self.pcohptr != NULL:
            del self.pcohptr

    def persistence(self, homology_coeff_field=11, min_persistence=0):
        if self.pcohptr != NULL:
            del self.pcohptr
        self.pcohptr = new Cubical_complex_persistence_interface(self.thisptr)
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
