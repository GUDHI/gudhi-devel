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

cdef extern from "Cubical_complex_interface.h" namespace "Gudhi":
    cdef cppclass Bitmap_cubical_complex_base_interface "Gudhi::Cubical_complex::Cubical_complex_interface<>":
        Bitmap_cubical_complex_base_interface(vector[unsigned] dimensions, vector[double] top_dimensional_cells)
        Bitmap_cubical_complex_base_interface(string perseus_file)
        int num_simplices()
        int dimension()

cdef extern from "Persistent_cohomology_interface.h" namespace "Gudhi":
    cdef cppclass Cubical_complex_persistence_interface "Gudhi::Persistent_cohomology_interface<Gudhi::Cubical_complex::Cubical_complex_interface<>>":
        Cubical_complex_persistence_interface(Bitmap_cubical_complex_base_interface * st, bool persistence_dim_max)
        vector[pair[int, pair[double, double]]] get_persistence(int homology_coeff_field, double min_persistence)
        vector[int] betti_numbers()
        vector[int] persistent_betti_numbers(double from_value, double to_value)
        vector[pair[double,double]] intervals_in_dimension(int dimension)

# CubicalComplex python interface
cdef class CubicalComplex:
    """The CubicalComplex is an example of a structured complex useful in
    computational mathematics (specially rigorous numerics) and image
    analysis.
    """
    cdef Bitmap_cubical_complex_base_interface * thisptr

    cdef Cubical_complex_persistence_interface * pcohptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, dimensions=None, top_dimensional_cells=None,
                  perseus_file=''):
        """CubicalComplex constructor from dimensions and
        top_dimensional_cells or from a perseus file style name.

        :param dimensions: A list of number of top dimensional cells.
        :type dimensions: list of int
        :param top_dimensional_cells: A list of top dimensional cells.
        :type top_dimensional_cells: list of double

        Or

        :param perseus_file: A perseus file style name.
        :type perseus_file: string
        """

    # The real cython constructor
    def __cinit__(self, dimensions=None, top_dimensional_cells=None,
                  perseus_file=''):
        if (dimensions is not None) and (top_dimensional_cells is not None) and (perseus_file is ''):
            self.thisptr = new Bitmap_cubical_complex_base_interface(dimensions, top_dimensional_cells)
        elif (dimensions is None) and (top_dimensional_cells is None) and (perseus_file is not ''):
            if os.path.isfile(perseus_file):
                self.thisptr = new Bitmap_cubical_complex_base_interface(perseus_file)
            else:
                print("file " + perseus_file + " not found.")
        else:
            print("CubicalComplex can be constructed from dimensions and "
              "top_dimensional_cells or from a perseus file style name.")

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr
        if self.pcohptr != NULL:
            del self.pcohptr

    def __is_defined(self):
        """Returns true if CubicalComplex pointer is not NULL.
         """
        return self.thisptr != NULL

    def __is_persistence_defined(self):
        """Returns true if Persistence pointer is not NULL.
         """
        return self.pcohptr != NULL

    def num_simplices(self):
        """This function returns the number of simplices of the simplicial
        complex.

        :returns:  int -- the simplicial complex number of simplices.
        """
        return self.thisptr.num_simplices()

    def dimension(self):
        """This function returns the dimension of the simplicial complex.

        :returns:  int -- the simplicial complex dimension.
        """
        return self.thisptr.dimension()

    def persistence(self, homology_coeff_field=11, min_persistence=0):
        """This function returns the persistence of the simplicial complex.

        :param homology_coeff_field: The homology coefficient field. Must be a
            prime number
        :type homology_coeff_field: int.
        :param min_persistence: The minimum persistence value to take into
            account (strictly greater than min_persistence). Default value is
            0.0.
            Sets min_persistence to -1.0 to see all values.
        :type min_persistence: float.
        :returns: list of pairs(dimension, pair(birth, death)) -- the
            persistence of the simplicial complex.
        """
        if self.pcohptr != NULL:
            del self.pcohptr
        if self.thisptr != NULL:
            self.pcohptr = new Cubical_complex_persistence_interface(self.thisptr, True)
        cdef vector[pair[int, pair[double, double]]] persistence_result
        if self.pcohptr != NULL:
            persistence_result = self.pcohptr.get_persistence(homology_coeff_field, min_persistence)
        return persistence_result

    def betti_numbers(self):
        """This function returns the Betti numbers of the simplicial complex.

        :returns: list of int -- The Betti numbers ([B0, B1, ..., Bn]).

        :note: betti_numbers function requires persistence function to be
            launched first.
        """
        cdef vector[int] bn_result
        if self.pcohptr != NULL:
            bn_result = self.pcohptr.betti_numbers()
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

        :returns: list of int -- The persistent Betti numbers ([B0, B1, ...,
            Bn]).

        :note: persistent_betti_numbers function requires persistence
            function to be launched first.
        """
        cdef vector[int] pbn_result
        if self.pcohptr != NULL:
            pbn_result = self.pcohptr.persistent_betti_numbers(<double>from_value, <double>to_value)
        return pbn_result

    def intervals_in_dim(self, dimension):
        """This function returns the persistence intervals of the simplicial
        complex in a specific dimension.

        :param dimension: The specific dimension.
        :type from_value: int.
        :returns: The persistence intervals.
        :rtype:  list of pair of float

        :note: betti_numbers function requires persistence function to be
            launched first.
        """
        cdef vector[pair[double,double]] intervals_result
        if self.pcohptr != NULL:
            intervals_result = self.pcohptr.intervals_in_dimension(dimension)
        else:
            print("intervals_in_dim function requires persistence function"
                  " to be launched first.")
        return intervals_result
