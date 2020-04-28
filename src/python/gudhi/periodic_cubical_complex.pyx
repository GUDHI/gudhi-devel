# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from __future__ import print_function
from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp cimport bool
import sys
import os

import numpy as np

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

cdef extern from "Cubical_complex_interface.h" namespace "Gudhi":
    cdef cppclass Periodic_cubical_complex_base_interface "Gudhi::Cubical_complex::Cubical_complex_interface<Gudhi::cubical_complex::Bitmap_cubical_complex_periodic_boundary_conditions_base<double>>":
        Periodic_cubical_complex_base_interface(vector[unsigned] dimensions, vector[double] top_dimensional_cells, vector[bool] periodic_dimensions)
        Periodic_cubical_complex_base_interface(string perseus_file)
        int num_simplices()
        int dimension()

cdef extern from "Persistent_cohomology_interface.h" namespace "Gudhi":
    cdef cppclass Periodic_cubical_complex_persistence_interface "Gudhi::Persistent_cohomology_interface<Gudhi::Cubical_complex::Cubical_complex_interface<Gudhi::cubical_complex::Bitmap_cubical_complex_periodic_boundary_conditions_base<double>>>":
        Periodic_cubical_complex_persistence_interface(Periodic_cubical_complex_base_interface * st, bool persistence_dim_max)
        void compute_persistence(int homology_coeff_field, double min_persistence)
        vector[pair[int, pair[double, double]]] get_persistence()
        vector[int] betti_numbers()
        vector[int] persistent_betti_numbers(double from_value, double to_value)
        vector[pair[double,double]] intervals_in_dimension(int dimension)

# PeriodicCubicalComplex python interface
cdef class PeriodicCubicalComplex:
    """The PeriodicCubicalComplex is an example of a structured complex useful
    in computational mathematics (specially rigorous numerics) and image
    analysis.
    """
    cdef Periodic_cubical_complex_base_interface * thisptr

    cdef Periodic_cubical_complex_persistence_interface * pcohptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, dimensions=None, top_dimensional_cells=None,
                 periodic_dimensions=None, perseus_file=''):
        """PeriodicCubicalComplex constructor from dimensions and
        top_dimensional_cells or from a Perseus-style file name.

        :param dimensions: A list of number of top dimensional cells.
        :type dimensions: list of int
        :param top_dimensional_cells: A list of cells filtration values.
        :type top_dimensional_cells: list of double
        :param periodic_dimensions: A list of top dimensional cells periodicity value.
        :type periodic_dimensions: list of boolean

        Or

        :param top_dimensional_cells: A multidimensional array of cells
            filtration values.
        :type top_dimensional_cells: anything convertible to a numpy ndarray
        :param periodic_dimensions: A list of top dimensional cells periodicity value.
        :type periodic_dimensions: list of boolean

        Or

        :param perseus_file: A Perseus-style file name.
        :type perseus_file: string
        """

    # The real cython constructor
    def __cinit__(self, dimensions=None, top_dimensional_cells=None,
                  periodic_dimensions=None, perseus_file=''):
        if ((dimensions is not None) and (top_dimensional_cells is not None)
            and (periodic_dimensions is not None) and (perseus_file == '')):
            self.thisptr = new Periodic_cubical_complex_base_interface(dimensions,
                                                                       top_dimensional_cells,
                                                                       periodic_dimensions)
        elif ((dimensions is None) and (top_dimensional_cells is not None)
            and (periodic_dimensions is not None) and (perseus_file == '')):
            top_dimensional_cells = np.array(top_dimensional_cells,
                                             copy = False,
                                             order = 'F')
            dimensions = top_dimensional_cells.shape
            top_dimensional_cells = top_dimensional_cells.ravel(order='F')
            self.thisptr = new Periodic_cubical_complex_base_interface(dimensions,
                                                                       top_dimensional_cells,
                                                                       periodic_dimensions)
        elif ((dimensions is None) and (top_dimensional_cells is None)
            and (periodic_dimensions is None) and (perseus_file != '')):
            if os.path.isfile(perseus_file):
                self.thisptr = new Periodic_cubical_complex_base_interface(perseus_file.encode('utf-8'))
            else:
                print("file " + perseus_file + " not found.", file=sys.stderr)
        else:
            print("CubicalComplex can be constructed from dimensions, "
              "top_dimensional_cells and periodic_dimensions, or from "
              "top_dimensional_cells and periodic_dimensions or from "
              "a Perseus-style file name.", file=sys.stderr)

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr
        if self.pcohptr != NULL:
            del self.pcohptr

    def __is_defined(self):
        """Returns true if PeriodicCubicalComplex pointer is not NULL.
         """
        return self.thisptr != NULL

    def __is_persistence_defined(self):
        """Returns true if Persistence pointer is not NULL.
         """
        return self.pcohptr != NULL

    def num_simplices(self):
        """This function returns the number of all cubes in the complex.

        :returns:  int -- the number of all cubes in the complex.
        """
        return self.thisptr.num_simplices()

    def dimension(self):
        """This function returns the dimension of the complex.

        :returns:  int -- the complex dimension.
        """
        return self.thisptr.dimension()

    def compute_persistence(self, homology_coeff_field=11, min_persistence=0):
        """This function computes the persistence of the complex, so it can be
        accessed through :func:`persistent_betti_numbers`,
        :func:`persistence_intervals_in_dimension`, etc. This function is
        equivalent to :func:`persistence` when you do not want the list
        :func:`persistence` returns.

        :param homology_coeff_field: The homology coefficient field. Must be a
            prime number
        :type homology_coeff_field: int.
        :param min_persistence: The minimum persistence value to take into
            account (strictly greater than min_persistence). Default value is
            0.0.
            Sets min_persistence to -1.0 to see all values.
        :type min_persistence: float.
        :returns: Nothing.
        """
        if self.pcohptr != NULL:
            del self.pcohptr
        assert self.__is_defined()
        self.pcohptr = new Periodic_cubical_complex_persistence_interface(self.thisptr, True)
        self.pcohptr.compute_persistence(homology_coeff_field, min_persistence)

    def persistence(self, homology_coeff_field=11, min_persistence=0):
        """This function computes and returns the persistence of the complex.

        :param homology_coeff_field: The homology coefficient field. Must be a
            prime number
        :type homology_coeff_field: int.
        :param min_persistence: The minimum persistence value to take into
            account (strictly greater than min_persistence). Default value is
            0.0.
            Sets min_persistence to -1.0 to see all values.
        :type min_persistence: float.
        :returns: list of pairs(dimension, pair(birth, death)) -- the
            persistence of the complex.
        """
        self.compute_persistence(homology_coeff_field, min_persistence)
        return self.pcohptr.get_persistence()

    def betti_numbers(self):
        """This function returns the Betti numbers of the complex.

        :returns: list of int -- The Betti numbers ([B0, B1, ..., Bn]).

        :note: betti_numbers function requires :func:`compute_persistence` function to be
            launched first.

        :note: This function always returns the Betti numbers of a torus as infinity
            filtration cubes are not removed from the complex.
        """
        assert self.pcohptr != NULL, "compute_persistence() must be called before betti_numbers()"
        return self.pcohptr.betti_numbers()

    def persistent_betti_numbers(self, from_value, to_value):
        """This function returns the persistent Betti numbers of the complex.

        :param from_value: The persistence birth limit to be added in the
            numbers (persistent birth <= from_value).
        :type from_value: float.
        :param to_value: The persistence death limit to be added in the
            numbers (persistent death > to_value).
        :type to_value: float.

        :returns: list of int -- The persistent Betti numbers ([B0, B1, ...,
            Bn]).

        :note: persistent_betti_numbers function requires :func:`compute_persistence`
            function to be launched first.
        """
        assert self.pcohptr != NULL, "compute_persistence() must be called before persistent_betti_numbers()"
        return self.pcohptr.persistent_betti_numbers(<double>from_value, <double>to_value)

    def persistence_intervals_in_dimension(self, dimension):
        """This function returns the persistence intervals of the complex in a
        specific dimension.

        :param dimension: The specific dimension.
        :type dimension: int.
        :returns: The persistence intervals.
        :rtype:  numpy array of dimension 2

        :note: intervals_in_dim function requires :func:`compute_persistence` function to be
            launched first.
        """
        assert self.pcohptr != NULL, "compute_persistence() must be called before persistence_intervals_in_dimension()"
        return np.array(self.pcohptr.intervals_in_dimension(dimension))
