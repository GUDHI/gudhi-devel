# This file is part of the Gudhi Library - https://gudhi.inria.fr/ -
# which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full
# license details.
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
from libcpp cimport bool
import errno
import os
import sys

import numpy as np
cimport numpy as np

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

# Necessary because of PyArray_SimpleNewFromData
np.import_array()

cdef extern from "Cubical_complex_interface.h" namespace "Gudhi":
    cdef cppclass Bitmap_cubical_complex_interface "Gudhi::Cubical_complex::Cubical_complex_interface":
        Bitmap_cubical_complex_interface(vector[unsigned] dimensions, vector[double] cells, bool input_top_cells) nogil except +
        Bitmap_cubical_complex_interface(const char* perseus_file) nogil except +
        int num_simplices() nogil
        int dimension() nogil
        vector[unsigned] shape() nogil
        vector[double] data

cdef extern from "Persistent_cohomology_interface.h" namespace "Gudhi":
    cdef cppclass Cubical_complex_persistence_interface "Gudhi::Persistent_cohomology_interface<Gudhi::Cubical_complex::Cubical_complex_interface>":
        Cubical_complex_persistence_interface(Bitmap_cubical_complex_interface * st, bool persistence_dim_max) nogil
        void compute_persistence(int homology_coeff_field, double min_persistence) nogil except +
        vector[pair[int, pair[double, double]]] get_persistence() nogil
        vector[vector[int]] cofaces_of_cubical_persistence_pairs() nogil
        vector[vector[int]] vertices_of_cubical_persistence_pairs() nogil
        vector[int] betti_numbers() nogil
        vector[int] persistent_betti_numbers(double from_value, double to_value) nogil
        vector[pair[double,double]] intervals_in_dimension(int dimension) nogil

# CubicalComplex python interface
cdef class CubicalComplex:
    """The CubicalComplex is an example of a structured complex useful in
    computational mathematics (specially rigorous numerics) and image
    analysis.
    """
    cdef Bitmap_cubical_complex_interface * thisptr
    cdef Cubical_complex_persistence_interface * pcohptr
    cdef bool _built_from_vertices

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, *, top_dimensional_cells=None, vertices=None, dimensions=None, perseus_file=''):
        """CubicalComplex constructor from the filtration values of either
        the top-dimensional cells, or the vertices.

        Note that in case `top_dimensional_cells` or `vertices` is passed as a flat array,
        its true shape can be given in `dimensions`, which is considered using the
        fortran ordering (column-major order). However, we recommend passing directly an
        array of the right shape.

        :param top_dimensional_cells: A multidimensional array of
            top dimensional cells filtration values.
        :type top_dimensional_cells: anything convertible to a numpy.ndarray

        Or

        :param vertices: A multidimensional array of vertices
            filtration values.
        :type vertices: anything convertible to a numpy.ndarray

        Or

        :param top_dimensional_cells: Filtration values of the top-dimensional cells.
        :type top_dimensional_cells: Iterable[float]
        :param dimensions: Shape of the array of top dimensional cells (Fortran order).
        :type dimensions: Iterable[int]

        Or

        :param vertices: Filtration values of the vertices.
        :type vertices: Iterable[float]
        :param dimensions: Shape of the array of vertices (Fortran order).
        :type dimensions: Iterable[int]

        Or

        :param perseus_file: A `Perseus-style <fileformats.html#perseus>`_ file name, giving the filtration values
            of top-dimensional cells.
        :type perseus_file: str
        """

    # The real cython constructor
    def __cinit__(self, *, top_dimensional_cells=None, vertices=None, dimensions=None, perseus_file=''):
        cdef const char* file
        self._built_from_vertices = False
        if perseus_file:
            if top_dimensional_cells is not None or vertices is not None or dimensions is not None:
                raise ValueError("The Perseus file contains all the information, do not specify anything else")
            perseus_file = perseus_file.encode('utf-8')
            file = perseus_file
            self._construct_from_file(file)
            return
        if top_dimensional_cells is not None:
            if vertices is not None:
                raise ValueError("Can only specify the top dimensional cells OR the vertices, not both")
            array = top_dimensional_cells
        elif vertices is not None:
            array = vertices
            self._built_from_vertices = True
        else:
            raise ValueError("Must specify one of top_dimensional_cells, vertices, or perseus_file")
        if dimensions is None:
            array = np.array(array, copy=False, order='F')
            dimensions = array.shape
            array = array.ravel(order='F')
        self._construct_from_cells(dimensions, array, vertices is None)

    def _construct_from_cells(self, vector[unsigned] dimensions, vector[double] cells, bool input_top_cells):
        with nogil:
            self.thisptr = new Bitmap_cubical_complex_interface(dimensions, cells, input_top_cells)

    def _construct_from_file(self, const char* filename):
        with nogil:
            self.thisptr = new Bitmap_cubical_complex_interface(filename)

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr
        if self.pcohptr != NULL:
            del self.pcohptr

    def _is_defined(self):
        """Returns true if CubicalComplex pointer is not NULL.
         """
        return self.thisptr != NULL

    def _is_persistence_defined(self):
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

    def all_cells(self):
        """Array with the filtration values of all the cells of the complex.
        Modifying the values is strongly discouraged.

        :returns:  numpy.ndarray
        """
        cdef np.npy_intp dim = self.thisptr.data.size()
        a = np.PyArray_SimpleNewFromData(1, &dim, np.NPY_DOUBLE, self.thisptr.data.data())
        np.set_array_base(a, self)
        return a.reshape([2*d+1 for d in self.thisptr.shape()], order='F')

    def top_dimensional_cells(self):
        """Array with the filtration values of the top-dimensional cells of the complex.
        Modifying the values is strongly discouraged.

        :returns:  numpy.ndarray
        """
        return self.all_cells()[(slice(1, None, 2),) * self.thisptr.dimension()]

    def vertices(self):
        """Array with the filtration values of the vertices of the complex.
        Modifying the values is strongly discouraged.

        :returns:  numpy.ndarray
        """
        return self.all_cells()[(slice(0, None, 2),) * self.thisptr.dimension()]

    def compute_persistence(self, homology_coeff_field=11, min_persistence=0):
        """This function computes the persistence of the complex, so it can be
        accessed through :func:`persistent_betti_numbers`,
        :func:`persistence_intervals_in_dimension`, etc. This function is
        equivalent to :func:`persistence` when you do not want the list
        :func:`persistence` returns.

        :param homology_coeff_field: The homology coefficient field. Must be a
            prime number. Default value is 11. Max is 46337.
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
        assert self._is_defined()
        cdef int field = homology_coeff_field
        cdef double minp = min_persistence
        with nogil:
            self.pcohptr = new Cubical_complex_persistence_interface(self.thisptr, 1)
            self.pcohptr.compute_persistence(field, minp)

    def persistence(self, homology_coeff_field=11, min_persistence=0):
        """This function computes and returns the persistence of the complex.

        :param homology_coeff_field: The homology coefficient field. Must be a
            prime number. Default value is 11. Max is 46337.
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

    def cofaces_of_persistence_pairs(self):
        """A persistence interval is described by a pair of cells, one that creates the
        feature and one that kills it. The filtration values of those 2 cells give coordinates
        for a point in a persistence diagram, or a bar in a barcode. Structurally, in the
        cubical complexes provided here, the filtration value of any cell is the minimum of the
        filtration values of the maximal cells that contain it. Connecting persistence diagram
        coordinates to the corresponding value in the input (i.e. the filtration values of
        the top-dimensional cells) is useful for differentiation purposes.

        Since the cubical complex construction from vertices is different from the top dimensional one,
        the former may not have an equivalent with the second construction, and vice versa.
        Therefore, using :func:`cofaces_of_persistence_pairs` with a cubical complex constructed from vertices
        can lead to an undefined behavior.

        This function returns a list of pairs of top-dimensional cells corresponding to
        the persistence birth and death cells of the filtration. The cells are represented by
        their indices in the input list of top-dimensional cells (and not their indices in the
        internal datastructure that includes non-maximal cells). Note that when two adjacent
        top-dimensional cells have the same filtration value, we arbitrarily return one of the two
        when calling the function on one of their common faces.

        :returns: The top-dimensional cells/cofaces of the positive and negative cells,
            together with the corresponding homological dimension, in two lists of numpy arrays of integers.
            The first list contains the regular persistence pairs, grouped by dimension.
            It contains numpy arrays of shape (number_of_persistence_points, 2).
            The indices of the arrays in the list correspond to the homological dimensions, and the
            integers of each row in each array correspond to: (index of positive top-dimensional cell,
            index of negative top-dimensional cell).
            The second list contains the essential features, grouped by dimension.
            It contains numpy arrays of shape (number_of_persistence_points,).
            The indices of the arrays in the list correspond to the homological dimensions, and the
            integers of each row in each array correspond to: (index of positive top-dimensional cell).
        """

        assert self.pcohptr != NULL, "compute_persistence() must be called before cofaces_of_persistence_pairs()"
        assert not self._built_from_vertices, (
                "cofaces_of_persistence_pairs() only makes sense for a complex"
                " initialized from the values of the top-dimensional cells"
        )

        cdef vector[vector[int]] persistence_result
        output = [[],[]]
        with nogil:
            persistence_result = self.pcohptr.cofaces_of_cubical_persistence_pairs()
        pr = np.array(persistence_result)

        ess_ind = np.argwhere(pr[:,2] == -1)[:,0]
        ess = pr[ess_ind]
        max_h = np.max(ess[:,0])+1 if len(ess) > 0 else 0
        for h in range(max_h):
            hidxs = np.argwhere(ess[:,0] == h)[:,0]
            output[1].append(ess[hidxs][:,1])

        reg_ind = np.setdiff1d(np.arange(len(pr)), ess_ind)
        reg = pr[reg_ind]
        max_h = np.max(reg[:,0])+1 if len(reg) > 0 else 0
        for h in range(max_h):
            hidxs = np.argwhere(reg[:,0] == h)[:,0]
            output[0].append(reg[hidxs][:,1:])

        return output

    def vertices_of_persistence_pairs(self):
        """This function returns a list of pairs of vertices corresponding to
        the persistence birth and death cells of the filtration. The cells are represented by
        their indices in the input list of vertices (and not their indices in the
        internal datastructure that includes non-minimal cells). Note that when two adjacent
        vertices have the same filtration value, we arbitrarily return one of the two
        when calling the function on one of their common faces.

        Since the cubical complex construction from vertices is different from the top dimensional one,
        the former may not have an equivalent with the second construction, and vice versa.
        Therefore, using :func:`vertices_of_persistence_pairs` with a cubical complex constructed from
        top_dimensional_cells can lead to an undefined behavior.

        :returns: The vertices of the positive and negative cells,
            together with the corresponding homological dimension, in two lists of numpy arrays of integers.
            The first list contains the regular persistence pairs, grouped by dimension.
            It contains numpy arrays of shape (number_of_persistence_points, 2).
            The indices of the arrays in the list correspond to the homological dimensions, and the
            integers of each row in each array correspond to: (index of positive vertex,
            index of negative vertex).
            The second list contains the essential features, grouped by dimension.
            It contains numpy arrays of shape (number_of_persistence_points,).
            The indices of the arrays in the list correspond to the homological dimensions, and the
            integers of each row in each array correspond to: (index of positive vertex).
        """

        assert self.pcohptr != NULL, "compute_persistence() must be called before vertices_of_persistence_pairs()"
        assert self._built_from_vertices, (
                "vertices_of_persistence_pairs() only makes sense for a complex"
                " initialized from the values of the vertices"
        )

        cdef vector[vector[int]] persistence_result
        output = [[],[]]
        with nogil:
            persistence_result = self.pcohptr.vertices_of_cubical_persistence_pairs()
        pr = np.array(persistence_result)

        ess_ind = np.argwhere(pr[:,2] == -1)[:,0]
        ess = pr[ess_ind]
        max_h = max(ess[:,0])+1 if len(ess) > 0 else 0
        for h in range(max_h):
            hidxs = np.argwhere(ess[:,0] == h)[:,0]
            output[1].append(ess[hidxs][:,1])

        reg_ind = np.setdiff1d(np.array(range(len(pr))), ess_ind)
        reg = pr[reg_ind]
        max_h = max(reg[:,0])+1 if len(reg) > 0 else 0
        for h in range(max_h):
            hidxs = np.argwhere(reg[:,0] == h)[:,0]
            output[0].append(reg[hidxs][:,1:])

        return output

    def betti_numbers(self):
        """This function returns the Betti numbers of the complex.

        :returns: list of int -- The Betti numbers ([B0, B1, ..., Bn]).

        :note: betti_numbers function requires :func:`compute_persistence` function to be
            launched first.

        :note: betti_numbers function always returns [1, 0, 0, ...] as infinity
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
        piid = np.array(self.pcohptr.intervals_in_dimension(dimension))
        # Workaround https://github.com/GUDHI/gudhi-devel/issues/507
        if len(piid) == 0:
            return np.empty(shape = [0, 2])
        return piid
