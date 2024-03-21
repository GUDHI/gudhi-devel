# This file is part of the Gudhi Library - https://gudhi.inria.fr/ -
# which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full
# license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2024 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from __future__ import print_function
from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp cimport bool
from libc.stdint cimport intptr_t
import warnings

from gudhi.simplex_tree cimport *
from gudhi.simplex_tree import SimplexTree
from gudhi import read_points_from_off_file

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2024 Inria"
__license__ = "GPL v3"

cdef extern from "Delaunay_complex_interface.h" namespace "Gudhi":
    cdef cppclass Delaunay_complex_interface "Gudhi::delaunay_complex::Delaunay_complex_interface":
        Delaunay_complex_interface(vector[vector[double]] points, vector[double] weights, bool fast_version, bool exact_version) nogil except +
        vector[double] get_point(int vertex) nogil except +
        void create_simplex_tree(Simplex_tree_python_interface* simplex_tree, double max_alpha_square, bool default_filtration_value, bool assign_meb_filtration) nogil except +
        @staticmethod
        void set_float_relative_precision(double precision) nogil
        @staticmethod
        double get_float_relative_precision() nogil

# DelaunayCechComplex python interface
cdef class DelaunayCechComplex:
    """DelaunayCechComplex is a simplicial complex constructed from the finite cells of a Delaunay Triangulation.

    The filtration value of each simplex is computed as the square of the circumradius of the simplex if the
    circumsphere is empty (the simplex is then said to be Gabriel), and as the minimum of the filtration values of the
    codimension 1 cofaces that make it not Gabriel otherwise.

    All simplices that have a filtration value strictly greater than a given alpha squared value are not inserted into
    the complex.

    .. note::

        When DelaunayCechComplex is constructed with an infinite value of alpha, the complex is a Delaunay complex.
    """

    cdef Delaunay_complex_interface * this_ptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, points=[], precision='safe'):
        """DelaunayCechComplex constructor.

        :param points: A list of points in d-Dimension.
        :type points: Iterable[Iterable[float]]

        :param precision: Delaunay Cech complex precision can be 'fast', 'safe' or 'exact'. Default is 'safe'.
        :type precision: string
        """

    # The real cython constructor
    def __cinit__(self, points = [], precision = 'safe'):
        assert precision in ['fast', 'safe', 'exact'], "Delaunay Cech complex precision can only be 'fast', 'safe' or 'exact'"
        cdef bool fast = precision == 'fast'
        cdef bool exact = precision == 'exact'

        # need to copy the points to use them without the gil
        cdef vector[vector[double]] pts
        # use empty weights on purpose - Weighted DelaunayCechComplex is not available
        cdef vector[double] wgts
        pts = points
        with nogil:
            self.this_ptr = new Delaunay_complex_interface(pts, wgts, fast, exact)

    def __dealloc__(self):
        if self.this_ptr != NULL:
            del self.this_ptr

    def _is_defined(self):
        """Returns true if DelaunayCechComplex pointer is not NULL.
         """
        return self.this_ptr != NULL

    def get_point(self, vertex):
        """This function returns the point corresponding to a given vertex from the :class:`~gudhi.SimplexTree`.

        :param vertex: The vertex.
        :type vertex: int
        :rtype: list of float
        :returns: the point.
        """
        return self.this_ptr.get_point(vertex)

    def create_simplex_tree(self, max_alpha_square = float('inf'), default_filtration_value = False):
        """
        :param max_alpha_square: The maximum alpha square threshold the simplices shall not exceed. Default is set to
            infinity, and there is very little point using anything else since it does not save time.
        :type max_alpha_square: float
        :param default_filtration_value: Set this value to `True` if filtration values are not needed to be computed
            (will be set to `NaN`). Default value is `False` (which means compute the filtration values).
        :type default_filtration_value: bool
        :returns: A simplex tree created from the Delaunay Triangulation.
        :rtype: SimplexTree
        """
        stree = SimplexTree()
        cdef double mas = max_alpha_square
        cdef intptr_t stree_int_ptr=stree.thisptr
        cdef bool compute_filtration = default_filtration_value == True
        with nogil:
            self.this_ptr.create_simplex_tree(<Simplex_tree_python_interface*>stree_int_ptr,
                                              mas, compute_filtration, True)
        return stree

    @staticmethod
    def set_float_relative_precision(precision):
        """
        :param precision: When the DelaunayCechComplex is constructed with :code:`precision = 'safe'` (the default),
            one can set the float relative precision of filtration values computed in
            :func:`~gudhi.DelaunayCechComplex.create_simplex_tree`.
            Default is :code:`1e-5` (cf. :func:`~gudhi.DelaunayCechComplex.get_float_relative_precision`).
            For more details, please refer to
            `CGAL::Lazy_exact_nt<NT>::set_relative_precision_of_to_double <https://doc.cgal.org/latest/Number_types/classCGAL_1_1Lazy__exact__nt.html>`_
        :type precision: float
        """
        if precision <=0. or precision >= 1.:
            raise ValueError("Relative precision value must be strictly greater than 0 and strictly lower than 1")
        Delaunay_complex_interface.set_float_relative_precision(precision)
    
    @staticmethod
    def get_float_relative_precision():
        """
        :returns: The float relative precision of filtration values computation in
            :func:`~gudhi.DelaunayCechComplex.create_simplex_tree` when the DelaunayCechComplex is constructed with
            :code:`precision = 'safe'` (the default).
        :rtype: float
        """
        return Delaunay_complex_interface.get_float_relative_precision()
