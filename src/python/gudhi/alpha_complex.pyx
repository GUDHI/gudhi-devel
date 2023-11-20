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
from libcpp.string cimport string
from libcpp cimport bool
from libc.stdint cimport intptr_t
import warnings

from gudhi.simplex_tree cimport *
from gudhi.simplex_tree import SimplexTree
from gudhi import read_points_from_off_file

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "GPL v3"

cdef extern from "Alpha_complex_interface.h" namespace "Gudhi":
    cdef cppclass Alpha_complex_interface "Gudhi::alpha_complex::Alpha_complex_interface":
        Alpha_complex_interface(vector[vector[double]] points, vector[double] weights, bool fast_version, bool exact_version) nogil except +
        vector[double] get_point(int vertex) nogil except +
        void create_simplex_tree(Simplex_tree_python_interface* simplex_tree, double max_alpha_square, bool default_filtration_value) nogil except +
        @staticmethod
        void set_float_relative_precision(double precision) nogil
        @staticmethod
        double get_float_relative_precision() nogil

# AlphaComplex python interface
cdef class AlphaComplex:
    """AlphaComplex is a simplicial complex constructed from the finite cells of a Delaunay Triangulation.

    The filtration value of each simplex is computed as the square of the circumradius of the simplex if the
    circumsphere is empty (the simplex is then said to be Gabriel), and as the minimum of the filtration values of the
    codimension 1 cofaces that make it not Gabriel otherwise.

    All simplices that have a filtration value strictly greater than a given alpha squared value are not inserted into
    the complex.

    .. note::

        When Alpha_complex is constructed with an infinite value of alpha, the complex is a Delaunay complex.
    """

    cdef Alpha_complex_interface * this_ptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, points=[], off_file='', weights=None, precision='safe'):
        """AlphaComplex constructor.

        :param points: A list of points in d-Dimension.
        :type points: Iterable[Iterable[float]]

        :param off_file: **[deprecated]** An `OFF file style <fileformats.html#off-file-format>`_ name.
            If an `off_file` is given with `points` as arguments, only points from the file are taken into account.
        :type off_file: string

        :param weights: A list of weights. If set, the number of weights must correspond to the number of points.
        :type weights: Iterable[float]

        :param precision: Alpha complex precision can be 'fast', 'safe' or 'exact'. Default is 'safe'.
        :type precision: string

        :raises FileNotFoundError: **[deprecated]** If `off_file` is set but not found.
        :raises ValueError: In case of inconsistency between the number of points and weights.
        """

    # The real cython constructor
    def __cinit__(self, points = [], off_file = '', weights=None, precision = 'safe'):
        assert precision in ['fast', 'safe', 'exact'], "Alpha complex precision can only be 'fast', 'safe' or 'exact'"
        cdef bool fast = precision == 'fast'
        cdef bool exact = precision == 'exact'

        if off_file:
            warnings.warn("off_file is a deprecated parameter, please consider using gudhi.read_points_from_off_file",
                          DeprecationWarning)
            points = read_points_from_off_file(off_file = off_file)

        # weights are set but is inconsistent with the number of points
        if weights is not None and len(weights) != len(points):
            raise ValueError("Inconsistency between the number of points and weights")

        # need to copy the points to use them without the gil
        cdef vector[vector[double]] pts
        cdef vector[double] wgts
        pts = points
        if weights is not None:
            wgts = weights
        with nogil:
            self.this_ptr = new Alpha_complex_interface(pts, wgts, fast, exact)

    def __dealloc__(self):
        if self.this_ptr != NULL:
            del self.this_ptr

    def _is_defined(self):
        """Returns true if AlphaComplex pointer is not NULL.
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
                                              mas, compute_filtration)
        return stree

    @staticmethod
    def set_float_relative_precision(precision):
        """
        :param precision: When the AlphaComplex is constructed with :code:`precision = 'safe'` (the default),
            one can set the float relative precision of filtration values computed in
            :func:`~gudhi.AlphaComplex.create_simplex_tree`.
            Default is :code:`1e-5` (cf. :func:`~gudhi.AlphaComplex.get_float_relative_precision`).
            For more details, please refer to
            `CGAL::Lazy_exact_nt<NT>::set_relative_precision_of_to_double <https://doc.cgal.org/latest/Number_types/classCGAL_1_1Lazy__exact__nt.html>`_
        :type precision: float
        """
        if precision <=0. or precision >= 1.:
            raise ValueError("Relative precision value must be strictly greater than 0 and strictly lower than 1")
        Alpha_complex_interface.set_float_relative_precision(precision)
    
    @staticmethod
    def get_float_relative_precision():
        """
        :returns: The float relative precision of filtration values computation in
            :func:`~gudhi.AlphaComplex.create_simplex_tree` when the AlphaComplex is constructed with
            :code:`precision = 'safe'` (the default).
        :rtype: float
        """
        return Alpha_complex_interface.get_float_relative_precision()
