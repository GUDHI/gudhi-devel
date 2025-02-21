# This file is part of the Gudhi Library - https://gudhi.inria.fr/ -
# which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full
# license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - 2024/03 Vincent Rouvreau: Renamed AlphaComplex as DelaunayComplex. AlphaComplex inherits from it.
#   - 2024/10 Vincent Rouvreau: Add square root filtration values interface and new function interfaces (instead of
#                               classes)
#   - YYYY/MM Author: Description of the modification

from __future__ import print_function
from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp cimport bool
from libc.stdint cimport intptr_t
from typing import Literal, Optional
from collections.abc import Iterable

from gudhi.simplex_tree cimport *
from gudhi.simplex_tree import SimplexTree

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "GPL v3"

cdef extern from "Delaunay_complex_interface.h" namespace "Gudhi":
    cdef cppclass Delaunay_filtration "Gudhi::delaunay_complex::Delaunay_filtration":
        pass

cdef extern from "Delaunay_complex_interface.h" namespace "Gudhi::delaunay_complex::Delaunay_filtration":
    cdef Delaunay_filtration NONE
    cdef Delaunay_filtration CECH
    cdef Delaunay_filtration ALPHA

cdef extern from "Delaunay_complex_interface.h" namespace "Gudhi":
    cdef cppclass Delaunay_complex_interface "Gudhi::delaunay_complex::Delaunay_complex_interface":
        Delaunay_complex_interface(vector[vector[double]] points, vector[double] weights, bool fast_version, bool exact_version) nogil except +
        vector[double] get_point(int vertex) nogil except +
        void create_simplex_tree(Simplex_tree_python_interface* simplex_tree, double max_alpha_square, Delaunay_filtration filtration, bool output_squared_values) nogil except +
        @staticmethod
        void set_float_relative_precision(double precision) nogil
        @staticmethod
        double get_float_relative_precision() nogil

# DelaunayComplex python interface
cdef class DelaunayComplex:
    """DelaunayComplex is a simplicial complex constructed from the finite cells of a Delaunay Triangulation.

    When :paramref:`~gudhi.DelaunayComplex.create_simplex_tree.filtration` is:

    * `None` (default value) - The filtration value of each simplex is not computed (set to `NaN`)
    * `'alpha'`              - The filtration value of each simplex is computed as an :class:`~gudhi.AlphaComplex`
    * `'cech'`               - The filtration value of each simplex is computed as a :class:`~gudhi.DelaunayCechComplex`
    """

    cdef Delaunay_complex_interface * this_ptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, points: Iterable[Iterable[float]] = [], weights: Optional[Iterable[float]] = None,
                 precision: Literal['fast', 'safe', 'exact'] = 'safe'):
        """Args:
            points (Iterable[Iterable[float]]): A list of points in d-Dimension.
            weights (Optional[Iterable[float]]): A list of weights. If set, the number of weights must correspond to
                the number of points.
            precision (str): Complex precision can be `'fast'`, `'safe'` or `'exact'`. Default is `'safe'`.

        :raises ValueError: In case of inconsistency between the number of points and weights.
        """

    # The real cython constructor
    def __cinit__(self, points: Iterable[Iterable[float]] = [], weights: Optional[Iterable[float]] = None,
                  precision: Literal['fast', 'safe', 'exact'] = 'safe'):
        assert precision in [
            "fast",
            "safe",
            "exact",
        ], "Delaunay complex precision can only be 'fast', 'safe' or 'exact'"
        cdef bool fast = precision == 'fast'
        cdef bool exact = precision == 'exact'

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
            self.this_ptr = new Delaunay_complex_interface(pts, wgts, fast, exact)

    def __dealloc__(self):
        if self.this_ptr != NULL:
            del self.this_ptr

    def _is_defined(self):
        """Returns `True` if DelaunayComplex pointer is not NULL.
         """
        return self.this_ptr != NULL

    def create_simplex_tree(self, double max_alpha_square: float = float('inf'),
                            filtration: Optional[Literal['alpha', 'cech']] = None,
                            bool output_squared_values: bool = True) -> SimplexTree:
        """
        Args:
            max_alpha_square: The maximum alpha square threshold the simplices shall not exceed. Default is set to
                infinity, and there is very little point using anything else since it does not save time.
            filtration: Set this value to `None` (default value) if filtration values are not needed to be computed
                (will be set to `NaN`). Set it to `alpha` to compute the filtration values with the Alpha complex, or
                to `cech` to compute the Delaunay Cech complex.
            output_squared_values: Square filtration values when `True`. Default is `True`.
        Returns:
            SimplexTree: A simplex tree created from the Delaunay Triangulation. The vertex `k` corresponds to the k-th
                input point. The vertices may not be numbered contiguously as some points may be discarded in the
                triangulation (duplicate points, weighted hidden point, ...).
        """
        if not filtration in [None, 'alpha', 'cech']:
            raise ValueError(f"\'{filtration}\' is not a valid filtration value. Must be None, \'alpha\' or \'cech\'")

        cdef Delaunay_filtration filt = NONE
        if filtration == 'cech':
            filt = CECH
        elif filtration == 'alpha':
            filt = ALPHA

        stree = SimplexTree()
        cdef intptr_t stree_int_ptr=stree.thisptr

        with nogil:
            self.this_ptr.create_simplex_tree(<Simplex_tree_python_interface*>stree_int_ptr,
                                              max_alpha_square, filt, output_squared_values)
        return stree

    @staticmethod
    def set_float_relative_precision(precision: float):
        """
        Args:
            precision: When constructing :func:`~gudhi.delaunay_cech_complex`, :func:`~gudhi.alpha_complex`, or
                :func:`~gudhi.weighted_alpha_complex` with :code:`precision = 'safe'` (the default), one can
                set the float relative precision of filtration values computed. Default is :code:`1e-5` (cf.
                :func:`~gudhi.DelaunayComplex.get_float_relative_precision`). For more details, please refer to
                `CGAL::Lazy_exact_nt<NT>::set_relative_precision_of_to_double <https://doc.cgal.org/latest/Number_types/classCGAL_1_1Lazy__exact__nt.html>`_

        :raises ValueError: If precision is not in (0, 1).
        """
        if precision <=0. or precision >= 1.:
            raise ValueError("Relative precision value must be strictly greater than 0 and strictly lower than 1")
        Delaunay_complex_interface.set_float_relative_precision(precision)

    @staticmethod
    def get_float_relative_precision() -> float:
        """
        Returns:
            The float relative precision of filtration values computation when constructing with
            :code:`precision = 'safe'` (the default).
        """
        return Delaunay_complex_interface.get_float_relative_precision()


cdef class AlphaComplex(DelaunayComplex):
    """AlphaComplex is a simplicial complex constructed from the finite cells of a Delaunay Triangulation.

    The filtration value of each simplex is computed as the squared radius of the smallest empty sphere passing through
    all of its vertices.

    All simplices that have a filtration value strictly greater than a given alpha squared value are not inserted into
    the complex.

    For more details about the algorithm, please refer to the
    `Alpha complex C++ documentation <https://gudhi.inria.fr/doc/latest/group__alpha__complex.html>`_

    .. note::

        When Alpha Complex is constructed with an infinite value of alpha, the complex is a Delaunay complex.
    """
    def create_simplex_tree(self, double max_alpha_square: float = float('inf'),
                            default_filtration_value: bool = False,
                            bool output_squared_values: bool = True) -> SimplexTree:
        """
        Args:
            max_alpha_square: The maximum alpha square threshold the simplices shall not exceed. Default is set to
                infinity, and there is very little point using anything else since it does not save time.
            default_filtration_value: Default value is `False` (which means compute the filtration values). Set this
                value to `True` if filtration values are not needed to be computed (will be set to `NaN`), but please
                consider constructing a :class:`~gudhi.DelaunayComplex` instead.
            output_squared_values: Square filtration values when `True`. Default is `True` to keep backward
            compatibility.
        Returns:
            SimplexTree: A simplex tree created from the Delaunay Triangulation. The vertex `k` corresponds to the k-th
                input point. The vertices may not be numbered contiguously as some points may be discarded in the
                triangulation (duplicate points, weighted hidden point, ...).
        """
        filtration = 'alpha'
        if default_filtration_value:
            filtration = None
        return super().create_simplex_tree(max_alpha_square, filtration, output_squared_values)

    def get_point(self, vertex: int) -> list[float]:
        """This function returns the point corresponding to a given vertex from the :class:`~gudhi.SimplexTree` (the
        same as the k-th input point, where `k=vertex`)

        Args:
            vertex: The vertex.
        Returns:
            the point.

        :raises IndexError: In case the point has no associated vertex in the diagram (because of weights or because it
            is a duplicate).
        """
        return self.this_ptr.get_point(vertex)


cdef class DelaunayCechComplex(DelaunayComplex):
    """DelaunayCechComplex is a simplicial complex constructed from the finite cells of a Delaunay Triangulation.

    The filtration value of each simplex is equal to the squared radius of its minimal enclosing ball (MEB).

    All simplices that have a filtration value strictly greater than a given alpha squared value are not inserted into
    the complex.

    .. note::

        When DelaunayCechComplex is constructed with an infinite value of alpha, the complex is a Delaunay complex.
    """
    def __init__(self, points: Iterable[Iterable[float]] = [], precision: Literal['fast', 'safe', 'exact'] = 'safe'):
        """
        Args:
            points (Iterable[Iterable[float]]): A list of points in d-Dimension.
            precision (str): Complex precision can be `'fast'`, `'safe'` or `'exact'`. Default is `'safe'`.
        """
        super().__init__(points = points, weights = None, precision = precision)

    def create_simplex_tree(self, double max_alpha_square: float = float('inf'),
                            bool output_squared_values: bool = True) -> SimplexTree:
        """
        Args:
            max_alpha_square: The maximum alpha square threshold the simplices shall not exceed. Default is set to
                infinity, and there is very little point using anything else since it does not save time.
            output_squared_values: Square filtration values when `True`. Default is `True`.
        Returns:
            SimplexTree: A simplex tree created from the Delaunay Čech Triangulation. The vertex `k` corresponds to the
            k-th input point. The vertices may not be numbered contiguously as some points may be discarded in the
            triangulation(duplicate points, ...).
        """
        return super().create_simplex_tree(max_alpha_square, 'cech', output_squared_values)
