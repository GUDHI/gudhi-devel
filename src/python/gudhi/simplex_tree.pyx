# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from cython.operator import dereference, preincrement
from libc.stdint cimport intptr_t
from numpy import array as np_array
cimport simplex_tree

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

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
    # unfortunately 'cdef public Simplex_tree_interface_full_featured* thisptr' is not possible
    # Use intptr_t instead to cast the pointer
    cdef public intptr_t thisptr

    # Get the pointer casted as it should be
    cdef Simplex_tree_interface_full_featured* get_ptr(self):
        return <Simplex_tree_interface_full_featured*>(self.thisptr)

    cdef Simplex_tree_persistence_interface * pcohptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self):
        """SimplexTree constructor.
        """

    # The real cython constructor
    def __cinit__(self):
        self.thisptr = <intptr_t>(new Simplex_tree_interface_full_featured())

    def __dealloc__(self):
        cdef Simplex_tree_interface_full_featured* ptr = self.get_ptr()
        if ptr != NULL:
            del ptr
        if self.pcohptr != NULL:
            del self.pcohptr

    def __is_defined(self):
        """Returns true if SimplexTree pointer is not NULL.
         """
        return self.get_ptr() != NULL

    def __is_persistence_defined(self):
        """Returns true if Persistence pointer is not NULL.
         """
        return self.pcohptr != NULL

    def filtration(self, simplex):
        """This function returns the filtration value for a given N-simplex in
        this simplicial complex, or +infinity if it is not in the complex.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int.
        :returns:  The simplicial complex filtration value.
        :rtype:  float
        """
        return self.get_ptr().simplex_filtration(simplex)

    def assign_filtration(self, simplex, filtration):
        """This function assigns a new filtration value to a
        given N-simplex.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int.
        :param filtration:  The new filtration value.
        :type filtration:  float

        .. note::
            Beware that after this operation, the structure may not be a valid
            filtration anymore, a simplex could have a lower filtration value
            than one of its faces. Callers are responsible for fixing this
            (with more :meth:`assign_filtration` or
            :meth:`make_filtration_non_decreasing` for instance) before calling
            any function that relies on the filtration property, like
            :meth:`persistence`.
        """
        self.get_ptr().assign_simplex_filtration(simplex, filtration)

    def initialize_filtration(self):
        """This function initializes and sorts the simplicial complex
        filtration vector.

        .. deprecated:: 3.2.0
        """
        self.get_ptr().initialize_filtration()

    def num_vertices(self):
        """This function returns the number of vertices of the simplicial
        complex.

        :returns:  The simplicial complex number of vertices.
        :rtype:  int
        """
        return self.get_ptr().num_vertices()

    def num_simplices(self):
        """This function returns the number of simplices of the simplicial
        complex.

        :returns:  the simplicial complex number of simplices.
        :rtype:  int
        """
        return self.get_ptr().num_simplices()

    def dimension(self):
        """This function returns the dimension of the simplicial complex.

        :returns:  the simplicial complex dimension.
        :rtype:  int

        .. note::

            This function is not constant time because it can recompute
            dimension if required (can be triggered by
            :func:`remove_maximal_simplex`
            or
            :func:`prune_above_filtration`
            methods).
        """
        return self.get_ptr().dimension()

    def upper_bound_dimension(self):
        """This function returns a valid dimension upper bound of the
        simplicial complex.

        :returns:  an upper bound on the dimension of the simplicial complex.
        :rtype:  int
        """
        return self.get_ptr().upper_bound_dimension()

    def set_dimension(self, dimension):
        """This function sets the dimension of the simplicial complex.

        :param dimension: The new dimension value.
        :type dimension: int.

        .. note::

            This function must be used with caution because it disables
            dimension recomputation when required
            (this recomputation can be triggered by
            :func:`remove_maximal_simplex`
            or
            :func:`prune_above_filtration`
            ).
        """
        self.get_ptr().set_dimension(<int>dimension)

    def find(self, simplex):
        """This function returns if the N-simplex was found in the simplicial
        complex or not.

        :param simplex: The N-simplex to find, represented by a list of vertex.
        :type simplex: list of int.
        :returns:  true if the simplex was found, false otherwise.
        :rtype:  bool
        """
        return self.get_ptr().find_simplex(simplex)

    def insert(self, simplex, filtration=0.0):
        """This function inserts the given N-simplex and its subfaces with the
        given filtration value (default value is '0.0'). If some of those
        simplices are already present with a higher filtration value, their
        filtration value is lowered.

        :param simplex: The N-simplex to insert, represented by a list of
            vertex.
        :type simplex: list of int.
        :param filtration: The filtration value of the simplex.
        :type filtration: float.
        :returns:  true if the simplex was not yet in the complex, false
            otherwise (whatever its original filtration value).
        :rtype:  bool
        """
        return self.get_ptr().insert(simplex, <double>filtration)

    def get_simplices(self):
        """This function returns a generator with simplices and their given
        filtration values.

        :returns:  The simplices.
        :rtype:  generator with tuples(simplex, filtration)
        """
        cdef Simplex_tree_simplices_iterator it = self.get_ptr().get_simplices_iterator_begin()
        cdef Simplex_tree_simplices_iterator end = self.get_ptr().get_simplices_iterator_end()
        cdef Simplex_tree_simplex_handle sh = dereference(it)

        while it != end:
            yield self.get_ptr().get_simplex_and_filtration(dereference(it))
            preincrement(it)

    def get_filtration(self):
        """This function returns a generator with simplices and their given
        filtration values sorted by increasing filtration values.

        :returns:  The simplices sorted by increasing filtration values.
        :rtype:  generator with tuples(simplex, filtration)
        """
        cdef vector[Simplex_tree_simplex_handle].const_iterator it = self.get_ptr().get_filtration_iterator_begin()
        cdef vector[Simplex_tree_simplex_handle].const_iterator end = self.get_ptr().get_filtration_iterator_end()

        while it != end:
            yield self.get_ptr().get_simplex_and_filtration(dereference(it))
            preincrement(it)

    def get_skeleton(self, dimension):
        """This function returns the (simplices of the) skeleton of a maximum
        given dimension.

        :param dimension: The skeleton dimension value.
        :type dimension: int.
        :returns:  The (simplices of the) skeleton of a maximum dimension.
        :rtype:  list of tuples(simplex, filtration)
        """
        cdef Simplex_tree_skeleton_iterator it = self.get_ptr().get_skeleton_iterator_begin(dimension)
        cdef Simplex_tree_skeleton_iterator end = self.get_ptr().get_skeleton_iterator_end(dimension)

        while it != end:
            yield self.get_ptr().get_simplex_and_filtration(dereference(it))
            preincrement(it)

    def get_star(self, simplex):
        """This function returns the star of a given N-simplex.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int.
        :returns:  The (simplices of the) star of a simplex.
        :rtype:  list of tuples(simplex, filtration)
        """
        cdef vector[int] csimplex
        for i in simplex:
            csimplex.push_back(i)
        cdef vector[pair[vector[int], double]] star \
            = self.get_ptr().get_star(csimplex)
        ct = []
        for filtered_simplex in star:
            v = []
            for vertex in filtered_simplex.first:
                v.append(vertex)
            ct.append((v, filtered_simplex.second))
        return ct

    def get_cofaces(self, simplex, codimension):
        """This function returns the cofaces of a given N-simplex with a
        given codimension.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int.
        :param codimension: The codimension. If codimension = 0, all cofaces
            are returned (equivalent of get_star function)
        :type codimension: int.
        :returns:  The (simplices of the) cofaces of a simplex
        :rtype:  list of tuples(simplex, filtration)
        """
        cdef vector[int] csimplex
        for i in simplex:
            csimplex.push_back(i)
        cdef vector[pair[vector[int], double]] cofaces \
            = self.get_ptr().get_cofaces(csimplex, <int>codimension)
        ct = []
        for filtered_simplex in cofaces:
            v = []
            for vertex in filtered_simplex.first:
                v.append(vertex)
            ct.append((v, filtered_simplex.second))
        return ct

    def remove_maximal_simplex(self, simplex):
        """This function removes a given maximal N-simplex from the simplicial
        complex.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int.

        .. note::

            The dimension of the simplicial complex may be lower after calling
            remove_maximal_simplex than it was before. However,
            :func:`upper_bound_dimension`
            method will return the old value, which
            remains a valid upper bound. If you care, you can call
            :func:`dimension`
            to recompute the exact dimension.
        """
        self.get_ptr().remove_maximal_simplex(simplex)

    def prune_above_filtration(self, filtration):
        """Prune above filtration value given as parameter.

        :param filtration: Maximum threshold value.
        :type filtration: float.
        :returns: The filtration modification information.
        :rtype: bool


        .. note::

            Note that the dimension of the simplicial complex may be lower
            after calling
            :func:`prune_above_filtration`
            than it was before. However,
            :func:`upper_bound_dimension`
            will return the old value, which remains a
            valid upper bound. If you care, you can call
            :func:`dimension`
            method to recompute the exact dimension.
        """
        return self.get_ptr().prune_above_filtration(filtration)

    def expansion(self, max_dim):
        """Expands the Simplex_tree containing only its one skeleton
        until dimension max_dim.

        The expanded simplicial complex until dimension :math:`d`
        attached to a graph :math:`G` is the maximal simplicial complex of
        dimension at most :math:`d` admitting the graph :math:`G` as
        :math:`1`-skeleton.
        The filtration value assigned to a simplex is the maximal filtration
        value of one of its edges.

        The Simplex_tree must contain no simplex of dimension bigger than
        1 when calling the method.

        :param max_dim: The maximal dimension.
        :type max_dim: int.
        """
        self.get_ptr().expansion(max_dim)

    def make_filtration_non_decreasing(self):
        """This function ensures that each simplex has a higher filtration
        value than its faces by increasing the filtration values.

        :returns: True if any filtration value was modified,
            False if the filtration was already non-decreasing.
        :rtype: bool
        """
        return self.get_ptr().make_filtration_non_decreasing()

    def extend_filtration(self):
        """ Extend filtration for computing extended persistence. This function only uses the 
        filtration values at the 0-dimensional simplices, and computes the extended persistence 
        diagram induced by the lower-star filtration computed with these values. 

        .. note::

            Note that after calling this function, the filtration 
            values are actually modified within the Simplex_tree. 
            The function :func:`extended_persistence`
            retrieves the original values.

        .. note::

            Note that this code creates an extra vertex internally, so you should make sure that
            the Simplex_tree does not contain a vertex with the largest possible value (i.e., 4294967295). 
        """
        self.get_ptr().compute_extended_filtration()

    def extended_persistence(self, homology_coeff_field=11, min_persistence=0):
        """This function retrieves good values for extended persistence, and separate the diagrams 
        into the Ordinary, Relative, Extended+ and Extended- subdiagrams.

        :param homology_coeff_field: The homology coefficient field. Must be a
            prime number. Default value is 11.
        :type homology_coeff_field: int.
        :param min_persistence: The minimum persistence value (i.e., the absolute value of the difference between the persistence diagram point coordinates) to take into
            account (strictly greater than min_persistence). Default value is
            0.0.
            Sets min_persistence to -1.0 to see all values.
        :type min_persistence: float.
        :returns: A list of four persistence diagrams in the format described in :func:`persistence`. The first one is Ordinary, the second one is Relative, the third one is Extended+ and the fourth one is Extended-. See https://link.springer.com/article/10.1007/s10208-008-9027-z and/or section 2.2 in https://link.springer.com/article/10.1007/s10208-017-9370-z for a description of these subtypes.

        .. note::

            This function should be called only if :func:`extend_filtration` has been called first!

        .. note::

            The coordinates of the persistence diagram points might be a little different than the
            original filtration values due to the internal transformation (scaling to [-2,-1]) that is 
            performed on these values during the computation of extended persistence.
        """
        cdef vector[pair[int, pair[double, double]]] persistence_result
        if self.pcohptr != NULL:
            del self.pcohptr
        self.pcohptr = new Simplex_tree_persistence_interface(self.get_ptr(), False)
        self.pcohptr.compute_persistence(homology_coeff_field, -1.)
        persistence_result = self.pcohptr.get_persistence()
        return self.get_ptr().compute_extended_persistence_subdiagrams(persistence_result, min_persistence)


    def persistence(self, homology_coeff_field=11, min_persistence=0, persistence_dim_max = False):
        """This function computes and returns the persistence of the simplicial complex.

        :param homology_coeff_field: The homology coefficient field. Must be a
            prime number. Default value is 11.
        :type homology_coeff_field: int.
        :param min_persistence: The minimum persistence value to take into
            account (strictly greater than min_persistence). Default value is
            0.0.
            Sets min_persistence to -1.0 to see all values.
        :type min_persistence: float.
        :param persistence_dim_max: If true, the persistent homology for the
            maximal dimension in the complex is computed. If false, it is
            ignored. Default is false.
        :type persistence_dim_max: bool
        :returns: The persistence of the simplicial complex.
        :rtype:  list of pairs(dimension, pair(birth, death))
        """
        self.compute_persistence(homology_coeff_field, min_persistence, persistence_dim_max)
        return self.pcohptr.get_persistence()

    def compute_persistence(self, homology_coeff_field=11, min_persistence=0, persistence_dim_max = False):
        """This function computes the persistence of the simplicial complex, so it can be accessed through
        :func:`persistent_betti_numbers`, :func:`persistence_pairs`, etc. This function is equivalent to :func:`persistence`
        when you do not want the list :func:`persistence` returns.

        :param homology_coeff_field: The homology coefficient field. Must be a
            prime number. Default value is 11.
        :type homology_coeff_field: int.
        :param min_persistence: The minimum persistence value to take into
            account (strictly greater than min_persistence). Default value is
            0.0.
            Sets min_persistence to -1.0 to see all values.
        :type min_persistence: float.
        :param persistence_dim_max: If true, the persistent homology for the
            maximal dimension in the complex is computed. If false, it is
            ignored. Default is false.
        :type persistence_dim_max: bool
        :returns: Nothing.
        """
        if self.pcohptr != NULL:
            del self.pcohptr
        self.pcohptr = new Simplex_tree_persistence_interface(self.get_ptr(), persistence_dim_max)
        self.pcohptr.compute_persistence(homology_coeff_field, min_persistence)

    def betti_numbers(self):
        """This function returns the Betti numbers of the simplicial complex.

        :returns: The Betti numbers ([B0, B1, ..., Bn]).
        :rtype:  list of int

        :note: betti_numbers function requires
            :func:`compute_persistence`
            function to be launched first.
        """
        assert self.pcohptr != NULL, "compute_persistence() must be called before betti_numbers()"
        return self.pcohptr.betti_numbers()

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

        :note: persistent_betti_numbers function requires
            :func:`compute_persistence`
            function to be launched first.
        """
        assert self.pcohptr != NULL, "compute_persistence() must be called before persistent_betti_numbers()"
        return self.pcohptr.persistent_betti_numbers(<double>from_value, <double>to_value)

    def persistence_intervals_in_dimension(self, dimension):
        """This function returns the persistence intervals of the simplicial
        complex in a specific dimension.

        :param dimension: The specific dimension.
        :type dimension: int.
        :returns: The persistence intervals.
        :rtype:  numpy array of dimension 2

        :note: intervals_in_dim function requires
            :func:`compute_persistence`
            function to be launched first.
        """
        assert self.pcohptr != NULL, "compute_persistence() must be called before persistence_intervals_in_dimension()"
        return np_array(self.pcohptr.intervals_in_dimension(dimension))

    def persistence_pairs(self):
        """This function returns a list of persistence birth and death simplices pairs.

        :returns: A list of persistence simplices intervals.
        :rtype:  list of pair of list of int

        :note: persistence_pairs function requires
            :func:`compute_persistence`
            function to be launched first.
        """
        assert self.pcohptr != NULL, "compute_persistence() must be called before persistence_pairs()"
        return self.pcohptr.persistence_pairs()

    def write_persistence_diagram(self, persistence_file):
        """This function writes the persistence intervals of the simplicial
        complex in a user given file name.

        :param persistence_file: Name of the file.
        :type persistence_file: string.

        :note: intervals_in_dim function requires
            :func:`compute_persistence`
            function to be launched first.
        """
        assert self.pcohptr != NULL, "compute_persistence() must be called before write_persistence_diagram()"
        self.pcohptr.write_output_diagram(persistence_file.encode('utf-8'))
