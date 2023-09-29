# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - 2023/02 Vincent Rouvreau: Add serialize/deserialize for pickle feature
#   - YYYY/MM Author: Description of the modification

from cython.operator import dereference, preincrement
from libc.stdint cimport intptr_t, int32_t, int64_t
import numpy as np
cimport gudhi.simplex_tree
cimport cython
from numpy.math cimport INFINITY

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

ctypedef fused some_int:
    int32_t
    int64_t

ctypedef fused some_float:
    float
    double

cdef bool callback(vector[int] simplex, void *blocker_func):
    return (<object>blocker_func)(simplex)

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
    # unfortunately 'cdef public Simplex_tree_python_interface* thisptr' is not possible
    # Use intptr_t instead to cast the pointer
    cdef public intptr_t thisptr

    # Get the pointer casted as it should be
    cdef Simplex_tree_python_interface* get_ptr(self) nogil:
        return <Simplex_tree_python_interface*>(self.thisptr)

    cdef Simplex_tree_persistence_interface * pcohptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, other = None):
        """SimplexTree constructor.

        :param other: If `other` is `None` (default value), an empty `SimplexTree` is created.
            If `other` is a `SimplexTree`, the `SimplexTree` is constructed from a deep copy of `other`.
        :type other: SimplexTree (Optional)
        :returns: An empty or a copy simplex tree.
        :rtype: SimplexTree

        :raises TypeError: In case `other` is neither `None`, nor a `SimplexTree`.
        :note: If the `SimplexTree` is a copy, the persistence information is not copied. If you need it in the clone,
            you have to call :func:`compute_persistence` on it even if you had already computed it in the original.
        """

    # The real cython constructor
    def __cinit__(self, other = None):
        if other:
            if isinstance(other, SimplexTree):
                self.thisptr = _get_copy_intptr(other)
            else:
                raise TypeError("`other` argument requires to be of type `SimplexTree`, or `None`.")
        else:
            self.thisptr = <intptr_t>(new Simplex_tree_python_interface())

    def __dealloc__(self):
        cdef Simplex_tree_python_interface* ptr = self.get_ptr()
        if ptr != NULL:
            del ptr
        if self.pcohptr != NULL:
            del self.pcohptr

    def _is_defined(self):
        """Returns true if SimplexTree pointer is not NULL.
         """
        return self.get_ptr() != NULL

    def _is_persistence_defined(self):
        """Returns true if Persistence pointer is not NULL.
         """
        return self.pcohptr != NULL

    def copy(self):
        """ 
        :returns: A simplex tree that is a deep copy of itself.
        :rtype: SimplexTree

        :note: The persistence information is not copied. If you need it in the clone, you have to call
            :func:`compute_persistence` on it even if you had already computed it in the original.
        """
        stree = SimplexTree()
        stree.thisptr = _get_copy_intptr(self)
        return stree

    def __deepcopy__(self):
        return self.copy()

    def filtration(self, simplex):
        """This function returns the filtration value for a given N-simplex in
        this simplicial complex, or +infinity if it is not in the complex.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int
        :returns:  The simplicial complex filtration value.
        :rtype:  float
        """
        return self.get_ptr().simplex_filtration(simplex)

    def assign_filtration(self, simplex, filtration):
        """This function assigns a new filtration value to a
        given N-simplex.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int
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
        import warnings
        warnings.warn("Since Gudhi 3.2, calling SimplexTree.initialize_filtration is unnecessary.", DeprecationWarning)
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

    def is_empty(self):
        """This function returns whether the simplicial complex is empty.

        :returns:  True if the simplicial complex is empty.
        :rtype:  bool
        """
        return self.get_ptr().is_empty()

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
        :type dimension: int

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
        :type simplex: list of int
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
        :type simplex: list of int
        :param filtration: The filtration value of the simplex.
        :type filtration: float
        :returns:  true if the simplex was not yet in the complex, false
            otherwise (whatever its original filtration value).
        :rtype:  bool
        """
        return self.get_ptr().insert(simplex, <double>filtration)

    @staticmethod
    @cython.boundscheck(False)
    def create_from_array(filtrations, double max_filtration=INFINITY):
        """Creates a new, empty complex and inserts vertices and edges. The vertices are numbered from 0 to n-1, and
        the filtration values are encoded in the array, with the diagonal representing the vertices. It is the
        caller's responsibility to ensure that this defines a filtration, which can be achieved with either::

            filtrations[np.diag_indices_from(filtrations)] = filtrations.min(axis=1)

        or::

            diag = filtrations.diagonal()
            filtrations = np.fmax(np.fmax(filtrations, diag[:, None]), diag[None, :])

        :param filtrations: the filtration values of the vertices and edges to insert. The matrix is assumed to be symmetric.
        :type filtrations: numpy.ndarray of shape (n,n)
        :param max_filtration: only insert vertices and edges with filtration values no larger than max_filtration
        :type max_filtration: float
        :returns: the new complex
        :rtype: SimplexTree
        """
        # TODO: document which half of the matrix is actually read?
        filtrations = np.asanyarray(filtrations, dtype=float)
        cdef double[:,:] F = filtrations
        ret = SimplexTree()
        cdef int n = F.shape[0]
        assert n == F.shape[1], 'create_from_array() expects a square array'
        with nogil:
            ret.get_ptr().insert_matrix(&F[0,0], n, F.strides[0], F.strides[1], max_filtration)
        return ret

    def insert_edges_from_coo_matrix(self, edges):
        """Inserts edges given by a sparse matrix in `COOrdinate format
        <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html>`_.
        If an edge is repeated, the smallest filtration value is used. Missing entries are not inserted.
        Diagonal entries are currently interpreted as vertices, although we do not guarantee this behavior
        in the future, and this is only useful if you want to insert vertices with a smaller filtration value
        than the smallest edge containing it, since vertices are implicitly inserted together with the edges.

        :param edges: the edges to insert and their filtration values.
        :type edges: scipy.sparse.coo_matrix of shape (n,n)

        .. seealso:: :func:`insert_batch`
        """
        # Without this, it could be slow if we end up inserting vertices in a bad order (flat_map).
        self.get_ptr().insert_batch_vertices(np.unique(np.stack((edges.row, edges.col))), INFINITY)
        # TODO: optimize this?
        for edge in zip(edges.row, edges.col, edges.data):
            self.get_ptr().insert((edge[0], edge[1]), edge[2])

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def insert_batch(self, some_int[:,:] vertex_array, some_float[:] filtrations):
        """Inserts k-simplices given by a sparse array in a format similar
        to `torch.sparse <https://pytorch.org/docs/stable/sparse.html>`_.
        The n-th simplex has vertices `vertex_array[0,n]`, ...,
        `vertex_array[k,n]` and filtration value `filtrations[n]`.
        If a simplex is repeated, the smallest filtration value is used.
        Simplices with a repeated vertex are currently interpreted as lower
        dimensional simplices, but we do not guarantee this behavior in the
        future. Any time a simplex is inserted, its faces are inserted as well
        if needed to preserve a simplicial complex.

        :param vertex_array: the k-simplices to insert.
        :type vertex_array: numpy.array of shape (k+1,n)
        :param filtrations: the filtration values.
        :type filtrations: numpy.array of shape (n,)
        """
        cdef vector[int] vertices = np.unique(vertex_array)
        cdef Py_ssize_t k = vertex_array.shape[0]
        cdef Py_ssize_t n = vertex_array.shape[1]
        assert filtrations.shape[0] == n, 'inconsistent sizes for vertex_array and filtrations'
        cdef Py_ssize_t i
        cdef Py_ssize_t j
        cdef vector[int] v
        with nogil:
            # Without this, it could be slow if we end up inserting vertices in a bad order (flat_map).
            # NaN currently does the wrong thing
            self.get_ptr().insert_batch_vertices(vertices, INFINITY)
            for i in range(n):
                for j in range(k):
                    v.push_back(vertex_array[j, i])
                self.get_ptr().insert(v, filtrations[i])
                v.clear()

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
        """This function returns a generator with the (simplices of the) skeleton of a maximum given dimension.

        :param dimension: The skeleton dimension value.
        :type dimension: int
        :returns:  The (simplices of the) skeleton of a maximum dimension.
        :rtype:  generator with tuples(simplex, filtration)
        """
        cdef Simplex_tree_skeleton_iterator it = self.get_ptr().get_skeleton_iterator_begin(dimension)
        cdef Simplex_tree_skeleton_iterator end = self.get_ptr().get_skeleton_iterator_end(dimension)

        while it != end:
            yield self.get_ptr().get_simplex_and_filtration(dereference(it))
            preincrement(it)

    def get_star(self, simplex):
        """This function returns the star of a given N-simplex.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int
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
        :type simplex: list of int
        :param codimension: The codimension. If codimension = 0, all cofaces
            are returned (equivalent of get_star function)
        :type codimension: int
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

    def get_boundaries(self, simplex):
        """This function returns a generator with the boundaries of a given N-simplex.
        If you do not need the filtration values, the boundary can also be obtained as
        :code:`itertools.combinations(simplex,len(simplex)-1)`.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int.
        :returns:  The (simplices of the) boundary of a simplex
        :rtype:  generator with tuples(simplex, filtration)
        """
        cdef pair[Simplex_tree_boundary_iterator, Simplex_tree_boundary_iterator] it =  self.get_ptr().get_boundary_iterators(simplex)

        while it.first != it.second:
            yield self.get_ptr().get_simplex_and_filtration(dereference(it.first))
            preincrement(it.first)

    def remove_maximal_simplex(self, simplex):
        """This function removes a given maximal N-simplex from the simplicial
        complex.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int

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
        :type filtration: float
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

    def prune_above_dimension(self, dimension):
        """Remove all simplices of dimension greater than a given value.

        :param dimension: Maximum dimension value.
        :type dimension: int
        :returns: The modification information.
        :rtype: bool
        """
        return self.get_ptr().prune_above_dimension(dimension)

    def expansion(self, max_dimension):
        """Expands the simplex tree containing only its one skeleton
        until dimension max_dim.

        The expanded simplicial complex until dimension :math:`d`
        attached to a graph :math:`G` is the maximal simplicial complex of
        dimension at most :math:`d` admitting the graph :math:`G` as
        :math:`1`-skeleton.
        The filtration value assigned to a simplex is the maximal filtration
        value of one of its edges.

        The simplex tree must contain no simplex of dimension bigger than
        1 when calling the method.

        :param max_dimension: The maximal dimension.
        :type max_dimension: int
        """
        cdef int maxdim = max_dimension
        with nogil:
            self.get_ptr().expansion(maxdim)

    def make_filtration_non_decreasing(self):
        """This function ensures that each simplex has a higher filtration
        value than its faces by increasing the filtration values.

        :returns: True if any filtration value was modified,
            False if the filtration was already non-decreasing.
        :rtype: bool
        """
        return self.get_ptr().make_filtration_non_decreasing()

    def reset_filtration(self, filtration, min_dim = 0):
        """This function resets the filtration value of all the simplices of dimension at least min_dim. Resets all the
        simplex tree when `min_dim = 0`.
        `reset_filtration` may break the filtration property with `min_dim > 0`, and it is the user's responsibility to
        make it a valid filtration (using a large enough `filt_value`, or calling `make_filtration_non_decreasing`
        afterwards for instance).

        :param filtration: New threshold value.
        :type filtration: float.
        :param min_dim: The minimal dimension. Default value is 0.
        :type min_dim: int.
        """
        self.get_ptr().reset_filtration(filtration, min_dim)

    def extend_filtration(self):
        """ Extend filtration for computing extended persistence. This function only uses the filtration values at the
        0-dimensional simplices, and computes the extended persistence diagram induced by the lower-star filtration
        computed with these values.

        .. note::

            Note that after calling this function, the filtration values are actually modified within the simplex tree.
            The function :func:`extended_persistence` retrieves the original values.

        .. note::

            Note that this code creates an extra vertex internally, so you should make sure that the simplex tree does
            not contain a vertex with the largest possible value (i.e., 4294967295).

        This `notebook <https://github.com/GUDHI/TDA-tutorial/blob/master/Tuto-GUDHI-extended-persistence.ipynb>`_
        explains how to compute an extension of persistence called extended persistence.
        """
        self.get_ptr().compute_extended_filtration()

    def extended_persistence(self, homology_coeff_field=11, min_persistence=0):
        """This function retrieves good values for extended persistence, and separate the diagrams into the Ordinary,
        Relative, Extended+ and Extended- subdiagrams.

        :param homology_coeff_field: The homology coefficient field. Must be a prime number. Default value is 11. Max is 46337.
        :type homology_coeff_field: int
        :param min_persistence: The minimum persistence value (i.e., the absolute value of the difference between the
            persistence diagram point coordinates) to take into account (strictly greater than min_persistence).
            Default value is 0.0. Sets min_persistence to -1.0 to see all values.
        :type min_persistence: float
        :returns: A list of four persistence diagrams in the format described in :func:`persistence`. The first one is
            Ordinary, the second one is Relative, the third one is Extended+ and the fourth one is Extended-.
            See https://link.springer.com/article/10.1007/s10208-008-9027-z and/or section 2.2 in
            https://link.springer.com/article/10.1007/s10208-017-9370-z for a description of these subtypes.

        .. note::

            This function should be called only if :func:`extend_filtration` has been called first!

        .. note::

            The coordinates of the persistence diagram points might be a little different than the
            original filtration values due to the internal transformation (scaling to [-2,-1]) that is 
            performed on these values during the computation of extended persistence.

        This `notebook <https://github.com/GUDHI/TDA-tutorial/blob/master/Tuto-GUDHI-extended-persistence.ipynb>`_
        explains how to compute an extension of persistence called extended persistence.
        """
        cdef vector[pair[int, pair[double, double]]] persistence_result
        if self.pcohptr != NULL:
            del self.pcohptr
        self.pcohptr = new Simplex_tree_persistence_interface(self.get_ptr(), False)
        self.pcohptr.compute_persistence(homology_coeff_field, -1.)
        return self.pcohptr.compute_extended_persistence_subdiagrams(min_persistence)

    def expansion_with_blocker(self, max_dim, blocker_func):
        """Expands the Simplex_tree containing only a graph. Simplices corresponding to cliques in the graph are added
        incrementally, faces before cofaces, unless the simplex has dimension larger than `max_dim` or `blocker_func`
        returns `True` for this simplex.

        The function identifies a candidate simplex whose faces are all already in the complex, inserts it with a
        filtration value corresponding to the maximum of the filtration values of the faces, then calls `blocker_func`
        with this new simplex (represented as a list of int). If `blocker_func` returns `True`, the simplex is removed,
        otherwise it is kept. The algorithm then proceeds with the next candidate.

        .. warning::
            Several candidates of the same dimension may be inserted simultaneously before calling `blocker_func`, so
            if you examine the complex in `blocker_func`, you may hit a few simplices of the same dimension that have
            not been vetted by `blocker_func` yet, or have already been rejected but not yet removed.

        :param max_dim: Expansion maximal dimension value.
        :type max_dim: int
        :param blocker_func: Blocker oracle.
        :type blocker_func: Callable[[List[int]], bool]
        """
        self.get_ptr().expansion_with_blockers_callback(max_dim, callback, <void*>blocker_func)

    def persistence(self, homology_coeff_field=11, min_persistence=0, persistence_dim_max = False):
        """This function computes and returns the persistence of the simplicial complex.

        :param homology_coeff_field: The homology coefficient field. Must be a
            prime number. Default value is 11. Max is 46337.
        :type homology_coeff_field: int
        :param min_persistence: The minimum persistence value to take into
            account (strictly greater than min_persistence). Default value is
            0.0.
            Set min_persistence to -1.0 to see all values.
        :type min_persistence: float
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
            prime number. Default value is 11. Max is 46337.
        :type homology_coeff_field: int
        :param min_persistence: The minimum persistence value to take into
            account (strictly greater than min_persistence). Default value is
            0.0.
            Sets min_persistence to -1.0 to see all values.
        :type min_persistence: float
        :param persistence_dim_max: If true, the persistent homology for the
            maximal dimension in the complex is computed. If false, it is
            ignored. Default is false.
        :type persistence_dim_max: bool
        :returns: Nothing.
        """
        if self.pcohptr != NULL:
            del self.pcohptr
        cdef bool pdm = persistence_dim_max
        cdef int coef = homology_coeff_field
        cdef double minp = min_persistence
        with nogil:
            self.pcohptr = new Simplex_tree_persistence_interface(self.get_ptr(), pdm)
            self.pcohptr.compute_persistence(coef, minp)

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
        :type from_value: float
        :param to_value: The persistence death limit to be added in the
            numbers (persistent death > to_value).
        :type to_value: float

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
        :type dimension: int
        :returns: The persistence intervals.
        :rtype:  numpy array of dimension 2

        :note: intervals_in_dim function requires
            :func:`compute_persistence`
            function to be launched first.
        """
        assert self.pcohptr != NULL, "compute_persistence() must be called before persistence_intervals_in_dimension()"
        piid = np.array(self.pcohptr.intervals_in_dimension(dimension))
        # Workaround https://github.com/GUDHI/gudhi-devel/issues/507
        if len(piid) == 0:
            return np.empty(shape = [0, 2])
        return piid

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
        :type persistence_file: string

        :note: intervals_in_dim function requires
            :func:`compute_persistence`
            function to be launched first.
        """
        assert self.pcohptr != NULL, "compute_persistence() must be called before write_persistence_diagram()"
        self.pcohptr.write_output_diagram(persistence_file.encode('utf-8'))

    def lower_star_persistence_generators(self):
        """Assuming this is a lower-star filtration, this function returns the persistence pairs,
        where each simplex is replaced with the vertex that gave it its filtration value.

        :returns: First the regular persistence pairs, grouped by dimension, with one vertex per extremity,
            and second the essential features, grouped by dimension, with one vertex each
        :rtype: Tuple[List[numpy.array[int] of shape (n,2)], List[numpy.array[int] of shape (m,)]]

        :note: lower_star_persistence_generators requires that `persistence()` be called first.
        """
        assert self.pcohptr != NULL, "lower_star_persistence_generators() requires that persistence() be called first."
        gen = self.pcohptr.lower_star_generators()
        normal = [np.array(d).reshape(-1,2) for d in gen.first]
        infinite = [np.array(d) for d in gen.second]
        return (normal, infinite)

    def flag_persistence_generators(self):
        """Assuming this is a flag complex, this function returns the persistence pairs,
        where each simplex is replaced with the vertices of the edges that gave it its filtration value.

        :returns: First the regular persistence pairs of dimension 0, with one vertex for birth and two for death;
            then the other regular persistence pairs, grouped by dimension, with 2 vertices per extremity;
            then the connected components, with one vertex each;
            finally the other essential features, grouped by dimension, with 2 vertices for birth.
        :rtype: Tuple[numpy.array[int] of shape (n,3), List[numpy.array[int] of shape (m,4)], numpy.array[int] of shape (l,), List[numpy.array[int] of shape (k,2)]]

        :note: flag_persistence_generators requires that `persistence()` be called first.
        """
        assert self.pcohptr != NULL, "flag_persistence_generators() requires that persistence() be called first."
        gen = self.pcohptr.flag_generators()
        if len(gen.first) == 0:
            normal0 = np.empty((0,3))
            normals = []
        else:
            l = iter(gen.first)
            normal0 = np.array(next(l)).reshape(-1,3)
            normals = [np.array(d).reshape(-1,4) for d in l]
        if len(gen.second) == 0:
            infinite0 = np.empty(0)
            infinites = []
        else:
            l = iter(gen.second)
            infinite0 = np.array(next(l))
            infinites = [np.array(d).reshape(-1,2) for d in l]
        return (normal0, normals, infinite0, infinites)

    def collapse_edges(self, nb_iterations = 1):
        """Assuming the complex is a graph (simplices of higher dimension are ignored), this method implicitly
        interprets it as the 1-skeleton of a flag complex, and replaces it with another (smaller) graph whose
        expansion has the same persistent homology, using a technique known as edge collapses
        (see :cite:`edgecollapsearxiv`).

        A natural application is to get a simplex tree of dimension 1 from :class:`~gudhi.RipsComplex`,
        then collapse edges, perform :meth:`expansion()` and finally compute persistence
        (cf. :download:`rips_complex_edge_collapse_example.py <../example/rips_complex_edge_collapse_example.py>`).

        :param nb_iterations: The number of edge collapse iterations to perform. Default is 1.
        :type nb_iterations: int
        """
        # Backup old pointer
        cdef Simplex_tree_python_interface* ptr = self.get_ptr()
        cdef int nb_iter = nb_iterations
        with nogil:
            # New pointer is a new collapsed simplex tree
            self.thisptr = <intptr_t>(ptr.collapse_edges(nb_iter))
            # Delete old pointer
            del ptr

    def __eq__(self, other:SimplexTree):
        """:returns: True if the 2 complexes have the same simplices with the same filtration values, False otherwise.
        :rtype: bool
        """
        return dereference(self.get_ptr()) == dereference(other.get_ptr())
    
    def __getstate__(self):
        """:returns: Serialized (or flattened) SimplexTree data structure in order to pickle SimplexTree.
        :rtype: numpy.array of shape (n,)
        """
        cdef size_t buffer_size = self.get_ptr().get_serialization_size()
        # Let's use numpy to allocate a buffer. Will be deleted automatically
        np_buffer = np.empty(buffer_size, dtype='B')
        cdef char[:] buffer = np_buffer
        cdef char* buffer_start = &buffer[0]
        with nogil:
            self.get_ptr().serialize(buffer_start, buffer_size)
        
        return np_buffer

    def __setstate__(self, state):
        """Construct the SimplexTree data structure from a Numpy Array (cf. :func:`~gudhi.SimplexTree.__getstate__`)
        in order to unpickle a SimplexTree.
        
        :param state: Serialized SimplexTree data structure
        :type state: numpy.array of shape (n,)
        """
        cdef char[:] buffer = state
        cdef size_t buffer_size = state.shape[0]
        cdef char* buffer_start = &buffer[0]
        # Delete pointer, just in case, as deserialization requires an empty SimplexTree
        cdef Simplex_tree_python_interface* ptr = self.get_ptr()
        del ptr
        self.thisptr = <intptr_t>(new Simplex_tree_python_interface())
        with nogil:
            # New pointer is a deserialized simplex tree
            self.get_ptr().deserialize(buffer_start, buffer_size)


cdef intptr_t _get_copy_intptr(SimplexTree stree) nogil:
    return <intptr_t>(new Simplex_tree_python_interface(dereference(stree.get_ptr())))
