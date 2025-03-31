# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - 2023/02 Vincent Rouvreau: Add serialize/deserialize for pickle feature
#   - 2025/03 Alexis Gobé & Hannah Schreiber: Use nanobind instead of Cython for python bindings.
#   - YYYY/MM Author: Description of the modification

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"


from gudhi import _simplex_tree_ext as t

import numpy as np

# SimplexTree python interface
class SimplexTree(t._Simplex_tree_python_interface):
    """The simplex tree is an efficient and flexible data structure for
    representing general (filtered) simplicial complexes. The data structure
    is described in Jean-Daniel Boissonnat and Clément Maria. The Simplex
    Tree: An Efficient Data Structure for General Simplicial Complexes.
    Algorithmica, pages 1–22, 2014.

    This class is a filtered, with keys, and non contiguous vertices version
    of the simplex tree.
    """

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
        self._pers = None
        if other:
            if isinstance(other, SimplexTree):
                super().__init__(other)
            else:
                raise TypeError("`other` argument requires to be of type `SimplexTree`, or `None`.")
        else:
            super().__init__()

    def _is_persistence_defined(self):
        """Returns true if Persistence pointer is not None.
        """
        return self._pers != None

    def copy(self):
        """ 
        :returns: A simplex tree that is a deep copy of itself.
        :rtype: SimplexTree

        :note: The persistence information is not copied. If you need it in the clone, you have to call
            :func:`compute_persistence` on it even if you had already computed it in the original.
        """
        simplex_tree = SimplexTree(self)
        return simplex_tree

    def __deepcopy__(self):
        return self.copy()

    def initialize_filtration(self):
        """This function initializes and sorts the simplicial complex
        filtration vector.

        .. deprecated:: 3.2.0
        """
        import warnings
        warnings.warn("Since Gudhi 3.2, calling SimplexTree.initialize_filtration is unnecessary.", DeprecationWarning)
        super().initialize_filtration()

    @staticmethod
    def create_from_array(filtrations, max_filtration:float = float('inf')):
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
        assert filtrations.shape[0] == filtrations.shape[1], 'create_from_array() expects a square array'
        ret = SimplexTree()
        ret.insert_matrix(filtrations, max_filtration)
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
        super().insert_batch_vertices(np.unique(np.stack((edges.row, edges.col))), float('inf'))
        # TODO: optimize this?
        for edge in zip(edges.row, edges.col, edges.data):
            super().insert((edge[0], edge[1]), edge[2])

    def insert_batch(self, vertex_array: np.ndarray[np.int32] | np.ndarray[np.int64], filtrations: np.ndarray[np.float32] | np.ndarray[np.float64]):
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
        vertices = np.unique(vertex_array)
        k = vertex_array.shape[0]
        n = vertex_array.shape[1]
        assert filtrations.shape[0] == n, 'inconsistent sizes for vertex_array and filtrations'
        v = []
        # Without this, it could be slow if we end up inserting vertices in a bad order (flat_map).
        # NaN currently does the wrong thing
        super().insert_batch_vertices(vertices, float('inf'))
        for i in range(n):
            for j in range(k):
                v.append(vertex_array[j, i])
            super().insert(v, filtrations[i])
            v.clear()

    def get_simplices(self):
        """This function returns a generator with simplices and their given
        filtration values.

        :returns:  The simplices.
        :rtype:  generator with tuples(simplex, filtration)
        """
        for sh in super().simplex_iter():
            yield super().get_simplex_and_filtration(sh)

    def get_filtration(self):
        """This function returns a generator with simplices and their given
        filtration values sorted by increasing filtration values.

        :returns:  The simplices sorted by increasing filtration values.
        :rtype:  generator with tuples(simplex, filtration)
        """
        for sh in super().filtration_iter():
            yield super().get_simplex_and_filtration(sh)

    def get_skeleton(self, dimension):
        """This function returns a generator with the (simplices of the) skeleton of a maximum given dimension.

        :param dimension: The skeleton dimension value.
        :type dimension: int
        :returns:  The (simplices of the) skeleton of a maximum dimension.
        :rtype:  generator with tuples(simplex, filtration)
        """
        for sh in super().skeleton_iter(dimension):
            yield super().get_simplex_and_filtration(sh)

    def get_boundaries(self, simplex):
        """This function returns a generator with the boundaries of a given N-simplex.
        If you do not need the filtration values, the boundary can also be obtained as
        :code:`itertools.combinations(simplex,len(simplex)-1)`.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int.
        :returns:  The (simplices of the) boundary of a simplex
        :rtype:  generator with tuples(simplex, filtration)
        """
        for sh in super().boundary_iter(simplex):
            yield super().get_simplex_and_filtration(sh)

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
        self._pers = t._Simplex_tree_persistence_interface(self, False)
        self._pers.compute_persistence(homology_coeff_field, -1.)
        return self._pers.compute_extended_persistence_subdiagrams(min_persistence)

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
        return self._pers.get_persistence()

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
        self._pers = t._Simplex_tree_persistence_interface(self, persistence_dim_max)
        self._pers.compute_persistence(homology_coeff_field, min_persistence)

    def betti_numbers(self):
        """This function returns the Betti numbers of the simplicial complex.

        :returns: The Betti numbers ([B0, B1, ..., Bn]).
        :rtype:  list of int

        :note: betti_numbers function requires
            :func:`compute_persistence`
            function to be launched first.
        """
        assert self._pers != None, "compute_persistence() must be called before betti_numbers()"
        return self._pers.betti_numbers()

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
        assert self._pers != None, "compute_persistence() must be called before persistent_betti_numbers()"
        return self._pers.persistent_betti_numbers(from_value, to_value)

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
        assert self._pers != None, "compute_persistence() must be called before persistence_intervals_in_dimension()"
        piid = np.array(self._pers.intervals_in_dimension(dimension))
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
        assert self._pers != None, "compute_persistence() must be called before persistence_pairs()"
        return self._pers.persistence_pairs()

    def write_persistence_diagram(self, persistence_file):
        """This function writes the persistence intervals of the simplicial
        complex in a user given file name.

        :param persistence_file: Name of the file.
        :type persistence_file: string

        :note: intervals_in_dim function requires
            :func:`compute_persistence`
            function to be launched first.
        """
        assert self._pers != None, "compute_persistence() must be called before write_persistence_diagram()"
        self._pers.write_output_diagram(persistence_file)

    def lower_star_persistence_generators(self):
        """Assuming this is a lower-star filtration, this function returns the persistence pairs,
        where each simplex is replaced with the vertex that gave it its filtration value.

        :returns: First the regular persistence pairs, grouped by dimension, with one vertex per extremity,
            and second the essential features, grouped by dimension, with one vertex each
        :rtype: Tuple[List[numpy.array[int] of shape (n,2)], List[numpy.array[int] of shape (m,)]]

        :note: lower_star_persistence_generators requires that `persistence()` be called first.
        """
        assert self._pers != None, "lower_star_persistence_generators() requires that persistence() be called first."
        gen = self._pers.lower_star_generators()
        normal = [np.array(d).reshape(-1,2) for d in gen[0]]
        infinite = [np.array(d) for d in gen[1]]
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
        assert self._pers != None, "flag_persistence_generators() requires that persistence() be called first."
        gen = self._pers.flag_generators()
        if len(gen[0]) == 0:
            normal0 = np.empty((0,3))
            normals = []
        else:
            l = iter(gen[0])
            normal0 = np.array(next(l)).reshape(-1,3)
            normals = [np.array(d).reshape(-1,4) for d in l]
        if len(gen[1]) == 0:
            infinite0 = np.empty(0)
            infinites = []
        else:
            l = iter(gen[1])
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
        if nb_iterations < 1:
            return
        super().collapse_edges(nb_iterations)




