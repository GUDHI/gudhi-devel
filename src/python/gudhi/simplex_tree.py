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

from __future__ import annotations


__license__ = "MIT"


import numpy as np
from numpy.typing import ArrayLike
import warnings

from gudhi import _simplex_tree_ext as t


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

    def __init__(self, other=None):
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
                raise TypeError(
                    "`other` argument requires to be of type `SimplexTree`, or `None`."
                )
        else:
            super().__init__()

    def _is_persistence_defined(self) -> bool:
        """Returns `True` if Persistence pointer is not `None`."""
        return self._pers != None

    def copy(self) -> SimplexTree:
        """
        :returns: A simplex tree that is a deep copy of itself.
        :rtype: SimplexTree

        :note: The persistence information is not copied. If you need it in the clone, you have to call
            :func:`compute_persistence` on it even if you had already computed it in the original.
        """
        simplex_tree = SimplexTree(self)
        return simplex_tree

    def __deepcopy__(self) -> SimplexTree:
        return self.copy()

    def initialize_filtration(self):
        """.. deprecated:: 3.2.0
            The initialization is done automatically if necessary; there is no need to call this method since.
        
        This function initializes and sorts the simplicial complex
        filtration vector.
        """
        warnings.warn(
            "Since Gudhi 3.2, calling SimplexTree.initialize_filtration is unnecessary.",
            DeprecationWarning,
        )
        super()._initialize_filtration()

    @staticmethod
    def create_from_array(filtrations, max_filtration: float = float("inf")) -> SimplexTree:
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
        if filtrations.ndim != 2:
            raise ValueError(f"`filtrations` has to be a 2d array. Got {filtrations.ndim=}")
        if filtrations.shape[0] != filtrations.shape[1]:
            raise ValueError(
                f"`filtrations` has to be a square array. Got {filtrations.shape=}"
            )
        ret = SimplexTree()
        ret._insert_matrix(filtrations, max_filtration)
        return ret

    def insert_edges_from_coo_matrix(self, edges) -> SimplexTree:
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
        super()._insert_batch_vertices(
            np.unique(np.stack((edges.row, edges.col))), float("inf")
        )
        # TODO: optimize this?
        for edge in zip(edges.row, edges.col, edges.data):
            super().insert((edge[0], edge[1]), edge[2])
        return self

    def insert_batch(self, vertex_array: ArrayLike, filtrations: ArrayLike) -> SimplexTree:
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
        simplices = np.asarray(vertex_array, dtype=np.intc)
        vertices = np.unique(simplices)
        fil = np.asarray(filtrations, dtype=np.double)
        super()._insert_batch(vertices, simplices, fil)
        return self

    def extended_persistence(
        self, homology_coeff_field=11, min_persistence=0
    ) -> list[list[tuple[int, tuple[float, float]]]]:
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
        self._pers._compute_persistence(homology_coeff_field, -1.0)
        return self._pers._compute_extended_persistence_subdiagrams(min_persistence)

    def persistence(
        self, homology_coeff_field=11, min_persistence=0, persistence_dim_max=False
    ) -> list[tuple[int, tuple[float, float]]]:
        """This function computes and returns the persistence of the simplicial complex.

        :param homology_coeff_field: The homology coefficient field. Must be a
            prime number. Default value is 11. Max is 46337.
        :type homology_coeff_field: int
        :param min_persistence: The minimum persistence value to take into
            account (strictly greater than min_persistence). Default value is
            0.0.
            Set min_persistence to -1.0 to see all values.
        :type min_persistence: float
        :param persistence_dim_max: If `True`, the persistent homology for the
            maximal dimension in the complex is computed. If `False`, it is
            ignored. Default is `False`.
        :type persistence_dim_max: bool
        :returns: The persistence of the simplicial complex.
        :rtype:  list of pairs(dimension, pair(birth, death))
        """
        self.compute_persistence(homology_coeff_field, min_persistence, persistence_dim_max)
        return self._pers._get_persistence()

    def compute_persistence(
        self, homology_coeff_field=11, min_persistence=0, persistence_dim_max=False
    ) -> SimplexTree:
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
        :param persistence_dim_max: If `True`, the persistent homology for the
            maximal dimension in the complex is computed. If `False`, it is
            ignored. Default is `False`.
        :type persistence_dim_max: bool
        :returns: Nothing.
        """
        # bool(persistence_dim_max) because of numpy bool_ is not recognized as the right type
        self._pers = t._Simplex_tree_persistence_interface(self, bool(persistence_dim_max))
        self._pers._compute_persistence(homology_coeff_field, min_persistence)
        return self

    def betti_numbers(self) -> list[int]:
        """This function returns the Betti numbers of the simplicial complex.

        :returns: The Betti numbers ([B0, B1, ..., Bn]).
        :rtype:  list of int

        :note: betti_numbers function requires
            :func:`compute_persistence`
            function to be launched first.
        """
        if self._pers == None:
            raise RuntimeError("compute_persistence() must be called before betti_numbers()")
        return self._pers._betti_numbers()

    def persistent_betti_numbers(self, from_value, to_value) -> list[int]:
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
        if self._pers == None:
            raise RuntimeError(
                "compute_persistence() must be called before persistent_betti_numbers()"
            )
        return self._pers._persistent_betti_numbers(from_value, to_value)

    def persistence_intervals_in_dimension(self, dimension) -> np.ndarray:
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
        if self._pers == None:
            raise RuntimeError(
                "compute_persistence() must be called before persistence_intervals_in_dimension()"
            )
        piid = np.array(self._pers._intervals_in_dimension(dimension))
        # Workaround https://github.com/GUDHI/gudhi-devel/issues/507
        if len(piid) == 0:
            return np.empty(shape=[0, 2])
        return piid

    def persistence_pairs(self) -> list[tuple[list[int], list[int]]]:
        """This function returns a list of persistence birth and death simplices pairs.

        :returns: A list of persistence simplices intervals.
        :rtype:  list of pair of list of int

        :note: persistence_pairs function requires
            :func:`compute_persistence`
            function to be launched first.
        """
        if self._pers == None:
            raise RuntimeError(
                "compute_persistence() must be called before persistence_pairs()"
            )
        return self._pers._persistence_pairs()

    def write_persistence_diagram(self, persistence_file) -> SimplexTree:
        """This function writes the persistence intervals of the simplicial
        complex in a user given file name.

        :param persistence_file: Name of the file.
        :type persistence_file: string

        :note: intervals_in_dim function requires
            :func:`compute_persistence`
            function to be launched first.
        """
        if self._pers == None:
            raise RuntimeError(
                "compute_persistence() must be called before write_persistence_diagram()"
            )
        self._pers._write_output_diagram(persistence_file)
        return self

    def lower_star_persistence_generators(self) -> tuple[list[np.ndarray], list[np.ndarray]]:
        """Assuming this is a lower-star filtration, this function returns the persistence pairs,
        where each simplex is replaced with the vertex that gave it its filtration value.

        :returns: First the regular persistence pairs, grouped by dimension, with one vertex per extremity,
            and second the essential features, grouped by dimension, with one vertex each
        :rtype: Tuple[List[numpy.array[int] of shape (n,2)], List[numpy.array[int] of shape (m,)]]

        :note: lower_star_persistence_generators requires that `persistence()` be called first.
        """
        if self._pers == None:
            raise RuntimeError(
                "lower_star_persistence_generators() requires that persistence() be called first"
            )
        gen = self._pers._lower_star_generators()
        normal = [np.array(d).reshape(-1, 2) for d in gen[0]]
        infinite = [np.array(d) for d in gen[1]]
        return (normal, infinite)

    def flag_persistence_generators(
        self,
    ) -> tuple[np.ndarray, list[np.ndarray], np.ndarray, list[np.ndarray]]:
        """Assuming this is a flag complex, this function returns the persistence pairs,
        where each simplex is replaced with the vertices of the edges that gave it its filtration value.

        :returns: First the regular persistence pairs of dimension 0, with one vertex for birth and two for death;
            then the other regular persistence pairs, grouped by dimension, with 2 vertices per extremity;
            then the connected components, with one vertex each;
            finally the other essential features, grouped by dimension, with 2 vertices for birth.
        :rtype: Tuple[numpy.array[int] of shape (n,3), List[numpy.array[int] of shape (m,4)], numpy.array[int] of shape (l,), List[numpy.array[int] of shape (k,2)]]

        :note: flag_persistence_generators requires that `persistence()` be called first.
        """
        if self._pers == None:
            raise RuntimeError(
                "flag_persistence_generators() requires that persistence() be called first"
            )
        gen = self._pers._flag_generators()
        if len(gen[0]) == 0:
            normal0 = np.empty((0, 3))
            normals = []
        else:
            l = iter(gen[0])
            normal0 = np.array(next(l)).reshape(-1, 3)
            normals = [np.array(d).reshape(-1, 4) for d in l]
        if len(gen[1]) == 0:
            infinite0 = np.empty(0)
            infinites = []
        else:
            l = iter(gen[1])
            infinite0 = np.array(next(l))
            infinites = [np.array(d).reshape(-1, 2) for d in l]
        return (normal0, normals, infinite0, infinites)

    def collapse_edges_of_graph(
        self, nb_iterations: int = 1, inplace: bool = True
    ) -> SimplexTree:
        """Assuming the complex is a graph (simplices of higher dimension are ignored), this method implicitly
        interprets it as the 1-skeleton of a flag complex, and replaces it with another (smaller) graph whose
        expansion has the same persistent homology, using a technique known as edge collapses
        (see :cite:`edgecollapsearxiv`).

        A natural application is to get a simplex tree of dimension 1 from :class:`~gudhi.RipsComplex`,
        then collapse edges, perform :meth:`expansion()` and finally compute persistence
        (cf. :download:`rips_complex_edge_collapse_example.py <../example/rips_complex_edge_collapse_example.py>`).

        :param nb_iterations: The number of edge collapse iterations to perform. If `nb_iterations` is strictly less
            than 1, the method does nothing and returns `self` (independently of the value of `inplace`). Default is 1.
        :type nb_iterations: int
        :param inplace: If `True`, the collapse is done on this simplex tree. Otherwise, the collapse is done on a new
            tree which is then returned. Default is `True`.
        :type inplace: bool
        :returns: `self` (after modifications from collapses) if `inplace` is set to `True` and a new tree resulting
            from the collapses otherwise.
        :rtype: SimplexTree

        .. warning::
            The current simplex tree is assumed to be a graph, that is of maximal dimension 1.
            If it is not the case, all higher dimensional simplices will get lost during the
            reduction process and not be reinserted. In the case of a flag complex, they can be regained by calling
            `expansion(max_dim)` afterwards.
        """
        if nb_iterations < 1:
            return self
        if inplace:
            super()._collapse_edges_inplace(nb_iterations)
            return self
        return super()._collapse_edges(nb_iterations)

    def collapse_edges_of_flag_complex(
        self, nb_iterations: int = 1, max_expansion_dim: int = None, inplace: bool = True
    ) -> SimplexTree:
        """Assuming the complex is a flag complex (e.g. a Rips complex), applies the edge collapses technique
        (see :cite:`edgecollapsearxiv`) to the 1-skeleton of the complex and rebuilds it by extending the skeleton to a
        new flag complex (with the given maximal dimension). The new complex will be smaller but with same persistent
        homology as the original one (at least in the existing dimensions: all cycle classes with higher or equal
        dimension than `max_expansion_dim` will be ignored even if they existed before).

        If `max_expansion_dim` is set to `1` or less, the method is equivalent to :meth:`collapse_edges_of_graph()`.

        :param nb_iterations: The number of edge collapse iterations to perform. If `nb_iterations` is strictly less
            than 1, the method does nothing and returns `self` (independently of the value of `inplace`). Default is 1.
        :type nb_iterations: int
        :param max_expansion_dim: The maximal dimension to which the new complex has to be expended to. If `None`, the
            current dimension is chosen. Note that the final dimension of the new complex can be smaller if no
            higher-dimensional simplex can exist. Default is `None`.
        :type max_expansion_dim: int
        :param inplace: If `True`, the collapse is done on this simplex tree. Otherwise, the collapse is done on a new
            tree which is then returned. Default is `True`.
        :type inplace: bool
        :returns: `self` (after modifications from collapses) if `inplace` is set to `True` and a new tree resulting
            from the collapses otherwise.
        :rtype: SimplexTree

        .. warning::
            The current simplex tree is assumed to be a flag complex. If it is not the case, some information may be
            lost during the re-expansion, i.e., the persistent homology may differ from before.
        """
        if nb_iterations < 1:
            return self
        if max_expansion_dim is None:
            max_expansion_dim = self.dimension()
        if inplace:
            super()._collapse_edges_inplace(nb_iterations)
            if max_expansion_dim > 1:
                super().expansion(max_expansion_dim)
            return self
        collapsed_complex = super()._collapse_edges(nb_iterations)
        if max_expansion_dim > 1:
            collapsed_complex.expansion(max_expansion_dim)
        return collapsed_complex

    def collapse_edges(self, nb_iterations: int = 1) -> SimplexTree:
        """.. deprecated:: 3.12.0 
            The method was renamed :meth:`collapse_edges_of_graph()` to clarify the expected
            pre-conditions. Please use this new name instead.

        :param nb_iterations: The number of edge collapse iterations to perform. If `nb_iterations` is strictly less
            than 1, the method does nothing. Default is 1.
        :type nb_iterations: int
        :returns: `self` (after modifications from the collapses)
        :rtype: SimplexTree
        """
        warnings.warn(
            "Since Gudhi 3.12, `collapse_edges_of_graph(nb_iterations)` should be called instead of"
            + " `collapse_edges(nb_iterations)`. The method was renamed to clarify the pre-condition that the used"
            + " simplex tree should be a graph. Note also the existence of the new method"
            + " `collapse_edges_of_flag_complex`.",
            DeprecationWarning,
        )
        return self.collapse_edges_of_graph(nb_iterations, inplace=True)
