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


def _truncate_per_dim(finites: list, infinites: list, max_dim: int) -> tuple:
    """Slice per-dim ``(finites, infinites)`` lists down to ``0..max_dim``."""
    cutoff = max_dim + 1
    return finites[:cutoff], infinites[:cutoff]


def _copy_per_dim(finites: list, infinites: list) -> tuple:
    """Defensive deep-copy used by the Sq^0 = identity short-circuit.

    The numpy arrays are duplicated so that the ordinary and Steenrod
    halves of the returned pair don't alias.
    """
    return ([f.copy() for f in finites], [i.copy() for i in infinites])


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
        """Returns true if Persistence pointer is not None."""
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
        """This function initializes and sorts the simplicial complex
        filtration vector.

        .. deprecated:: 3.2.0
        """
        import warnings

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
        :param persistence_dim_max: If true, the persistent homology for the
            maximal dimension in the complex is computed. If false, it is
            ignored. Default is false.
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
        :param persistence_dim_max: If true, the persistent homology for the
            maximal dimension in the complex is computed. If false, it is
            ignored. Default is false.
        :type persistence_dim_max: bool
        :returns: Nothing.
        """
        # bool(persistence_dim_max) because of numpy bool_ is not recognized as the right type
        self._pers = t._Simplex_tree_persistence_interface(self, bool(persistence_dim_max))
        self._pers._compute_persistence(homology_coeff_field, min_persistence)
        return self

    def compute_steenrod_barcodes(self, k: int = 1, absolute: bool = False, max_dim: int | None = None, n_jobs: int = -1):
        """Compute the ordinary persistence barcode together with the
        Sq\\ :sup:`k` Steenrod barcode of the filtered simplicial complex.

        Coefficients are in :math:`\\mathbb{F}_2` — Steenrod squares are only
        defined there, so no field parameter is exposed. This method does
        **not** require a previous call to :func:`compute_persistence`; the
        reduction needed to obtain cocycle representatives is run internally.

        :param k: the Steenrod square exponent (non-negative integer).
            Default ``1``.  ``k = 0`` is supported and returns the ordinary
            barcode in both outputs (Sq⁰ is the identity); negative values
            raise :class:`ValueError`.
        :type k: int
        :param absolute: if ``False`` (default), return bars in the
            **relative cohomology convention** of Lupo, Medina-Mardones,
            Tauzin (2022) §2.4 — see ``:returns:`` for the on-the-wire
            format.  If ``True``, return bars in the **absolute cohomology
            convention**: finite bars at relative dimension ``d`` shift
            down to absolute dimension ``d - 1``; essentials stay where
            they are; the numerical values do not change.  The absolute
            interpretation rests on a duality bijection that requires the
            relative ordinary barcode to have no essential bars at degrees
            ``[1, max_dim]``; a :class:`UserWarning` is emitted when this
            condition is not satisfied.
        :type absolute: bool
        :param max_dim: if not ``None``, truncate the returned barcodes to
            dimensions ``0..max_dim``.  Higher-dimensional simplices are
            still used internally for the reduction (so that ``H^max_dim``
            classes can be killed by simplices in dimension ``max_dim + 1``),
            but no bars in those higher dimensions are reported.  Default
            ``None`` returns all dimensions of the complex.
        :type max_dim: int or None
        :param n_jobs: number of OpenMP threads to use for the parallelised
            stages (``compute_steenrod_matrix`` and
            ``compute_steenrod_barcode``). ``-1`` (default) uses all available
            cores. Has no effect if the library was built without OpenMP.
        :type n_jobs: int
        :returns: A pair ``(ordinary, steenrod)``.  Each element is itself a
            pair ``(finites, infinites)`` of per-dimension lists of numpy
            arrays:

            * ``finites[d]`` has shape ``(n_bars, 2)``; rows are
              ``(death_value, birth_value)`` with ``death < birth``,
              encoding the relative-cohomology bar ``[a_p, a_{q+1})`` of
              Lupo, Medina-Mardones, Tauzin (2022) §2.4.
            * ``infinites[d]`` has shape ``(n_bars,)``; entries are the
              birth values of essential bars (the implicit lower endpoint
              of the relative-cohomology interval is ``-inf``).

            For Steenrod, ``finites[d]`` and ``infinites[d]`` are empty for
            ``d < k + 1`` — the first ``k`` relative dimensions are
            produced empty by the algorithm.
        :rtype: tuple(tuple(list, list), tuple(list, list))

        Example (:math:`\\mathbb{R}P^2`, Sq¹)::

            from gudhi import SimplexTree

            # Minimal triangulation of RP^2: 6 vertices, 15 edges, 10 triangles.
            rp2_top = [
                [1, 2, 4], [2, 3, 4], [1, 3, 5], [2, 3, 5], [1, 4, 5],
                [1, 2, 6], [1, 3, 6], [3, 4, 6], [2, 5, 6], [4, 5, 6],
            ]

            st = SimplexTree()
            for tri in rp2_top:
                st.insert(tri, filtration=0.0)

            ordinary, steenrod = st.compute_steenrod_barcodes(k=1)
            (st_finites, st_infinites) = steenrod
            # Sq^1: H^1(RP^2) -> H^2(RP^2) is an isomorphism over F_2, so
            # st_infinites[2] holds exactly one essential bar (the upper
            # endpoint of the relative interval; the implicit lower
            # endpoint is -inf).
            assert st_infinites[2].shape == (1,)
        """
        iface = t._Steenrod_barcode_interface(self, int(k))

        if absolute:
            ordinary, steenrod, problematic_dims = iface._compute_absolute(
                int(n_jobs), -1 if max_dim is None else int(max_dim))
            if problematic_dims:
                dims_str = ", ".join(f"H^{d}" for d in problematic_dims)
                warnings.warn(
                    f"absolute=True: the duality condition is not satisfied "
                    f"(ordinary {dims_str} has essential bars in relative "
                    f"convention).  The absolute interpretation of the "
                    f"Steenrod barcode is not theoretically guaranteed.  "
                    f"Use absolute=False for the unambiguous relative "
                    f"convention, or pass a smaller max_dim to limit the "
                    f"check to the dimensions you care about.",
                    UserWarning,
                    stacklevel=2,
                )
        else:
            ordinary, steenrod = iface._compute(int(n_jobs))

        ord_fin, ord_inf = ordinary
        st_fin,  st_inf  = steenrod

        # Sq^0 is the identity, so the Steenrod barcode equals the ordinary
        # one.  The C++ pipeline (Steenrod_barcode.h) already short-circuits
        # k == 0; this Python branch is defensive and also gives the user
        # genuinely distinct arrays to mutate if they want to.
        if int(k) == 0:
            st_fin, st_inf = _copy_per_dim(ord_fin, ord_inf)

        # Relative path: C++ returns every dimension.  Truncate here.
        # Absolute path: C++ already truncated (and used the truncation for
        # the duality check), so skip the slice.
        if max_dim is not None and not absolute:
            ord_fin, ord_inf = _truncate_per_dim(ord_fin, ord_inf, int(max_dim))
            st_fin,  st_inf  = _truncate_per_dim(st_fin,  st_inf,  int(max_dim))

        return (ord_fin, ord_inf), (st_fin, st_inf)

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

    def collapse_edges(self, nb_iterations=1) -> SimplexTree:
        """Assuming the complex is a graph (simplices of higher dimension are ignored), this method implicitly
        interprets it as the 1-skeleton of a flag complex, and replaces it with another (smaller) graph whose
        expansion has the same persistent homology, using a technique known as edge collapses
        (see :cite:`edgecollapsearxiv`).

        A natural application is to get a simplex tree of dimension 1 from :class:`~gudhi.RipsComplex`,
        then collapse edges, perform :meth:`expansion()` and finally compute persistence
        (cf. :download:`rips_complex_edge_collapse_example.py <../example/rips_complex_edge_collapse_example.py>`).

        :param nb_iterations: The number of edge collapse iterations to perform. Default is 1.
        :type nb_iterations: int

        .. warning::
            The current simplex tree is assumed to be a graph, that is of maximal dimension 1.
            If it is not the case, all higher dimensional simplices will get lost during the
            reduction process and not be reinserted. To regain them, call `expansion(max_dim)` afterwards.
        """
        if nb_iterations < 1:
            return
        if self.dimension() > 1:
            message = "collapse_edges() ignores all the simplices of dimension 2 or more in this complex."
            # Always print this specific warning
            warnings.filterwarnings("always", category=RuntimeWarning, message=message)
            warnings.warn(message, RuntimeWarning)
        super()._collapse_edges(nb_iterations)
        return self
