# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - 2025/03 Hannah Schreiber: Use nanobind instead of Cython for python bindings.
#   - YYYY/MM Author: Description of the modification

__license__ = "MIT"


import numpy as np

from gudhi._cubical_complex_ext import (
    _Periodic_cubical_complex_interface,
    _Periodic_cubical_complex_persistence_interface,
)


# PeriodicCubicalComplex python interface
class PeriodicCubicalComplex(_Periodic_cubical_complex_interface):
    """The PeriodicCubicalComplex is an example of a structured complex useful
    in computational mathematics (specially rigorous numerics) and image
    analysis.
    """

    def __init__(
        self,
        dimensions=None,
        top_dimensional_cells=None,
        vertices=None,
        periodic_dimensions=None,
        perseus_file="",
    ):
        """PeriodicCubicalComplex constructor from dimensions and cells
        (either vertices or top dimensional) or from a Perseus-style file name.

        Note that in case the periodic cubical complex is constructed from
        'top_dimensional_cells' or 'vertices' using a flat array,
        the given 'dimensions' will be considered using the
        fortran ordering (column-major order).

        :param dimensions: A list of number of top dimensional cells.
        :type dimensions: list of int
        :param top_dimensional_cells: A list of top dimensional cells filtration values.
        :type top_dimensional_cells: list of double
        :param periodic_dimensions: A list of cells periodicity value.
        :type periodic_dimensions: list of boolean

        Or

        :param dimensions: A list of number of vertices.
        :type dimensions: list of int
        :param vertices: A list of vertices filtration values.
        :type vertices: list of double
        :param periodic_dimensions: A list of cells periodicity value.
        :type periodic_dimensions: list of boolean

        Or

        :param top_dimensional_cells: A multidimensional array of top
            dimensional cells filtration values.
        :type top_dimensional_cells: anything convertible to a numpy ndarray
        :param periodic_dimensions: A list of cells periodicity value.
        :type periodic_dimensions: list of boolean

        Or

        :param vertices: A multidimensional array vertices filtration values.
        :type vertices: anything convertible to a numpy ndarray
        :param periodic_dimensions: A list of cells periodicity value.
        :type periodic_dimensions: list of boolean

        Or

        :param perseus_file: A `Perseus-style <fileformats.html#perseus>`_ file name, giving the filtration values
            of top-dimensional cells and the periodicity.
        :type perseus_file: string
        """
        self._pers = None
        self._built_from_vertices = False
        if perseus_file:
            if (
                top_dimensional_cells is not None
                or vertices is not None
                or dimensions is not None
                or periodic_dimensions is not None
            ):
                raise ValueError(
                    "The Perseus file contains all the information, do not specify anything else"
                )
            super().__init__(perseus_file)
            return
        if periodic_dimensions is None:
            raise ValueError("The periodic_dimensions must be specified")
        if top_dimensional_cells is not None:
            if vertices is not None:
                raise ValueError(
                    "Can only specify the top dimensional cells OR the vertices, not both"
                )
            array = top_dimensional_cells
        elif vertices is not None:
            array = vertices
            self._built_from_vertices = True
        else:
            raise ValueError(
                "Must specify one of top_dimensional_cells, vertices, or perseus_file"
            )
        if dimensions is None:
            array = np.asarray(array, order="F")
            dimensions = array.shape
            array = array.ravel(order="F")
        super().__init__(dimensions, array, periodic_dimensions, vertices is None)

    def _is_persistence_defined(self):
        """Returns true if Persistence pointer is not None."""
        return self._pers != None

    def all_cells(self):
        """Array with the filtration values of all the cells of the complex.
        Modifying the values is strongly discouraged.

        :returns:  numpy.ndarray
        """
        a = super()._get_numpy_array()
        return a.reshape(
            [2 * d + (not p) for (d, p) in zip(super().shape(), super().periodicities())],
            order="F",
        )

    def top_dimensional_cells(self):
        """Array with the filtration values of the top-dimensional cells of the complex.
        Modifying the values is strongly discouraged.

        :returns:  numpy.ndarray
        """
        return self.all_cells()[(slice(1, None, 2),) * super().dimension()]

    def vertices(self):
        """Array with the filtration values of the vertices of the complex.
        Modifying the values is strongly discouraged.

        :returns:  numpy.ndarray
        """
        return self.all_cells()[(slice(0, None, 2),) * super().dimension()]

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
        self._pers = _Periodic_cubical_complex_persistence_interface(self, True)
        self._pers.compute_persistence(homology_coeff_field, min_persistence)

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
        return self._pers.get_persistence()

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
        Therefore, using cofaces_of_persistence_pairs with a cubical complex constructed from vertices
        can lead to an undefined behavior.

        This function returns a list of pairs of top-dimensional cells corresponding to
        the persistence birth and death cells of the filtration. The cells are represented by
        their indices in the input list of top-dimensional cells (and not their indices in the
        internal data structure that includes non-maximal cells). Note that when two adjacent
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
        if self._pers == None:
            raise RuntimeError(
                "compute_persistence() must be called before cofaces_of_persistence_pairs()"
            )
        if self._built_from_vertices:
            raise RuntimeError(
                "cofaces_of_persistence_pairs() only makes sense for a complex"
                " initialized from the values of the top-dimensional cells"
            )

        output = [[], []]
        # TODO: verify the return type of cofaces_of_cubical_persistence_pairs() by nanobind
        # a copy is perhaps avoidable?
        pr = np.array(self._pers.cofaces_of_cubical_persistence_pairs())

        ess_ind = np.argwhere(pr[:, 2] == -1)[:, 0]
        ess = pr[ess_ind]
        max_h = np.max(ess[:, 0]) + 1 if len(ess) > 0 else 0
        for h in range(max_h):
            hidxs = np.argwhere(ess[:, 0] == h)[:, 0]
            output[1].append(ess[hidxs][:, 1])

        reg_ind = np.setdiff1d(np.arange(len(pr)), ess_ind)
        reg = pr[reg_ind]
        max_h = np.max(reg[:, 0]) + 1 if len(reg) > 0 else 0
        for h in range(max_h):
            hidxs = np.argwhere(reg[:, 0] == h)[:, 0]
            output[0].append(reg[hidxs][:, 1:])
        return output

    def vertices_of_persistence_pairs(self):
        """This function returns a list of pairs of vertices corresponding to
        the persistence birth and death cells of the filtration. The cells are represented by
        their indices in the input list of vertices (and not their indices in the
        internal data structure that includes non-minimal cells). Note that when two adjacent
        vertices have the same filtration value, we arbitrarily return one of the two
        when calling the function on one of their common faces.

        Since the cubical complex construction from vertices is different from the top dimensional one,
        the former may not have an equivalent with the second construction, and vice versa.
        Therefore, using vertices_of_persistence_pairs with a cubical complex constructed from top_dimensional_cells
        can lead to an undefined behavior.

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
        if self._pers == None:
            raise RuntimeError(
                "compute_persistence() must be called before vertices_of_persistence_pairs()"
            )
        if not self._built_from_vertices:
            raise RuntimeError(
                "vertices_of_persistence_pairs() only makes sense for a complex"
                " initialized from the values of the vertices"
            )

        output = [[], []]
        # TODO: verify the return type of cofaces_of_cubical_persistence_pairs() by nanobind
        # a copy is perhaps avoidable?
        pr = np.array(self._pers.vertices_of_cubical_persistence_pairs())

        # except pr, is there any difference from the method above? If not, factorize?
        ess_ind = np.argwhere(pr[:, 2] == -1)[:, 0]
        ess = pr[ess_ind]
        max_h = max(ess[:, 0]) + 1 if len(ess) > 0 else 0
        for h in range(max_h):
            hidxs = np.argwhere(ess[:, 0] == h)[:, 0]
            output[1].append(ess[hidxs][:, 1])

        reg_ind = np.setdiff1d(np.array(range(len(pr))), ess_ind)
        reg = pr[reg_ind]
        max_h = max(reg[:, 0]) + 1 if len(reg) > 0 else 0
        for h in range(max_h):
            hidxs = np.argwhere(reg[:, 0] == h)[:, 0]
            output[0].append(reg[hidxs][:, 1:])
        return output

    def betti_numbers(self):
        """This function returns the Betti numbers of the complex.

        :returns: list of int -- The Betti numbers ([B0, B1, ..., Bn]).

        :note: betti_numbers function requires :func:`compute_persistence` function to be
            launched first.

        :note: This function always returns the Betti numbers of a torus as infinity
            filtration cubes are not removed from the complex.
        """
        if self._pers == None:
            raise RuntimeError("compute_persistence() must be called before betti_numbers()")
        return self._pers.betti_numbers()

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
        if self._pers == None:
            raise RuntimeError(
                "compute_persistence() must be called before persistent_betti_numbers()"
            )
        return self._pers.persistent_betti_numbers(from_value, to_value)

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
        if self._pers == None:
            raise RuntimeError(
                "compute_persistence() must be called before persistence_intervals_in_dimension()"
            )
        piid = np.array(self._pers.intervals_in_dimension(dimension))
        # Workaround https://github.com/GUDHI/gudhi-devel/issues/507
        if len(piid) == 0:
            return np.empty(shape=[0, 2])
        return piid
