# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Hannah Schreiber
#
# Copyright (C) 2025 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

__license__ = "MIT"


import os
import glob
import numpy as np
import warnings
import matplotlib as mpl
import matplotlib.pyplot as plt

from gudhi import _vineyard_ext as t
from gudhi.rips_complex import RipsComplex


class Vineyard(t.Vineyard_interface):
    """The data structure is a one skeleton graph, or Rips graph, containing edges when the edge length is less or
    equal to a given threshold. Edge length is computed from a user given point cloud with a given distance function,
    or a distance matrix.
    """

    def __init__(self):
        """RipsComplex constructor.

        :param points: A list of points in d-Dimension.
        :type points: Sequence[Sequence[float]] or any array like object of ndim 2 and dtype convertible to float.

        Or

        :param distance_matrix: A distance matrix (full square or lower triangular).
        :type distance_matrix: Sequence[Sequence[float]] (square or just the lower triangle) or any square array like
            object of dtype convertible to float.

        And in both cases

        :param max_edge_length: Maximal edge length. All edges of the graph strictly greater than `threshold` are not
            inserted in the graph.
        :type max_edge_length: float
        :param sparse: If this is not None, it switches to building a sparse Rips and represents the approximation
            parameter epsilon.
        :type sparse: float
        """
        super().__init__()

    def build_rips_vines_from_point_cloud_files(
        self,
        path_prefix: str,
        path_suffix: str = ".txt",
        number_of_updates: int = None,
        max_edge_length: float = float("inf"),
    ) -> list[np.ndarray]:
        """
        :param max_dimension: graph expansion for Rips until this given maximal dimension.
        :type max_dimension: int
        :returns: A simplex tree encoding the Vietoris–Rips filtration.
        :rtype: SimplexTree
        """

        if number_of_updates is None:
            # path_prefix should be absolute to avoid any problems
            number_of_updates = len(glob.glob(path_prefix + "*" + path_suffix))

        # assumes file counting starts with 0
        # TODO: allow to start with any number
        points = np.loadtxt(path_prefix + "0" + path_suffix)
        # careful with max_edge_length
        rips_complex = RipsComplex(points=points, max_edge_length=max_edge_length)
        self.initialize(
            rips_complex.create_simplex_tree(max_dimension=2),
            number_of_updates=number_of_updates,
        )

        for step in range(1, number_of_updates):
            points = np.loadtxt(path_prefix + str(step) + path_suffix)
            rips_complex = RipsComplex(points=points, max_edge_length=max_edge_length)
            self.update(rips_complex.create_simplex_tree(max_dimension=2))

        return self.get_current_vineyard_view()

    def initialize(
        self,
        filtered_cpx=None,
        boundaries: list[np.ndarray] = None,
        dimensions: np.ndarray = None,
        filtration_values: np.ndarray = None,
        number_of_updates: int = 0,
    ) -> list[np.ndarray]:
        """
        :param max_dimension: graph expansion for Rips until this given maximal dimension.
        :type max_dimension: int
        :returns: A simplex tree encoding the Vietoris–Rips filtration.
        :rtype: SimplexTree
        """
        if filtered_cpx is not None:
            if boundaries != None or dimensions != None or filtration_values != None:
                raise ValueError(
                    "Either a filtered complex or (boundaries/dimensions/filtration values) should be given, but not both."
                )
            super()._initialize_from_complex(filtered_cpx, number_of_updates)
            return self.get_current_vineyard_view()

        if boundaries is None or dimensions is None or filtration_values is None:
            raise ValueError(
                "Either a filtered complex or (boundaries/dimensions/filtration values) should be not None."
            )
        super()._initialize(boundaries, dimensions, filtration_values, number_of_updates)
        return self.get_current_vineyard_view()

    def update(
        self,
        filtered_cpx=None,
        filtration_values: np.ndarray = None,
    ) -> list[np.ndarray]:
        """
        :param max_dimension: graph expansion for Rips until this given maximal dimension.
        :type max_dimension: int
        :returns: A simplex tree encoding the Vietoris–Rips filtration.
        :rtype: SimplexTree
        """
        if filtered_cpx is not None:
            if filtration_values != None:
                raise ValueError(
                    "Either a filtered complex or new filtration values should be given, but not both."
                )
            super()._update_from_complex(filtered_cpx)
            return self.get_current_vineyard_view()

        if filtration_values is None:
            raise ValueError(
                "Either a filtered complex or filtration values should be not None."
            )
        super()._update(filtration_values)
        return self.get_current_vineyard_view()

    def get_current_vineyard_view(self, dim: int = None) -> list[np.ndarray] | np.ndarray:
        # "view" because everything is read only and should not trigger a copy
        vineyard = super()._get_current_vineyard_view()
        # shape is update number x vine number x (birth, death)
        # so we change it to vine number x update number x (birth, death)
        # which seems more intuitive
        if dim is None:
            return [np.swapaxes(v, 0, 1) for v in vineyard]

        if dim >= len(vineyard):
            warnings.warn("No vine of given dimension was computed.", UserWarning)
            return np.empty((0, 0, 2))

        return np.swapaxes(vineyard[dim], 0, 1)

    def plot_vineyards(self, dim: int = None, max_dim: int = None):
        fig = plt.figure()
        ax = plt.axes(projection="3d")

        if dim is None:
            vineyard = self.get_current_vineyard_view()
            if max_dim is not None:
                vineyard = vineyard[:max_dim+1]
            max_death = max([np.max(vines[:,:,1], where=~np.isinf(vines[:,:,1]), initial=-1) for vines in vineyard])

            for d, vines in enumerate(vineyard):
                c = np.random.rand(3)   # very naive and non robust way to manage colors...
                num_vines = vines.shape[0]
                step = np.min(c / num_vines)
                for i in range(num_vines):
                    x = vines[i, :, 0]
                    y = vines[i, :, 1]
                    y[y >= np.inf] = max_death + 10     # arbitrary, just to distance them a bit from the finite bars
                    z = np.asarray(range(vines.shape[1]))
                    if np.any(y - x):
                        ax.plot3D(x, y, z, c = c)
                        c = c - step
        else:
            vines = self.get_current_vineyard_view(dim = dim)
            max_death = np.max(vines[:,:,1], where=~np.isinf(vines[:,:,1]), initial=-1)

            for i in range(vines.shape[0]):
                x = vines[i, :, 0]
                y = vines[i, :, 1]
                y[y >= np.inf] = max_death * 10     # arbitrary, just to distance them a bit from the finite bars
                z = np.asarray(range(vines.shape[1]))
                if np.any(y - x):
                    ax.plot3D(x, y, z)

        plt.show()
        return
