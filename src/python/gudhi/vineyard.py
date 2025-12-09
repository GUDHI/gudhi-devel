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
from numpy.typing import ArrayLike
import warnings
import matplotlib as mpl
import matplotlib.pyplot as plt

from gudhi import _vineyard_ext as t
from gudhi.rips_complex import RipsComplex
from gudhi.reader_utils import read_lower_triangular_matrix_from_csv_file


class Vineyard(t.Vineyard_interface):
    # TODO: representative cycles

    def __init__(self, store_latest_cycles: bool = False, cycles_dim: int = -1):
        super().__init__(store_latest_cycles, cycles_dim)

    def initialize(
        self,
        filtered_cpx=None,
        boundaries: list[np.ndarray] = None,
        dimensions: np.ndarray = None,
        filtration_values: np.ndarray = None,
        number_of_updates: int = 0,
    ) -> list[np.ndarray]:
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

    def _gray_on_diagonal(self, ax, x, y, z):
        gray = (0.75, 0.75, 0.75)
        mask = np.concatenate(([0], np.asarray(x == y), [0]))
        diff = np.diff(mask)
        starts = np.where(diff == 1)[0]
        ends = np.where(diff == -1)[0]
        for s, e in zip(starts, ends):
            ax.plot3D(x[s:e], y[s:e], z[s:e], c=gray)

    def _plot_vines(
        self, ax, vines, inf_v, min_bar_length, gray_diagonal, c=None, c_step=None
    ):
        for i in range(vines.shape[0]):
            x = vines[i, :, 0]
            y = vines[i, :, 1]
            y[y >= np.inf] = inf_v
            x[x >= np.inf] = inf_v
            z = np.asarray(range(vines.shape[1]))
            diff = y - x
            diff[diff < min_bar_length] = 0
            if np.any(diff):
                ax.plot3D(x, y, z, c=c)
                if gray_diagonal:
                    self._gray_on_diagonal(ax, x, y, z)
                if c_step is not None:
                    c = c - c_step

    def plot_vineyards(
        self,
        dim: int = None,
        max_dim: int = None,
        min_bar_length: np.number = -1,
        gray_diagonal: bool = True,
    ):
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")

        if dim is None:
            vineyard = self.get_current_vineyard_view()
            if max_dim is not None:
                vineyard = vineyard[: max_dim + 1]
            max_death = max(
                [
                    np.max(vines[:, :, 1], where=~np.isinf(vines[:, :, 1]), initial=-1)
                    for vines in vineyard
                ]
            )
            # arbitrary, just to distance them a bit from the finite bars
            inf_v = max_death + 10

            for d, vines in enumerate(vineyard):
                # very naive and non robust way to manage colors...
                c = np.random.rand(3)
                step = np.min(c / vines.shape[1])
                self._plot_vines(ax, vines, inf_v, min_bar_length, gray_diagonal, c, step)
        else:
            vines = self.get_current_vineyard_view(dim=dim)
            max_death = np.max(vines[:, :, 1], where=~np.isinf(vines[:, :, 1]), initial=-1)
            # arbitrary, just to distance them a bit from the finite bars
            inf_v = max_death + 10

            self._plot_vines(ax, vines, inf_v, min_bar_length, gray_diagonal)

        plt.show()
        return


class PointCloudRipsVineyard:
    # TODO: store representative cycles

    def __init__(self, store_point_coordinates: bool, store_cycles: bool):
        self._vineyard = Vineyard(store_cycles, cycles_dim=1)
        self._store_points = store_point_coordinates
        if store_point_coordinates:
            self._points = []
        self._store_cycles = store_cycles
        if store_cycles:
            self._cycles = []

    @classmethod
    def from_files(
        cls,
        path_prefix: str,
        path_suffix: str = ".txt",
        first_index: int = 0,
        delimiter: str = None,
        number_of_updates: int = None,
        file_type: str = "point_cloud",
        store_point_coordinates: bool = False,
        store_cycles: bool = False,
    ):
        res = cls(store_point_coordinates, store_cycles)

        # assumes counting starts with 0
        # will fail with "file not found" if not
        path = os.path.realpath(path_prefix + "0" + path_suffix, strict=True)
        path_prefix = path[: -len("0" + path_suffix)]

        if number_of_updates is None:
            number_of_updates = (
                len(glob.glob(path_prefix + "*" + path_suffix)) - first_index - 1
            )

        res.initialize(
            path=path_prefix + str(first_index) + path_suffix,
            delimiter=delimiter,
            number_of_updates=number_of_updates,
            data_type=file_type,
        )

        for step in range(first_index + 1, number_of_updates + 1):
            res.update(
                path=path_prefix + str(step) + path_suffix,
                delimiter=delimiter,
                data_type=file_type,
            )

        return res

    @classmethod
    def from_tensors(
        cls,
        data: ArrayLike,
        first_index: int = 0,
        data_type: str = "point_cloud",
        store_point_coordinates: bool = False,
        store_cycles: bool = False,
    ):
        res = cls(store_point_coordinates, store_cycles)

        if number_of_updates is None or number_of_updates > len(data) - first_index - 1:
            number_of_updates = len(data) - first_index - 1

        res.initialize(
            data=data[first_index],
            number_of_updates=number_of_updates,
            data_type=data_type,
        )

        for step in range(first_index + 1, number_of_updates + 1):
            res.update(
                data=data[step],
                data_type=data_type,
            )

        return res

    def _create_simplex_tree_from_file(self, path: str, delimiter: str, file_type: str):
        if file_type == "point_cloud":
            points = np.loadtxt(path, delimiter=delimiter)
            if self._store_points:
                self._points.append(points)
            return RipsComplex(points=points).create_simplex_tree(max_dimension=2)
        elif file_type == "distance_matrix":
            matrix = read_lower_triangular_matrix_from_csv_file(
                path_prefix + str(first_index) + path_suffix, separator=delimiter
            )
            return RipsComplex(distance_matrix=matrix).create_simplex_tree(max_dimension=2)
        else:
            raise ValueError(
                "'file_type' has to be either 'point_cloud' or 'distance_matrix'."
            )

    def _create_simplex_tree_from_data(self, data: ArrayLike, data_type: str):
        if data_type == "point_cloud":
            if self._store_points:
                self._points.append(data)
            return RipsComplex(points=data).create_simplex_tree(max_dimension=2)
        elif data_type == "distance_matrix":
            return RipsComplex(distance_matrix=data).create_simplex_tree(max_dimension=2)
        else:
            raise ValueError(
                "'data_type' has to be either 'point_cloud' or 'distance_matrix'."
            )

    def initialize(
        self,
        data: ArrayLike = None,
        path: str = None,
        delimiter: str = None,
        number_of_updates: int = 0,
        data_type: str = "point_cloud",
    ):
        if self._store_points:
            if data_type == "distance_matrix":
                warnings.warn(
                    "If `data_type` is `distance_matrix`, the point coordinates cannot be stored.",
                    UserWarning,
                )
            self._points = []

        if data is None:
            if path is None:
                raise ValueError("Either `data` or `path` have to be specified.")
            self._complex = self._create_simplex_tree_from_file(path, delimiter, data_type)
        elif path is None:
            self._complex = self._create_simplex_tree_from_data(data, data_type)
        else:
            raise ValueError("Both `data` and `path` cannot be specified.")

        self._vineyard.initialize(
            self._complex,
            number_of_updates=number_of_updates,
        )

        if self._store_cycles:
            self._cycles = []
            self._cycles.append(self._vineyard.get_latest_representative_cycles())

    def update(
        self,
        data: ArrayLike = None,
        path: str = None,
        delimiter: str = None,
        data_type: str = "point_cloud",
    ):
        if self._store_points and data_type == "distance_matrix":
            warnings.warn(
                "If `data_type` is `distance_matrix`, the point coordinates cannot be stored.",
                UserWarning,
            )

        if data is None:
            if path is None:
                raise ValueError("Either `data` or `path` have to be specified.")
            self._vineyard.update(
                self._create_simplex_tree_from_file(path, delimiter, data_type)
            )
        elif path is None:
            self._vineyard.update(self._create_simplex_tree_from_data(data, data_type))
        else:
            raise ValueError("Both `data` and `path` cannot be specified.")

        if self._store_cycles:
            self._cycles.append(self._vineyard.get_latest_representative_cycles())

    def get_current_vineyard_view(self, dim: int = None) -> list[np.ndarray] | np.ndarray:
        return self._vineyard.get_current_vineyard_view(dim)

    def get_complex(self) -> tuple[list[np.ndarray], np.ndarray]:
        return t._build_boundary_matrix_from_complex(self._complex)

    def get_points(self, step: int = None) -> np.ndarray:
        if not self._store_points:
            raise NotImplementedError(
                "Points cannot be retrieved if the store options is at False."
            )

        if step is None:
            return np.asarray(self._points)
        return self._points[step]

    def get_1D_representative_cycles(
        self, step: int = None
    ) -> list[tuple[list[np.ndarray], np.number]] | tuple[list[np.ndarray], np.number]:
        if not self._store_cycles:
            raise NotImplementedError(
                "Cycles cannot be retrieved if the store options is at False."
            )

        if step is None:
            return self._cycles
        return self._cycles[step]

    def plot_vineyards(
        self, dim: int = None, max_dim: int = None, min_bar_length: np.number = -1
    ):
        self._vineyard.plot_vineyards(dim, max_dim, min_bar_length)

    def _plot_cycle(self, axes, cycle, points, cpx):
        c = np.random.rand(3)
        for u, v in [(points[cpx[idx][0]], points[cpx[idx][1]]) for idx in cycle]:
            if points.shape[1] == 2:
                axes.plot([u[0], v[0]], [u[1], v[1]], color=c)
            else:
                axes.plot([u[0], v[0]], [u[1], v[1]], [u[2], v[2]], color=c)

    def plot_1D_representative_cycles(
        self, step: int, index: int = None, min_bar_length: np.number = -1
    ):
        points = self.get_points(step)
        cycles = self.get_1D_representative_cycles(step)
        cpx, _ = self.get_complex()

        fig = plt.figure()
        if points.shape[1] == 2:
            axes = fig.add_subplot()
        elif points.shape[1] == 3:
            axes = fig.add_subplot(projection="3d")
        else:
            raise ValueError("Plotting only possible in 2D and 3D.")

        if index is None:
            for cycle in cycles:
                if cycle[1] >= min_bar_length:
                    self._plot_cycle(axes, cycle[0], points, cpx)
        else:
            if min_bar_length != -1:
                warnings.warn(
                    "Specified argument `min_length` is ignored.",
                    UserWarning,
                )
            self._plot_cycle(axes, cycles[index][0], points, cpx)

        x = points[:, 0]
        y = points[:, 1]
        if points.shape[1] == 2:
            axes.plot(x, y, linestyle="none", markersize=3, marker="o")
        else:
            axes.plot(x, y, points[:, 2], linestyle="none", markersize=3, marker="o")

        plt.show()
        return
