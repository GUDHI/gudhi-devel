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
from typing import Literal
import warnings
import matplotlib as mpl
import matplotlib.pyplot as plt

from gudhi import _vineyard_ext as t
from gudhi.rips_complex import RipsComplex
from gudhi.reader_utils import read_lower_triangular_matrix_from_csv_file


# for debug only
def _verify_validity(
    vineyard: list[np.ndarray],
    path_prefix: str,
    path_suffix: str = ".txt",
    first_index: int = 0,
    delimiter: str = None,
    number_of_updates: int = None,
) -> bool:
    # assumes counting starts with 0
    # will fail with "file not found" if not
    path = os.path.realpath(path_prefix + "0" + path_suffix, strict=True)
    path_prefix = path[: -len("0" + path_suffix)]

    if number_of_updates is None:
        number_of_updates = len(glob.glob(path_prefix + "*" + path_suffix)) - first_index - 1

    if len(vineyard) == 0:
        print("Empty vineyard. There should be at least one connected component.")
        return False

    if vineyard[0].shape[1] != number_of_updates + 1:
        print(
            "Number of updates does not correspond.",
            vineyard[0].shape[1],
            number_of_updates + 1,
        )
        return False

    print("Comparing step:", end=" ")
    for step in range(0, number_of_updates + 1):
        print(step, end=" ")
        points = np.loadtxt(
            path_prefix + str(step + first_index) + path_suffix, delimiter=delimiter
        )
        st = RipsComplex(points=points).create_simplex_tree(max_dimension=2)
        pers = st.persistence(homology_coeff_field=2, min_persistence=-1)
        pers = [bar[1] for bar in pers if bar[0] == 1]
        if len(vineyard) == 1 and len(pers) != 0:
            print("\nNo 1-dimensional bars in vineyard while there should be some.", pers)
            return False
        v_pers = [(b, d) for b, d in vineyard[1][:, step, :]]
        if len(pers) != len(v_pers):
            print(
                "\nThere are not as many bars in the vineyard then there should be.",
                len(v_pers),
                len(pers),
            )
            return False
        pers.sort()
        v_pers.sort()
        if pers != v_pers:
            pers = np.asarray(pers)
            v_pers = np.asarray(v_pers)
            mask = ~(pers == v_pers)
            mask = np.any(mask, axis=1)
            print()
            print(pers[mask])
            print(v_pers[mask])
            return False

    print("\nBarcodes are valid.")
    return True


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

    def _denoise_vineyard(self, vineyard: np.ndarray, min_bar_length: np.number):
        mask = (vineyard[:, :, 1] - vineyard[:, :, 0]) >= min_bar_length
        mask = np.any(mask, axis=1)
        return vineyard[mask]

    def get_current_vineyard(
        self, dim: int = None, min_bar_length: np.number = -1
    ) -> list[np.ndarray] | np.ndarray:
        vineyard = self.get_current_vineyard_view(dim=dim)

        if dim is None:
            return [self._denoise_vineyard(v, min_bar_length) for v in vineyard]

        return self._denoise_vineyard(vineyard, min_bar_length)

    def _gray_on_band(self, ax, x, y, z, band):
        gray = (0.75, 0.75, 0.75)
        mask = np.concatenate(([0], np.asarray(y - x <= band), [0]))
        diff = np.diff(mask)
        starts = np.where(diff == 1)[0]
        ends = np.where(diff == -1)[0]
        for s, e in zip(starts, ends):
            ax.plot3D(x[s:e], y[s:e], z[s:e], c=gray)

    def _erase_on_band(self, ax, x, y, z, band, c):
        mask = np.concatenate(([0], np.asarray(y - x > band), [0]))
        diff = np.diff(mask)
        starts = np.where(diff == 1)[0]
        ends = np.where(diff == -1)[0]
        for s, e in zip(starts, ends):
            ax.plot3D(x[s:e], y[s:e], z[s:e], c=c)

    def _get_colors_for_all_dim(self, dim: int):
        plot_colors = [
            "Oranges",
            "Purples",
            "Greens",
            "Blues",
            "Reds",
            "YlOrBr",
            "PuRd",
            "YlOrRd",
            "OrRd",
            "RdPu",
            "BuPu",
            "GnBu",
            "PuBu",
            "YlGnBu",
            "PuBuGn",
            "BuGn",
            "YlGn",
        ]
        return mpl.colormaps[plot_colors[dim % len(plot_colors)]]

    def _plot_vines(
        self,
        ax,
        vines,
        inf_v,
        min_bar_length,
        noise_option: Literal[
            "none", "gray_diagonal", "gray_band", "erase_diagonal", "erase_band"
        ],
        cmap,
    ):
        real_vines = []
        for i in range(vines.shape[0]):
            x = vines[i, :, 0]
            y = vines[i, :, 1]
            y[y >= np.inf] = inf_v
            x[x >= np.inf] = inf_v
            z = np.asarray(range(vines.shape[1]))
            diff = y - x
            diff[diff < min_bar_length] = 0
            if np.any(diff):
                real_vines.append((x, y, z))

        cmap = list(reversed(cmap(np.linspace(0, 1, len(real_vines) + 2))))

        for i, (x, y, z) in enumerate(real_vines):
            match noise_option:
                case "none":
                    ax.plot3D(x, y, z, c=cmap[i])
                case "gray_diagonal":
                    ax.plot3D(x, y, z, c=cmap[i])
                    self._gray_on_band(ax, x, y, z, 0)
                case "gray_band":
                    ax.plot3D(x, y, z, c=cmap[i])
                    self._gray_on_band(ax, x, y, z, min_bar_length)
                case "erase_diagonal":
                    self._erase_on_band(ax, x, y, z, 0, cmap[i])
                case "erase_band":
                    self._erase_on_band(ax, x, y, z, min_bar_length, cmap[i])
                case _:
                    raise ValueError(
                        "argument `noise_option` does not contain a valid literal: "
                        + noise_option
                        + ". Possibilities are: `none`, `gray_diagonal`, `gray_band`, `erase_diagonal` and `erase_band`"
                    )

    def plot_vineyards(
        self,
        dim: int = None,
        max_dim: int = None,
        min_bar_length: np.number = -1,
        noise_option: Literal[
            "none", "gray_diagonal", "gray_band", "erase_diagonal", "erase_band"
        ] = "gray_diagonal",
        square_scaling: bool = True,
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
                cmap = self._get_colors_for_all_dim(d)
                self._plot_vines(ax, vines, inf_v, min_bar_length, noise_option, cmap)
        else:
            vines = self.get_current_vineyard_view(dim=dim)
            max_death = np.max(vines[:, :, 1], where=~np.isinf(vines[:, :, 1]), initial=-1)
            # arbitrary, just to distance them a bit from the finite bars
            inf_v = max_death + 10
            self._plot_vines(
                ax, vines, inf_v, min_bar_length, noise_option, mpl.colormaps["viridis"]
            )

        if square_scaling:
            scaling = np.array([getattr(ax, "get_{}lim".format(dim))() for dim in "xy"])
            ax.auto_scale_xyz(
                *[[np.min(scaling), np.max(scaling)]] * 2,
                [getattr(ax, "get_{}lim".format("z"))()],
            )
        ax.set_aspect("equalxy", "box")

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

        for step in range(first_index + 1, first_index + number_of_updates + 1):
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

    def get_current_vineyard(
        self, dim: int = None, min_bar_length: np.number = -1
    ) -> list[np.ndarray] | np.ndarray:
        return self._vineyard.get_current_vineyard(dim, min_bar_length)

    def _get_vertices(self, cpx: list[np.ndarray], i: int) -> set[int]:
        if cpx[i].shape[0] == 0:
            return {i}
        # not optimal for high dimensions, but enables no restrictions on cpx order
        return {v for idx in cpx[i] for v in self._get_vertices(cpx, idx)}

    def get_complex(self, as_vertices: bool = False) -> tuple[list[np.ndarray], np.ndarray]:
        cpx, dims = t._build_boundary_matrix_from_complex(self._complex)

        if not as_vertices:
            return cpx, dims
        cpx = [
            np.asarray(sorted([v for v in self._get_vertices(cpx, i)]))
            for i, _ in enumerate(cpx)
        ]
        return cpx, dims

    def get_points(self, step: int = None) -> np.ndarray:
        if not self._store_points:
            raise NotImplementedError(
                "Points cannot be retrieved if the store options is at False."
            )

        if step is None:
            return np.asarray(self._points)
        return self._points[step]

    def _denoise_cycle(
        self, cycle: list[tuple[np.ndarray, np.number]], min_bar_length: np.number
    ) -> list[tuple[np.ndarray, np.number]]:
        return [c for c in cycle if c[1] >= min_bar_length]

    def get_1D_representative_cycles(
        self, step: int = None, min_bar_length: np.number = 0
    ) -> list[list[tuple[np.ndarray, np.number]]] | list[tuple[np.ndarray, np.number]]:
        if not self._store_cycles:
            raise NotImplementedError(
                "Cycles cannot be retrieved if the store options is at False."
            )

        if min_bar_length <= 0:
            if step is None:
                return self._cycles
            return self._cycles[step]

        if step is None:
            return [self._denoise_cycle(cycle, min_bar_length) for cycle in self._cycles]
        return self._denoise_cycle(self._cycles[step], min_bar_length)

    def plot_vineyards(
        self,
        dim: int = None,
        max_dim: int = None,
        min_bar_length: np.number = -1,
        noise_option: Literal[
            "none", "gray_diagonal", "gray_band", "erase_diagonal", "erase_band"
        ] = "gray_diagonal",
        square_scaling: bool = True,
    ):
        self._vineyard.plot_vineyards(dim, max_dim, min_bar_length, noise_option, square_scaling)

    def _plot_cycle(self, axes, cycle, points, cpx, c):
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

        cmap_f = mpl.colormaps['tab10']
        cmap = cmap_f(np.arange(10))
        cmap_f = mpl.colormaps['Set2']
        cmap = np.concatenate((cmap, cmap_f(np.arange(8))))
        cmap_f = mpl.colormaps['Dark2']
        cmap = np.concatenate((cmap, cmap_f(np.arange(8))))
        csize = 10 + 8 + 8
        i = 0

        if index is None:
            for cycle in cycles:
                if cycle[1] >= min_bar_length:
                    self._plot_cycle(axes, cycle[0], points, cpx, cmap[i])
                    i = (i + 1) % csize
        else:
            if min_bar_length != -1:
                warnings.warn(
                    "Specified argument `min_length` is ignored.",
                    UserWarning,
                )
            self._plot_cycle(axes, cycles[index][0], points, cpx, cmap[i])

        gray = (0.25, 0.25, 0.25)
        x = points[:, 0]
        y = points[:, 1]
        if points.shape[1] == 2:
            axes.plot(x, y, linestyle="none", markersize=3, marker="o", c=gray)
        else:
            axes.plot(x, y, points[:, 2], linestyle="none", markersize=3, marker="o", c=gray)

        plt.show()
        return
