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
from gudhi.simplex_tree import SimplexTree
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
    """Class constructing and plotting vineyards from given inputs. Can also provide non-trivial representative cycles
    for the latest filtration given.

    A trivial representative cycle is a cycle representing a persistence bar of length 0. That is, the cycle is born
    and killed at the same time.
    """

    def __init__(self, store_latest_cycles: bool = False, cycles_dim: int = -1):
        """Constructs an empty class with given options. To construct the vineyard, call first :meth:`initialize` to
        build the initial persistent diagram and then a sequence of :meth:`update` to add a sequence of layers to
        the vineyard.

        :param store_latest_cycles: If `True`, stores the current non-trivial representative cycles of the current
            complex and enables :meth:`get_latest_representative_cycles`. Defaults to False.
        :type store_latest_cycles: bool, optional
        :param cycles_dim: If `store_latest_cycles` is set to `True`: if set to a positive number, only the cycle of
            dimension `cycles_dim` will be stored, otherwise stores all cycles. If `store_latest_cycles` is set to
            `False`: ignored. Defaults to -1.
        :type cycles_dim: int, optional
        """
        super().__init__(store_latest_cycles, cycles_dim)

    def initialize(
        self,
        filtered_cpx: SimplexTree = None,
        boundaries: list[np.ndarray] = None,
        dimensions: np.ndarray = None,
        filtration_values: np.ndarray = None,
        number_of_updates: int = 0,
    ) -> list[np.ndarray]:
        """Initializes the vineyard with the first persistence barcode. If another vineyard was initialized before,
        it will be completely replaced.

        :param filtered_cpx: A SimplexTree containing all simplices and filtration values of the initializing
            filtration. Alternative to providing `boundaries`, `dimensions` and `filtration_values`, so both should not
            be provided at the same time. Defaults to `None`.
        :type filtered_cpx: SimplexTree or `None`
        :param boundaries: List of the cell boundaries (do not have to be simplicial) of the initializing filtration.
            Alternative to providing `filtered_cpx`, so both should not be provide at the same time. If `boundaries` is
            provided, `dimensions` and `filtration_values` also have to be provided (such that the indices aligns).
            The cells in the boundaries have to be indexed by their own position in the list.
            E.g., `[[], [], [], [0, 1], [0, 2], [1, 2], [3, 4, 5]]`, represents the filtration containing a triangle
            and all its faces. Defaults to `None`.
        :type boundaries: list[np.ndarray] or `None`
        :param dimensions: Array of the dimensions of the cells represented in `boundaries`. Has to be provided if
            `boundaries` is provided, and such that `dimensions[i]` corresponds to the dimension of `boundaries[i]`.
            Defaults to `None`.
        :type dimensions: np.ndarray or `None`
        :param filtration_values: Array of the filtration values of the cells represented in `boundaries`. Has to be
            provided if `boundaries` is provided, and such that `filtration_values[i]` corresponds to the filtration
            value of `boundaries[i]`. Defaults to `None`.
        :type filtration_values: np.ndarray or `None`
        :param number_of_updates: Optional (for optimization purposes). Will allocate memory space for
            `number_of_updates` additional steps in the vineyard after initialization. Defaults to 0.
        :type number_of_updates: int, optional
        :raises ValueError: If both `filtered_cpx` and (`boundaries` or `dimensions` or `filtration_values`) are
            provided or if none of both are provided.
        :return: Current state of the vineyard. See :meth:`get_current_vineyard_view`.
        :rtype: list[read-only np.ndarray]
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
                "Either a filtered complex or (boundaries/dimensions/filtration values) should be not `None`."
            )
        super()._initialize(boundaries, dimensions, filtration_values, number_of_updates)
        return self.get_current_vineyard_view()

    def update(
        self,
        filtered_cpx: SimplexTree = None,
        filtration_values: np.ndarray = None,
    ) -> list[np.ndarray]:
        """Adds a layer to the current vineyard by updating the persistence diagram such that it corresponds to the
        given filtration.

        :param filtered_cpx: A SimplexTree whose simplices are identical (also label wise) to the SimplexTree provided
            to :meth:`initialize`, but with new filtration values corresponding to the new filtration. If none was
            provided for :meth:`initialize`, provide the argument `filtration_values` instead. Defaults to `None`.
        :type filtered_cpx: SimplexTree or `None`
        :param filtration_values: Array of new filtration values such that `filtration_values[i]` corresponds to
            the new value for `boundaries[i]` provided to :meth:`initialize`. If `boundaries` was not provided to
            :meth:`initialize`, provide the argument `filtered_cpx` instead. Defaults to `None`.
        :type filtration_values: np.ndarray or `None`
        :raises ValueError: If both `filtered_cpx` and `filtration_values` are provided or none of both.
        :return: Current state of the vineyard. See :meth:`get_current_vineyard_view`.
        :rtype: list[read-only np.ndarray]
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
                "Either a filtered complex or filtration values should be not `None`."
            )
        super()._update(filtration_values)
        return self.get_current_vineyard_view()

    def get_current_vineyard_view(self, dim: int = None) -> list[np.ndarray] | np.ndarray:
        """Returns the list of read-only and unfiltered vine views. See :meth:`get_current_vineyard` for a more flexible
        output. The format of the list is `dimension x vine number x update number x (birth, death)`, e.g.,
        `vineyard[a][b][c][0]` returns the birth value of the `a`-dimensional vine number `b` at step `c`
        (initialization is step 0), while `vineyard[a][b][c][1]` returns the corresponding death value.

        :param dim: Optional. If provided, the sub-array at index `dim` (corresponding to the vines of dimension `dim`)
            is returned instead of the whole vineyard. Defaults to `None`.
        :type dim: int, optional
        :return: List of read-only views.
        :rtype: list[read-only np.ndarray] (default) or read-only np.ndarray (if `dim` provided)

        .. note::
            As the returned vines are only views, they will get destroyed/invalidated if this class gets destroyed.
        """
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
        """Returns a copy of the current vineyard. If no copy is desired, see :meth:`get_current_vineyard_view`.
        The format of the returned list is `dimension x vine number x update number x (birth, death)`, e.g.,
        `vineyard[a][b][c][0]` returns the birth value of the `a`-dimensional vine number `b` at step `c`
        (initialization is step 0), while `vineyard[a][b][c][1]` returns the corresponding death value.

        :param dim: Optional. If provided, only the vines at given dimension are copied and the output format looses
            the first axes. Defaults to `None`.
        :type dim: int, optional
        :param min_bar_length: Optional. If provided, only vines with at least one coordinate corresponding to a bar of
            length equal or higher than `min_bar_length` are copied. Defaults to -1 (i.e. all vines are copied).
        :type min_bar_length: Any numerical type coercible to the filtration value type, optional
        :return: A copy of the current vineyard.
        :rtype: list[np.ndarray] (default) or np.ndarray (if `dim` provided)
        """
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
        """Plots the current vineyard, except for completely trivial vines (vines where all coordinates are on
        the diagonal).

        :param dim: Optional. If provided, only plots the vines of the given dimension. Defaults to `None`.
        :type dim: int, optional
        :param max_dim: Optional. If provided, only plots the vines with the given maximal dimension included. If `dim`
            was already provided, `max_dim` is ignored. Defaults to `None`.
        :type max_dim: int, optional
        :param min_bar_length: Optional. If provided, only vines with at least one coordinate corresponding to a bar of
            length equal or higher than `min_bar_length` are plotted. Defaults to -1 (i.e., all vines except trivial
            ones are plotted). Note that trivial vines are never plotted even if `min_bar_length` is set to 0.
        :type min_bar_length: Any numerical type coercible to the filtration value type, optional
        :param noise_option: Value describing how "noisy" parts of a vine should be treated. There are 5 options for
            now:
            - "none": Nothing particular is done.
            - "gray_diagonal": Every part of a vine contained in the diagonal will be plotted with a gray overline.
            - "gray_band": Every part of a vine where all coordinates differ by less than `min_bar_length` are plotted
                with a gray overline.
            - "erase_diagonal": Every part of a vine contained in the diagonal are not plotted.
            - "erase_band": Every part of a vine where all coordinates differ by less than `min_bar_length` are not
                plotted.
            Defaults to "gray_diagonal".
        :type noise_option: Literal["none", "gray_diagonal", "gray_band", "erase_diagonal", "erase_band"]
        :param square_scaling: If `True`, the min and max values of the birth (x) and death (y) axis of the plot
            are equalized. Defaults to `True`. In any case, the scale ration of both axis will be equal.
        :type square_scaling: bool, optional
        :raises ValueError: If the value provided for `noise_option` is not valid.
        """
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
    """Specialized overlay for :class:`Vineyard`. Computes the vineyard from a sequence of point clouds or distance
    matrices by computing the Rips Complex for each of them. Can also provide non-trivial representative 1-cycles for
    each step of the vineyard if enabled.

    A trivial representative cycle is a cycle representing a persistence bar of length 0. That is, the cycle is born
    and killed at the same time.
    """

    def __init__(self, store_point_coordinates: bool = False, store_cycles: bool = False):
        """Constructs an empty class with given options. To construct the vineyard, call first :meth:`initialize` to
        build the initial persistent diagram and then a sequence of :meth:`update` to add a sequence of layers to
        the vineyard. Other constructor alternatives are :meth:`from_files` and :meth:`from_tensors`.

        :param store_point_coordinates: Optional. If `True`, the given point clouds for :meth:`initialize` and
            :meth:`update` are copied and stored inside the class. Necessary for :meth:`get_points` and
            :meth:`plot_1D_representative_cycles`. Defaults to False.
        :type store_point_coordinates: bool, optional
        :param store_cycles: Optional. If `True`, the non-trivial representative 1-cycles will be computed and stored at
            each step. Necessary for :meth:`get_1D_representative_cycles` and :meth:`plot_1D_representative_cycles`.
            Defaults to False.
        :type store_cycles: bool, optional
        """
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
        number_of_updates: int = None,
        file_type: str = "point_cloud",
        delimiter: str = None,
        store_point_coordinates: bool = False,
        store_cycles: bool = False,
    ):
        """Construct the vineyard from a sequence of files containing either a point cloud or a distance matrix.
        All files have to be in the same directory and have the same name, except for a number representing their order.
        The numbers have to represent a continuous sequence of integers starting with 0 (even if the argument
        `first_index` is not set to 0). Do not prefix the numbers with various numbers of zeros (i.e., you can do
        `some_name_v009`, `some_name_v0010` but **not** `some_name_v009`, `some_name_v010`).
        From one file to the next, the order of points/distances have to be preserved. I.e., it is assumed that the
        point at line 2 in file `n` is the same point at line 2 in file `n + 1` just with different coordinates.

        File format:
        - Point cloud: plain text where each line represents a different point. A point is given by the ordered
        sequence of its coordinates separated by the same delimiter. Has to be readable by `numpy.loadtxt`.
        - Distance matrix: see file format for :meth:`gudhi.reader_utils.read_lower_triangular_matrix_from_csv_file`.

        :param path_prefix: Part of the file path before the file number. E.g., if files are named
            "./path/to/data_0*_v1.txt" with * being 0, 1, ..., `n`, then `path_prefix` should be set to
            "./path/to/simulation/data_0".
        :type path_prefix: str
        :param path_suffix: Part of the file path after the file number. E.g., if files are named
            "./path/to/data_0*_v1.txt" with * being 0, 1, ..., `n`, then `path_suffix` should be set to
            "_v1.txt". Defaults to ".txt".
        :type path_suffix: str
        :param first_index: Optional. Number of the file to start with, i.e. all files with number strictly lower than
            `first_index` are ignored. Defaults to 0.
        :type first_index: int, optional
        :param number_of_updates: Optional. Number of files to be token into account after `first_index`. That is, any
            file with number strictly greater than `(first_index + number_of_updates)` is ignored. Defaults to `None`
            (i.e., all files which can be found above `first_index`).
        :type number_of_updates: int, optional
        :param file_type: Indicates the content of the file. Has to be either `"point_cloud"` or `"distance_matrix"`.
            Defaults to `"point_cloud"`.
        :type file_type: str
        :param delimiter: Optional. If the files contain point clouds and the coordinates are separated with something
            else than a blank space, the delimiter should be given here. Defaults to `None`.
        :type delimiter: str, optional
        :param store_point_coordinates: Optional and only possible if `file_type` is `"point_cloud"`. If `True`, the
            given point clouds are copied and stored inside the class. Necessary for :meth:`get_points` and
            :meth:`plot_1D_representative_cycles`. Defaults to False.
        :type store_point_coordinates: bool, optional
        :param store_cycles: Optional. If `True`, the non-trivial representative 1-cycles will be computed and stored at
            each step. Necessary for :meth:`get_1D_representative_cycles` and :meth:`plot_1D_representative_cycles`.
            Defaults to False.
        :type store_cycles: bool, optional
        :raises FileNotFoundError: If the file `path_prefix + 0 + path_suffix` is not found and if any file with
            number between `first_index` and `first_index + number_of_updates` is not found.
        :raises ValueError: If `file_type` is invalid.
        :return: Class containing the constructed vineyard.
        :rtype: PointCloudRipsVineyard
        """
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
        number_of_updates: int = None,
        data_type: str = "point_cloud",
        store_point_coordinates: bool = False,
        store_cycles: bool = False,
    ):
        """Constructs the vineyard from a sequence of arrays containing either point clouds or distance matrices.

        :param data: A list or array of arrays corresponding to the format accepted by :class:`RipsComplex`. The first
            axis should correspond to the ordered steps of the vineyard.
        :type data: ArrayLike
        :param first_index: Optional. Index of `data` to start with, i.e. all arrays at index strictly lower than
            `first_index` in `data` are ignored. Defaults to 0.
        :type first_index: int, optional
        :param number_of_updates: Optional. Number of indices to be token into account after `first_index`. That is, any
            array at index strictly greater than `(first_index + number_of_updates)` in `data` is ignored. Defaults to
            `None` (i.e., all array which can be found above `first_index`).
        :type number_of_updates: int, optional
        :param data_type: Indicates the content of the arrays. Has to be either `"point_cloud"` or `"distance_matrix"`.
            Defaults to `"point_cloud"`.
        :type data_type: str
        :param store_point_coordinates: Optional and only possible if `data_type` is `"point_cloud"`. If `True`, the
            given point clouds are copied and stored inside the class. Necessary for :meth:`get_points` and
            :meth:`plot_1D_representative_cycles`. Defaults to False.
        :type store_point_coordinates: bool, optional
        :param store_cycles: Optional. If `True`, the non-trivial representative 1-cycles will be computed and stored at
            each step. Necessary for :meth:`get_1D_representative_cycles` and :meth:`plot_1D_representative_cycles`.
            Defaults to False.
        :type store_cycles: bool, optional
        :raises IndexError: If `first_index` is out of range in `data`.
        :raises ValueError: If `data_type` is invalid.
        :return: Class containing the constructed vineyard.
        :rtype: PointCloudRipsVineyard
        """
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
        data_type: str = "point_cloud",
        number_of_updates: int = 0,
    ):
        """Initializes the vineyard with the first persistence barcode. If another vineyard was initialized before,
        it will be completely replaced.

        :param data: Array corresponding to the format accepted by :class:`RipsComplex`. Alternative to `path`, so both
            should not be provided. Defaults to `None`.
        :type data: ArrayLike or `None`
        :param path: Path to a file containing either a point cloud or a distance matrix. The point cloud has to be
            readable by `numpy.loadtxt` and the distance matrix by
            :meth:`gudhi.reader_utils.read_lower_triangular_matrix_from_csv_file`. Alternative to `data`, so both
            should not be provided. Defaults to `None`.
        :type path: str or `None`
        :param delimiter: Optional. If `path` was provided for a point cloud and the coordinates of the points are
            separated by something else than a blank space, the delimiter should be indicated here. Defaults to `None`.
        :type delimiter: str, optional
        :param data_type: Indicates the content of the array/file. Has to be either `"point_cloud"` or
            `"distance_matrix"`. Defaults to `"point_cloud"`.
        :type data_type: str
        :param number_of_updates: Optional (for optimization purposes). Will allocate memory space for
            `number_of_updates` additional steps in the vineyard after initialization. Defaults to 0.
        :type number_of_updates: int, optional
        :raises ValueError: If both `data` and `path` are provided or none of both.

        .. note::
            If `data_type` is `"distance_matrix"` and `store_point_coordinates` was set to `True` at construction, an
            empty list will be initialized but nothing will be stored in it. So :meth:`get_points` will return an empty
            list and :meth:`plot_1D_representative_cycles` will raise an `IndexError`. A warning will be raised for
            this purpose.
        """
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
        """Adds a layer to the current vineyard by updating the persistence diagram such that it corresponds to the
        Rips filtration represented by the given point cloud or distance matrix.

        :param data: Array corresponding to the format accepted by :class:`RipsComplex`. Alternative to `path`, so both
            should not be provided. Defaults to `None`.
        :type data: ArrayLike or `None`
        :param path: Path to a file containing either a point cloud or a distance matrix. The point cloud has to be
            readable by `numpy.loadtxt` and the distance matrix by
            :meth:`gudhi.reader_utils.read_lower_triangular_matrix_from_csv_file`. Alternative to `data`, so both
            should not be provided. Defaults to `None`.
        :type path: str or `None`
        :param delimiter: Optional. If `path` was provided for a point cloud and the coordinates of the points are
            separated by something else than a blank space, the delimiter should be indicated here. Defaults to `None`.
        :type delimiter: str, optional
        :param data_type: Indicates the content of the array/file. Has to be either `"point_cloud"` or
            `"distance_matrix"`. Defaults to `"point_cloud"`.
        :type data_type: str
        :raises ValueError: If both `data` and `path` are provided or none of both.

        .. note::
            If `data_type` is `"distance_matrix"` and `store_point_coordinates` was set to `True` at construction, an
            empty list will be initialized but nothing will be stored in it. So :meth:`get_points` will return an empty
            list and :meth:`plot_1D_representative_cycles` will raise an `IndexError`. A warning will be raised for
            this purpose.
        """
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
        """Returns the list of read-only and unfiltered vine views. See :meth:`get_current_vineyard` for a more flexible
        output. The format of the list is `dimension x vine number x update number x (birth, death)`, e.g.,
        `vineyard[a][b][c][0]` returns the birth value of the `a`-dimensional vine number `b` at step `c`
        (initialization is step 0), while `vineyard[a][b][c][1]` returns the corresponding death value.

        :param dim: Optional. If provided, the sub-array at index `dim` (corresponding to the vines of dimension `dim`)
            is returned instead of the whole vineyard. Defaults to `None`.
        :type dim: int, optional
        :return: List of read-only views.
        :rtype: list[read-only np.ndarray] (default) or read-only np.ndarray (if `dim` provided)

        .. note::
            As the returned vines are only views, they will get destroyed/invalidated if this class gets destroyed.
        """
        return self._vineyard.get_current_vineyard_view(dim)

    def get_current_vineyard(
        self, dim: int = None, min_bar_length: np.number = -1
    ) -> list[np.ndarray] | np.ndarray:
        """Returns a copy of the current vineyard. If no copy is desired, see :meth:`get_current_vineyard_view`.
        The format of the returned list is `dimension x vine number x update number x (birth, death)`, e.g.,
        `vineyard[a][b][c][0]` returns the birth value of the `a`-dimensional vine number `b` at step `c`
        (initialization is step 0), while `vineyard[a][b][c][1]` returns the corresponding death value.

        :param dim: Optional. If provided, only the vines at given dimension are copied and the output format looses
            the first axes. Defaults to `None`.
        :type dim: int, optional
        :param min_bar_length: Optional. If provided, only vines with at least one coordinate corresponding to a bar of
            length equal or higher than `min_bar_length` are copied. Defaults to -1 (i.e. all vines are copied).
        :type min_bar_length: Any numerical type coercible to the filtration value type, optional
        :return: A copy of the current vineyard.
        :rtype: list[np.ndarray] (default) or np.ndarray (if `dim` provided)
        """
        return self._vineyard.get_current_vineyard(dim, min_bar_length)

    def _get_vertices(self, cpx: list[np.ndarray], i: int) -> set[int]:
        if cpx[i].shape[0] == 0:
            return {i}
        # not optimal for high dimensions, but enables no restrictions on cpx order
        return {v for idx in cpx[i] for v in self._get_vertices(cpx, idx)}

    def get_complex(self, as_vertices: bool = False) -> tuple[list[np.ndarray], np.ndarray]:
        """Returns the Rips complex used to construct the vineyard in the form of `(boundaries, dimensions)`.
        The simplices in the boundaries are represented by their own index in the boundaries and dimensions.
        E.g., `[[], [], [], [0, 1], [0, 2], [1, 2], [3, 4, 5]]`, represents the complex containing a triangle
        and all its faces and `[0, 0, 0, 1, 1, 1, 2]` would be the corresponding dimension array.

        :param as_vertices: Optional. If `True`, instead of their indices, the simplices in the boundaries are
            represented by their vertices. The label of the vertices correspond to the position of the point in the
            original point clouds or distance matrices. Defaults to False.
        :type as_vertices: bool, optional
        :return: Pair of boundaries and dimension array with matching indices.
        :rtype: tuple[list[np.ndarray], np.ndarray]

        .. note::
            Any vertex will have as index the position it had when the point cloud or distance matrix was provided.

        .. note::
            As the Rips complex is simplicial, every value at `dimensions[i]` will be equal to
            `max(len(boundaries[i]) - 1, 0)`. The array is just provided for consistency with other more general
            methods potentially asking for both.
        """
        cpx, dims = t._build_boundary_matrix_from_complex(self._complex)

        if not as_vertices:
            return cpx, dims
        cpx = [
            np.asarray(sorted([v for v in self._get_vertices(cpx, i)]))
            for i, _ in enumerate(cpx)
        ]
        return cpx, dims

    def get_points(self, step: int = None) -> np.ndarray:
        """If `store_point_coordinates` was set to `True` at construction and the provided data were point clouds,
        returns the points stored at each step. The order of the original point clouds is preserved. The array format
        is `step x point number x number of coordinates`.

        :param step: Optional. If provided, only the point cloud at given step is returned (first axis of the format
            is lost). Defaults to `None`.
        :type step: int, optional
        :raises NotImplementedError: If `store_point_coordinates` was set to `False` at construction.
        :return: Stored point cloud(s).
        :rtype: np.ndarray
        """
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

    # TODO: replacing (cycle, length) array with map {bar index, cycle} to match cycles and bars
    def get_1D_representative_cycles(
        self, step: int = None, min_bar_length: np.number = 0
    ) -> list[list[tuple[np.ndarray, np.number]]] | list[tuple[np.ndarray, np.number]]:
        """If `store_cycles` was set to `True` at construction, returns the stored non-trivial representative 1-cycles.
        The format is `step x cycle number x (cycle, bar length)`.

        :param step: Optional. If provided, only the cycles at given step are returned (first axis of the format
            is lost). Defaults to `None`.
        :type step: int, optional
        :param min_bar_length: Optional. If provided, only cycles corresponding to bars with length at least
            `min_bar_length` are returned. Defaults to 0.
        :type min_bar_length: np.number, optional
        :raises NotImplementedError: If `store_cycles` was set to `False` at construction.
        :return: Stored 1-cycles with bar length.
        :rtype: list[list[tuple[np.ndarray, np.number]]] (default) or list[tuple[np.ndarray, np.number]]
            (if `step` was provided)
        """
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
        """Plots the current vineyard, except for completely trivial vines (vines where all coordinates are on
        the diagonal).

        :param dim: Optional. If provided, only plots the vines of the given dimension. Defaults to `None`.
        :type dim: int, optional
        :param max_dim: Optional. If provided, only plots the vines with the given maximal dimension included. If `dim`
            was already provided, `max_dim` is ignored. Defaults to `None`.
        :type max_dim: int, optional
        :param min_bar_length: Optional. If provided, only vines with at least one coordinate corresponding to a bar of
            length equal or higher than `min_bar_length` are plotted. Defaults to -1 (i.e., all vines except trivial
            ones are plotted). Note that trivial vines are never plotted even if `min_bar_length` is set to 0.
        :type min_bar_length: Any numerical type coercible to the filtration value type, optional
        :param noise_option: Value describing how "noisy" parts of a vine should be treated. There are 5 options for
            now:
            - "none": Nothing particular is done.
            - "gray_diagonal": Every part of a vine contained in the diagonal will be plotted with a gray overline.
            - "gray_band": Every part of a vine where all coordinates differ by less than `min_bar_length` are plotted
                with a gray overline.
            - "erase_diagonal": Every part of a vine contained in the diagonal are not plotted.
            - "erase_band": Every part of a vine where all coordinates differ by less than `min_bar_length` are not
                plotted.
            Defaults to "gray_diagonal".
        :type noise_option: Literal["none", "gray_diagonal", "gray_band", "erase_diagonal", "erase_band"]
        :param square_scaling: If `True`, the min and max values of the birth (x) and death (y) axis of the plot
            are equalized. Defaults to `True`. In any case, the scale ration of both axis will be equal.
        :type square_scaling: bool, optional
        :raises ValueError: If the value provided for `noise_option` is not valid.
        """
        self._vineyard.plot_vineyards(
            dim, max_dim, min_bar_length, noise_option, square_scaling
        )

    def _plot_cycle(self, axes, cycle, points, cpx, c):
        for u, v in [(points[cpx[idx][0]], points[cpx[idx][1]]) for idx in cycle]:
            if points.shape[1] == 2:
                axes.plot([u[0], v[0]], [u[1], v[1]], color=c)
            else:
                axes.plot([u[0], v[0]], [u[1], v[1]], [u[2], v[2]], color=c)

    def plot_1D_representative_cycles(
        self, step: int, index: int = None, min_bar_length: np.number = -1
    ):
        """If `store_cycles` and `store_point_coordinates` were set to `True` at construction, plots the representative
        1-cycles at given step, except for those representing bars of length 0. The point coordinates must be 2 or
        3-dimensional.

        :param step: Vineyard step to plot.
        :type step: int
        :param index: Optional. If provided, only plots the cycle at given index at given step. Defaults to `None`.
        :type index: int, optional
        :param min_bar_length: Optional. If provided, only plots 1-cycles corresponding to bars with length at least
            `min_bar_length`. Defaults to -1 (i.e., all non-trivial ones).
        :type min_bar_length: Any numerical type coercible to the filtration value type, optional
        :raises NotImplementedError: If `store_cycles` or `store_point_coordinates` were set to `False` at construction.
        :raises ValueError: If the number of coordinates of a point is different from 2 or 3.
        :raises IndexError: If `index` is provided and is out of range.
        """
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

        cmap_f = mpl.colormaps["tab10"]
        cmap = cmap_f(np.arange(10))
        cmap_f = mpl.colormaps["Set2"]
        cmap = np.concatenate((cmap, cmap_f(np.arange(8))))
        cmap_f = mpl.colormaps["Dark2"]
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
