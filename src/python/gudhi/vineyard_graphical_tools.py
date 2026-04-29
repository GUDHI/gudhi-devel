# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Hannah Schreiber
#
# Copyright (C) 2025 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification


__license__ = "MIT"


import numpy as np
from typing import Literal, Optional
import warnings
import matplotlib as mpl
import matplotlib.pyplot as plt

from gudhi.vineyard import Vineyard, PointCloudRipsVineyard


def _gray_on_band(ax, x, y, z, band):
    gray = (0.75, 0.75, 0.75)
    mask = np.concatenate(([0], np.asarray(y - x <= band), [0]))
    diff = np.diff(mask)
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]
    for s, e in zip(starts, ends):
        ax.plot3D(x[s:e], y[s:e], z[s:e], c=gray)


def _erase_on_band(ax, x, y, z, band, c):
    mask = np.concatenate(([0], np.asarray(y - x > band), [0]))
    diff = np.diff(mask)
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]
    for s, e in zip(starts, ends):
        ax.plot3D(x[s:e], y[s:e], z[s:e], c=c)


def _get_colors_for_all_dim(dim: int):
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
                _gray_on_band(ax, x, y, z, 0)
            case "gray_band":
                ax.plot3D(x, y, z, c=cmap[i])
                _gray_on_band(ax, x, y, z, min_bar_length)
            case "erase_diagonal":
                _erase_on_band(ax, x, y, z, 0, cmap[i])
            case "erase_band":
                _erase_on_band(ax, x, y, z, min_bar_length, cmap[i])
            case _:
                raise ValueError(
                    "argument `noise_option` does not contain a valid literal: "
                    + noise_option
                    + ". Possibilities are: `none`, `gray_diagonal`, `gray_band`, `erase_diagonal` and `erase_band`"
                )


# TODO: option to handle points at infinity
def plot_vineyards(
    vineyard: Vineyard | PointCloudRipsVineyard,
    ax = None,
    dim: Optional[int] = None,
    max_dim: Optional[int] = None,
    min_bar_length: np.number = -1,
    noise_option: Literal[
        "none", "gray_diagonal", "gray_band", "erase_diagonal", "erase_band"
    ] = "gray_diagonal",
    square_scaling: Optional[bool] = True,
):
    """Plots the given vineyard, except for completely trivial vines (vines where all coordinates are on
    the diagonal). The points at infinity are mapped to a finite point a bit away from the other points.

    :param vineyard: Vineyard class to plot.
    :type vineyard: :class:`~gudhi.vineyard.Vineyard` or :class:`~gudhi.vineyard.PointCloudRipsVineyard`
    :param ax: Optional. Matplotlib 3D-axis to plot into. If not provided, a new one is created. Defaults to `None`.
    :type ax: `mpl_toolkits.mplot3d.axes3d.Axes3D \
        <https://matplotlib.org/stable/api/toolkits/mplot3d/axes3d.html#mpl_toolkits.mplot3d.axes3d.Axes3D>`_, optional
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
            - "gray_band": Every part of a vine where all coordinates differ by less than `min_bar_length` are \
                plotted with a gray overline.
            - "erase_diagonal": Every part of a vine contained in the diagonal are not plotted.
            - "erase_band": Every part of a vine where all coordinates differ by less than `min_bar_length` are \
                not plotted.

        Defaults to "gray_diagonal".
    :type noise_option: Literal["none", "gray_diagonal", "gray_band", "erase_diagonal", "erase_band"]
    :param square_scaling: If `True`, the min and max values of the birth (x) and death (y) axis of the plot
        are equalized. Defaults to `True`. In any case, the scale ration of both axis will be equal.
    :type square_scaling: bool, optional
    :raises ValueError: If the value provided for `noise_option` is not valid.
    :return: Matplotlib axis into which was plotted.
    :rtype: `mpl_toolkits.mplot3d.axes3d.Axes3D \
        <https://matplotlib.org/stable/api/toolkits/mplot3d/axes3d.html#mpl_toolkits.mplot3d.axes3d.Axes3D>`_
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")

    if dim is None:
        vy = vineyard.get_current_vineyard_view()
        if max_dim is not None:
            vy = vy[: max_dim + 1]
        max_death = max(
            [
                np.max(vines[:, :, 1], where=~np.isinf(vines[:, :, 1]), initial=-1)
                for vines in vy
            ]
        )
        # arbitrary, just to distance them a bit from the finite bars
        inf_v = max_death + 10

        for d, vines in enumerate(vy):
            cmap = _get_colors_for_all_dim(d)
            _plot_vines(ax, vines, inf_v, min_bar_length, noise_option, cmap)
    else:
        vines = vineyard.get_current_vineyard_view(dim=dim)
        max_death = np.max(vines[:, :, 1], where=~np.isinf(vines[:, :, 1]), initial=-1)
        # arbitrary, just to distance them a bit from the finite bars
        inf_v = max_death + 10
        _plot_vines(ax, vines, inf_v, min_bar_length, noise_option, mpl.colormaps["viridis"])

    if square_scaling:
        scaling = np.array([getattr(ax, "get_{}lim".format(dim))() for dim in "xy"])
        ax.auto_scale_xyz(
            *[[np.min(scaling), np.max(scaling)]] * 2,
            [getattr(ax, "get_{}lim".format("z"))()],
        )
    ax.set_aspect("equalxy", "box")

    return ax


def _plot_cycle(axes, cycle, points, cpx, c):
    for u, v in [(points[cpx[idx][0]], points[cpx[idx][1]]) for idx in cycle]:
        if points.shape[1] == 2:
            axes.plot([u[0], v[0]], [u[1], v[1]], color=c)
        else:
            axes.plot([u[0], v[0]], [u[1], v[1]], [u[2], v[2]], color=c)


def plot_1D_representative_cycles(
    vineyard: Vineyard | PointCloudRipsVineyard,
    step: int,
    ax = None,
    index: Optional[int] = None,
    min_bar_length: Optional[np.number] = None,
):
    """If `store_cycles` and `store_point_coordinates` were set to `True` at construction of `vineyard`, plots the
    representative 1-cycles at given step, except for those representing bars of length 0. The point coordinates must
    be 2 or 3-dimensional.

    :param vineyard: Vineyard class from where to plot the cycles.
    :type vineyard: :class:`~gudhi.vineyard.Vineyard` or :class:`~gudhi.vineyard.PointCloudRipsVineyard`
    :param step: Vineyard step to plot.
    :type step: int
    :param ax: Optional. Matplotlib 3D-axis to plot into. If not provided, a new one is created. Defaults to `None`.
    :type ax: `mpl_toolkits.mplot3d.axes3d.Axes3D \
        <https://matplotlib.org/stable/api/toolkits/mplot3d/axes3d.html#mpl_toolkits.mplot3d.axes3d.Axes3D>`_, optional
    :param index: Optional. If provided, only plots the cycle at given index at given step. Defaults to `None`.
    :type index: int, optional
    :param min_bar_length: Optional. If provided, only plots 1-cycles corresponding to bars with length at least
        `min_bar_length`. Defaults to None (i.e., all non-trivial ones).
    :type min_bar_length: Any numerical type coercible to the filtration value type, optional
    :raises NotImplementedError: If `store_cycles` or `store_point_coordinates` were set to `False` at construction.
    :raises ValueError: If the number of coordinates of a point is different from 2 or 3.
    :raises IndexError: If `index` is provided and is out of range.
    :return: Matplotlib axis into which was plotted.
    :rtype: `mpl_toolkits.mplot3d.axes3d.Axes3D \
        <https://matplotlib.org/stable/api/toolkits/mplot3d/axes3d.html#mpl_toolkits.mplot3d.axes3d.Axes3D>`_
    """
    points = vineyard.get_points(step)
    cycles = vineyard.get_1D_representative_cycles(step)
    cpx, _ = vineyard.get_complex()

    if ax is None:
        fig = plt.figure()
        if points.shape[1] == 2:
            ax = fig.add_subplot()
        elif points.shape[1] == 3:
            ax = fig.add_subplot(projection="3d")
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
        if min_bar_length is None:
            for _, c in cycles.items():
                _plot_cycle(ax, c, points, cpx, cmap[i])
                i = (i + 1) % csize
        else:
            vy = vineyard.get_current_vineyard_view(dim=1)
            for (_, idx), c in cycles.items():
                if vy[idx][step][1] - vy[idx][step][0] >= min_bar_length:
                    _plot_cycle(ax, c, points, cpx, cmap[i])
                    i = (i + 1) % csize
    else:
        if min_bar_length is not None:
            warnings.warn(
                "Specified argument `min_bar_length` is ignored.",
                UserWarning,
            )
        _plot_cycle(ax, cycles[index][0], points, cpx, cmap[i])

    gray = (0.25, 0.25, 0.25)
    x = points[:, 0]
    y = points[:, 1]
    if points.shape[1] == 2:
        ax.plot(x, y, linestyle="none", markersize=3, marker="o", c=gray)
    else:
        ax.plot(x, y, points[:, 2], linestyle="none", markersize=3, marker="o", c=gray)

    return ax
