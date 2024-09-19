# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau, Bertrand Michel
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - 2020/02 Theo Lacombe: Added more options for improved rendering and more flexibility.
#   - 2022/11 Vincent Rouvreau: "Automatic" legend display detected by _array_handler that returns if the persistence
#                               was a nx2 array.
#   - YYYY/MM Author: Description of the modification

from os import path
import numpy as np
from functools import lru_cache
import warnings
import errno
import os
import shutil

from gudhi.reader_utils import read_persistence_intervals_in_dimension
from gudhi.reader_utils import read_persistence_intervals_grouped_by_dimension

__author__ = "Vincent Rouvreau, Bertrand Michel, Theo Lacombe"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

_gudhi_matplotlib_use_tex = True


def _min_birth_max_death(persistence, band=0.0):
    """This function returns (min_birth, max_death) from the persistence.

    :param persistence: The persistence to plot.
    :type persistence: list of tuples(dimension, tuple(birth, death)).
    :param band: band
    :type band: float.
    :returns: (float, float) -- (min_birth, max_death).
    """
    # Look for minimum birth date and maximum death date for plot optimisation
    max_death = 0
    min_birth = persistence[0][1][0]
    for interval in reversed(persistence):
        if float(interval[1][1]) != float("inf"):
            if float(interval[1][1]) > max_death:
                max_death = float(interval[1][1])
        if float(interval[1][0]) > max_death:
            max_death = float(interval[1][0])
        if float(interval[1][0]) < min_birth:
            min_birth = float(interval[1][0])
    if band > 0.0:
        max_death += band
    # can happen if only points at inf death
    if min_birth == max_death:
        max_death = max_death + 1.0
    return (min_birth, max_death)


def _array_handler(a):
    """
    :param a: if array, assumes it is a (n x 2) np.array and returns a
                persistence-compatible list (padding with 0), so that the
                plot can be performed seamlessly.
    :returns: * List[dimension, [birth, death]] Persistence, compatible with plot functions, list.
              * boolean Modification status (True if output is different from input)
    """
    if isinstance(a[0][1], (np.floating, float)):
        return [[0, x] for x in a], True
    else:
        return a, False


def _limit_to_max_intervals(persistence, max_intervals, key):
    """This function returns truncated persistence if length is bigger than max_intervals.
    :param persistence: Persistence intervals values list. Can be grouped by dimension or not.
    :type persistence: an array of (dimension, (birth, death)) or an array of (birth, death).
    :param max_intervals: maximal number of intervals to display.
        Selected intervals are those with the longest life time. Set it
        to 0 to see all. Default value is 1000.
    :type max_intervals: int.
    :param key: key function for sort algorithm.
    :type key: function or lambda.
    """
    if max_intervals > 0 and max_intervals < len(persistence):
        warnings.warn(
            "There are %s intervals given as input, whereas max_intervals is set to %s."
            % (len(persistence), max_intervals)
        )
        # Sort by life time, then takes only the max_intervals elements
        return sorted(persistence, key=key, reverse=True)[:max_intervals]
    else:
        return persistence


@lru_cache(maxsize=1)
def _matplotlib_can_use_tex() -> bool:
    """This function returns True if matplotlib can deal with LaTeX, False otherwise.
    The returned value is cached.

    This code is taken
    https://github.com/matplotlib/matplotlib/blob/f975291a008f001047ad8964b15d7d64d2907f1e/lib/matplotlib/__init__.py#L454-L471
    deprecated from matplotlib 3.6 and removed in matplotlib 3.8.0
    """
    from matplotlib import _get_executable_info, ExecutableNotFoundError

    if not shutil.which("tex"):
        warnings.warn("usetex mode requires TeX.")
        return False
    try:
        _get_executable_info("dvipng")
    except ExecutableNotFoundError:
        warnings.warn("usetex mode requires dvipng.")
        return False
    try:
        _get_executable_info("gs")
    except ExecutableNotFoundError:
        warnings.warn("usetex mode requires ghostscript.")
        return False
    return True


def plot_persistence_barcode(
    persistence=[],
    persistence_file="",
    alpha=0.6,
    max_intervals=20000,
    inf_delta=0.1,
    legend=None,
    colormap=None,
    axes=None,
    fontsize=16,
):
    """This function plots the persistence bar code from persistence values list
    , a np.array of shape (N x 2) (representing a diagram
    in a single homology dimension),
    or from a `persistence diagram <fileformats.html#persistence-diagram>`_ file.

    :param persistence: Persistence intervals values list. Can be grouped by dimension or not.
    :type persistence: an array of (dimension, (birth, death)) or an array of (birth, death)
    :param persistence_file: A `persistence diagram <fileformats.html#persistence-diagram>`_ file style name
        (reset persistence if both are set).
    :type persistence_file: string
    :param alpha: barcode transparency value (0.0 transparent through 1.0
        opaque - default is 0.6).
    :type alpha: float
    :param max_intervals: maximal number of intervals to display.
        Selected intervals are those with the longest life time. Set it
        to 0 to see all. Default value is 20000.
    :type max_intervals: int
    :param inf_delta: Infinity is placed at :code:`((max_death - min_birth) x
        inf_delta)` above :code:`max_death` value. A reasonable value is
        between 0.05 and 0.5 - default is 0.1.
    :type inf_delta: float
    :param legend: Display the dimension color legend. Default is None, meaning the legend is displayed if dimension
        is specified in the persistence argument, and not displayed if dimension is not specified.
    :type legend: boolean or None
    :param colormap: A matplotlib-like qualitative colormaps. Default is None
        which means :code:`matplotlib.cm.Set1.colors`.
    :type colormap: tuple of colors (3-tuple of float between 0. and 1.)
    :param axes: A matplotlib-like subplot axes. If None, the plot is drawn on
        a new set of axes.
    :type axes: `matplotlib.axes.Axes`
    :param fontsize: Fontsize to use in axis.
    :type fontsize: int
    :returns: (`matplotlib.axes.Axes`): The axes on which the plot was drawn.
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    if _gudhi_matplotlib_use_tex and _matplotlib_can_use_tex():
        plt.rc("text", usetex=True)
        plt.rc("font", family="serif")
    else:
        plt.rc("text", usetex=False)
        plt.rc("font", family="DejaVu Sans")

    # By default, let's say the persistence is not an array of shape (N x 2) - Can be from a persistence file
    nx2_array = False
    if persistence_file != "":
        if path.isfile(persistence_file):
            # Reset persistence
            persistence = []
            diag = read_persistence_intervals_grouped_by_dimension(persistence_file=persistence_file)
            for key in diag.keys():
                for persistence_interval in diag[key]:
                    persistence.append((key, persistence_interval))
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), persistence_file)

    try:
        persistence, nx2_array = _array_handler(persistence)
        persistence = _limit_to_max_intervals(
            persistence, max_intervals, key=lambda life_time: life_time[1][1] - life_time[1][0]
        )
        (min_birth, max_death) = _min_birth_max_death(persistence)
        persistence = sorted(persistence, key=lambda birth: birth[1][0])
    except IndexError:
        min_birth, max_death = 0.0, 1.0
        pass

    delta = (max_death - min_birth) * inf_delta
    # Replace infinity values with max_death + delta for bar code to be more
    # readable
    infinity = max_death + delta
    axis_start = min_birth - delta

    if axes is None:
        _, axes = plt.subplots(1, 1)
    if colormap is None:
        colormap = plt.cm.Set1.colors

    x = [birth for (dim, (birth, death)) in persistence]
    y = [(death - birth) if death != float("inf") else (infinity - birth) for (dim, (birth, death)) in persistence]
    c = [colormap[dim] for (dim, (birth, death)) in persistence]

    axes.barh(range(len(x)), y, left=x, alpha=alpha, color=c, linewidth=0)

    if legend is None and not nx2_array:
        # By default, if persistence is an array of (dimension, (birth, death)), display the legend
        legend = True

    if legend:
        dimensions = {item[0] for item in persistence}
        axes.legend(
            handles=[mpatches.Patch(color=colormap[dim], label=str(dim)) for dim in dimensions],
            loc="best",
        )

    axes.set_title("Persistence barcode", fontsize=fontsize)
    axes.set_yticks([])
    axes.invert_yaxis()

    # Ends plot on infinity value and starts a little bit before min_birth
    if len(x) != 0:
        axes.set_xlim((axis_start, infinity))
    return axes


def plot_persistence_diagram(
    persistence=[],
    persistence_file="",
    alpha=0.6,
    band=0.0,
    max_intervals=1000000,
    inf_delta=0.1,
    legend=None,
    colormap=None,
    axes=None,
    fontsize=16,
    greyblock=True,
):
    r"""This function plots the persistence diagram from persistence values
    list, a np.array of shape (N x 2) representing a diagram in a single
    homology dimension, or from a `persistence diagram <fileformats.html#persistence-diagram>`_ file`.

    :param persistence: Persistence intervals values list. Can be grouped by dimension or not.
    :type persistence: an array of (dimension, (birth, death)) or an array of (birth, death)
    :param persistence_file: A `persistence diagram <fileformats.html#persistence-diagram>`_ file style name
        (reset persistence if both are set).
    :type persistence_file: string
    :param alpha: plot transparency value (0.0 transparent through 1.0
        opaque - default is 0.6).
    :type alpha: float
    :param band: band (not displayed if :math:`\leq` 0. - default is 0.)
    :type band: float
    :param max_intervals: maximal number of intervals to display.
        Selected intervals are those with the longest life time. Set it
        to 0 to see all. Default value is 1000000.
    :type max_intervals: int
    :param inf_delta: Infinity is placed at :code:`((max_death - min_birth) x
        inf_delta)` above :code:`max_death` value. A reasonable value is
        between 0.05 and 0.5 - default is 0.1.
    :type inf_delta: float
    :param legend: Display the dimension color legend. Default is None, meaning the legend is displayed if dimension
        is specified in the persistence argument, and not displayed if dimension is not specified.
    :type legend: boolean or None
    :param colormap: A matplotlib-like qualitative colormaps. Default is None
        which means :code:`matplotlib.cm.Set1.colors`.
    :type colormap: tuple of colors (3-tuple of float between 0. and 1.)
    :param axes: A matplotlib-like subplot axes. If None, the plot is drawn on
        a new set of axes.
    :type axes: `matplotlib.axes.Axes`
    :param fontsize: Fontsize to use in axis.
    :type fontsize: int
    :param greyblock: if we want to plot a grey patch on the lower half plane for nicer rendering. Default True.
    :type greyblock: boolean
    :returns: (`matplotlib.axes.Axes`): The axes on which the plot was drawn.
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    if _gudhi_matplotlib_use_tex and _matplotlib_can_use_tex():
        plt.rc("text", usetex=True)
        plt.rc("font", family="serif")
    else:
        plt.rc("text", usetex=False)
        plt.rc("font", family="DejaVu Sans")

    # By default, let's say the persistence is not an array of shape (N x 2) - Can be from a persistence file
    nx2_array = False
    if persistence_file != "":
        if path.isfile(persistence_file):
            # Reset persistence
            persistence = []
            diag = read_persistence_intervals_grouped_by_dimension(persistence_file=persistence_file)
            for key in diag.keys():
                for persistence_interval in diag[key]:
                    persistence.append((key, persistence_interval))
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), persistence_file)

    try:
        persistence, nx2_array = _array_handler(persistence)
        persistence = _limit_to_max_intervals(
            persistence, max_intervals, key=lambda life_time: life_time[1][1] - life_time[1][0]
        )
        min_birth, max_death = _min_birth_max_death(persistence, band)
    except IndexError:
        min_birth, max_death = 0.0, 1.0
        pass

    delta = (max_death - min_birth) * inf_delta
    # Replace infinity values with max_death + delta for diagram to be more
    # readable
    infinity = max_death + delta
    axis_end = max_death + delta / 2
    axis_start = min_birth - delta

    if axes is None:
        _, axes = plt.subplots(1, 1)
    if colormap is None:
        colormap = plt.cm.Set1.colors
    # bootstrap band
    if band > 0.0:
        x = np.linspace(axis_start, infinity, 1000)
        axes.fill_between(x, x, x + band, alpha=alpha, facecolor="red")
    # lower diag patch
    if greyblock:
        axes.add_patch(
            mpatches.Polygon(
                [[axis_start, axis_start], [axis_end, axis_start], [axis_end, axis_end]],
                fill=True,
                color="lightgrey",
            )
        )
    # line display of equation : birth = death
    axes.plot([axis_start, axis_end], [axis_start, axis_end], linewidth=1.0, color="k")

    x = [birth for (dim, (birth, death)) in persistence]
    y = [death if death != float("inf") else infinity for (dim, (birth, death)) in persistence]
    c = [colormap[dim] for (dim, (birth, death)) in persistence]

    axes.scatter(x, y, alpha=alpha, color=c)
    if float("inf") in (death for (dim, (birth, death)) in persistence):
        # infinity line and text
        axes.plot([axis_start, axis_end], [infinity, infinity], linewidth=1.0, color="k", alpha=alpha)
        # Infinity label
        yt = axes.get_yticks()
        yt = yt[np.where(yt < axis_end)]  # to avoid plotting ticklabel higher than infinity
        yt = np.append(yt, infinity)
        ytl = ["%.3f" % e for e in yt]  # to avoid float precision error
        ytl[-1] = r"$+\infty$"
        axes.set_yticks(yt)
        axes.set_yticklabels(ytl)

    if legend is None and not nx2_array:
        # By default, if persistence is an array of (dimension, (birth, death)), display the legend
        legend = True

    if legend:
        dimensions = list({item[0] for item in persistence})
        axes.legend(
            handles=[mpatches.Patch(color=colormap[dim], label=str(dim)) for dim in dimensions],
            loc="lower right",
        )

    axes.set_xlabel("Birth", fontsize=fontsize)
    axes.set_ylabel("Death", fontsize=fontsize)
    axes.set_title("Persistence diagram", fontsize=fontsize)
    # Ends plot on infinity value and starts a little bit before min_birth
    axes.axis([axis_start, axis_end, axis_start, infinity + delta / 2])
    return axes


def plot_persistence_density(
    persistence=[],
    persistence_file="",
    nbins=300,
    bw_method=None,
    max_intervals=1000,
    dimension=None,
    cmap=None,
    legend=True,
    axes=None,
    fontsize=16,
    greyblock=False,
):
    """This function plots the persistence density from persistence values list, np.array of shape (N x 2) representing
    a diagram in a single homology dimension, or from a `persistence diagram <fileformats.html#persistence-diagram>`_
    file. Be aware that this function does not distinguish the dimension, it is up to you to select the required one.
    This function also does not handle degenerate data set (scipy correlation matrix inversion can fail).

    :Requires: `SciPy <installation.html#scipy>`_

    :param persistence: Persistence intervals values list. Can be grouped by dimension or not.
    :type persistence: an array of (dimension, (birth, death)) or an array of (birth, death)
    :param persistence_file: A `persistence diagram <fileformats.html#persistence-diagram>`_
        file style name (reset persistence if both are set).
    :type persistence_file: string
    :param nbins: Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents (default is 300)
    :type nbins: int
    :param bw_method: The method used to calculate the estimator bandwidth. This can be 'scott', 'silverman', a scalar
        constant or a callable. If a scalar, this will be used directly as kde.factor. If a callable, it should take a
        gaussian_kde instance as only parameter and return a scalar. If None (default), 'scott' is used. See
        `scipy.stats.gaussian_kde documentation
        <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html>`_
        for more details.
    :type bw_method: str, scalar or callable, optional
    :param max_intervals: maximal number of points used in the density estimation. Selected intervals are those with
        the longest life time. Set it to 0 to see all. Default value is 1000.
    :type max_intervals: int
    :param dimension: the dimension to be selected in the intervals (default is None to mix all dimensions).
    :type dimension: int
    :param cmap: A matplotlib colormap (default is matplotlib.pyplot.cm.hot_r).
    :type cmap: cf. matplotlib colormap
    :param legend: Display the color bar values (default is True).
    :type legend: boolean
    :param axes: A matplotlib-like subplot axes. If None, the plot is drawn on a new set of axes.
    :type axes: `matplotlib.axes.Axes`
    :param fontsize: Fontsize to use in axis.
    :type fontsize: int
    :param greyblock: if we want to plot a grey patch on the lower half plane for nicer rendering. Default False.
    :type greyblock: boolean
    :returns: (`matplotlib.axes.Axes`): The axes on which the plot was drawn.
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from scipy.stats import kde

    if _gudhi_matplotlib_use_tex and _matplotlib_can_use_tex():
        plt.rc("text", usetex=True)
        plt.rc("font", family="serif")
    else:
        plt.rc("text", usetex=False)
        plt.rc("font", family="DejaVu Sans")

    if persistence_file != "":
        if dimension is None:
            # All dimension case
            dimension = -1
        if path.isfile(persistence_file):
            persistence_dim = read_persistence_intervals_in_dimension(
                persistence_file=persistence_file, only_this_dim=dimension
            )
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), persistence_file)

    # default cmap value cannot be done at argument definition level as matplotlib is not yet defined.
    if cmap is None:
        cmap = plt.cm.hot_r
    if axes is None:
        _, axes = plt.subplots(1, 1)

    try:
        # if not read from file but given by an argument
        persistence, _ = _array_handler(persistence)
        persistence_dim = np.array(
            [
                (dim_interval[1][0], dim_interval[1][1])
                for dim_interval in persistence
                if (dim_interval[0] == dimension) or (dimension is None)
            ]
        )
        persistence_dim = persistence_dim[np.isfinite(persistence_dim[:, 1])]
        persistence_dim = np.array(
            _limit_to_max_intervals(
                persistence_dim, max_intervals, key=lambda life_time: life_time[1] - life_time[0]
            )
        )

        # Set as numpy array birth and death (remove undefined values - inf and NaN)
        birth = persistence_dim[:, 0]
        death = persistence_dim[:, 1]
        birth_min = birth.min()
        birth_max = birth.max()
        death_min = death.min()
        death_max = death.max()

        # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
        k = kde.gaussian_kde([birth, death], bw_method=bw_method)
        xi, yi = np.mgrid[
            birth_min : birth_max : nbins * 1j,
            death_min : death_max : nbins * 1j,
        ]
        zi = k(np.vstack([xi.ravel(), yi.ravel()]))
        # Make the plot
        img = axes.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=cmap, shading="auto")
        plot_success = True

    # IndexError on empty diagrams, ValueError on only inf death values
    except (IndexError, ValueError):
        birth_min = 0.0
        birth_max = 1.0
        death_min = 0.0
        death_max = 1.0
        plot_success = False
        pass

    # line display of equation : birth = death
    x = np.linspace(death_min, birth_max, 1000)
    axes.plot(x, x, color="k", linewidth=1.0)

    if greyblock:
        axes.add_patch(
            mpatches.Polygon(
                [[birth_min, birth_min], [death_max, birth_min], [death_max, death_max]],
                fill=True,
                color="lightgrey",
            )
        )

    if plot_success and legend:
        plt.colorbar(img, ax=axes)

    axes.set_xlabel("Birth", fontsize=fontsize)
    axes.set_ylabel("Death", fontsize=fontsize)
    axes.set_title("Persistence density", fontsize=fontsize)

    return axes
