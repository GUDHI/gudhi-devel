# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau, Bertrand Michel
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - 2020/02 Theo Lacombe: Added more options for improved rendering and more flexibility.
#   - YYYY/MM Author: Description of the modification

from os import path
from math import isfinite
import numpy as np
from functools import lru_cache

from gudhi.reader_utils import read_persistence_intervals_in_dimension
from gudhi.reader_utils import read_persistence_intervals_grouped_by_dimension

__author__ = "Vincent Rouvreau, Bertrand Michel, Theo Lacombe"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

_gudhi_matplotlib_use_tex = True

def __min_birth_max_death(persistence, band=0.0):
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
    return (min_birth, max_death)


def _array_handler(a):
    '''
    :param a: if array, assumes it is a (n x 2) np.array and return a
                persistence-compatible list (padding with 0), so that the
                plot can be performed seamlessly.
    '''
    if isinstance(a[0][1], np.float64) or isinstance(a[0][1], float):
        return [[0, x] for x in a]
    else:
        return a

@lru_cache(maxsize=1)
def _matplotlib_can_use_tex():
    """This function returns True if matplotlib can deal with LaTeX, False otherwise.
    The returned value is cached.
    """
    try:
        from matplotlib import checkdep_usetex
        return checkdep_usetex(True)
    except ImportError:
        print("This function is not available, you may be missing matplotlib.")


def plot_persistence_barcode(
    persistence=[],
    persistence_file="",
    alpha=0.6,
    max_intervals=1000,
    max_barcodes=1000,
    inf_delta=0.1,
    legend=False,
    colormap=None,
    axes=None,
    fontsize=16,
):
    """This function plots the persistence bar code from persistence values list
    , a np.array of shape (N x 2) (representing a diagram 
    in a single homology dimension), 
    or from a `persistence diagram <fileformats.html#persistence-diagram>`_ file.

    :param persistence: Persistence intervals values list. Can be grouped by dimension or not.
    :type persistence: an array of (dimension, array of (birth, death)) or an array of (birth, death).
    :param persistence_file: A `persistence diagram <fileformats.html#persistence-diagram>`_ file style name
        (reset persistence if both are set).
    :type persistence_file: string
    :param alpha: barcode transparency value (0.0 transparent through 1.0
        opaque - default is 0.6).
    :type alpha: float.
    :param max_intervals: maximal number of intervals to display.
        Selected intervals are those with the longest life time. Set it
        to 0 to see all. Default value is 1000.
    :type max_intervals: int.
    :param inf_delta: Infinity is placed at :code:`((max_death - min_birth) x
        inf_delta)` above :code:`max_death` value. A reasonable value is
        between 0.05 and 0.5 - default is 0.1.
    :type inf_delta: float.
    :param legend: Display the dimension color legend (default is False).
    :type legend: boolean.
    :param colormap: A matplotlib-like qualitative colormaps. Default is None
        which means :code:`matplotlib.cm.Set1.colors`.
    :type colormap: tuple of colors (3-tuple of float between 0. and 1.).
    :param axes: A matplotlib-like subplot axes. If None, the plot is drawn on
        a new set of axes.
    :type axes: `matplotlib.axes.Axes`
    :param fontsize: Fontsize to use in axis.
    :type fontsize: int
    :returns: (`matplotlib.axes.Axes`): The axes on which the plot was drawn.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib import rc
        if _gudhi_matplotlib_use_tex and _matplotlib_can_use_tex():
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')
        else:
            plt.rc('text', usetex=False)
            plt.rc('font', family='DejaVu Sans')

        if persistence_file != "":
            if path.isfile(persistence_file):
                # Reset persistence
                persistence = []
                diag = read_persistence_intervals_grouped_by_dimension(
                    persistence_file=persistence_file
                )
                for key in diag.keys():
                    for persistence_interval in diag[key]:
                        persistence.append((key, persistence_interval))
            else:
                print("file " + persistence_file + " not found.")
                return None

        persistence = _array_handler(persistence)

        if max_barcodes != 1000:
            print("Deprecated parameter. It has been replaced by max_intervals")
            max_intervals = max_barcodes

        if max_intervals > 0 and max_intervals < len(persistence):
            # Sort by life time, then takes only the max_intervals elements
            persistence = sorted(
                persistence,
                key=lambda life_time: life_time[1][1] - life_time[1][0],
                reverse=True,
            )[:max_intervals]
        
        if colormap == None:
            colormap = plt.cm.Set1.colors
        if axes == None:
            fig, axes = plt.subplots(1, 1)

        persistence = sorted(persistence, key=lambda birth: birth[1][0])

        (min_birth, max_death) = __min_birth_max_death(persistence)
        ind = 0
        delta = (max_death - min_birth) * inf_delta
        # Replace infinity values with max_death + delta for bar code to be more
        # readable
        infinity = max_death + delta
        axis_start = min_birth - delta
        # Draw horizontal bars in loop
        for interval in reversed(persistence):
            if float(interval[1][1]) != float("inf"):
                # Finite death case
                axes.barh(
                    ind,
                    (interval[1][1] - interval[1][0]),
                    height=0.8,
                    left=interval[1][0],
                    alpha=alpha,
                    color=colormap[interval[0]],
                    linewidth=0,
                )
            else:
                # Infinite death case for diagram to be nicer
                axes.barh(
                    ind,
                    (infinity - interval[1][0]),
                    height=0.8,
                    left=interval[1][0],
                    alpha=alpha,
                    color=colormap[interval[0]],
                    linewidth=0,
                )
            ind = ind + 1

        if legend:
            dimensions = list(set(item[0] for item in persistence))
            axes.legend(
                handles=[
                    mpatches.Patch(color=colormap[dim], label=str(dim))
                    for dim in dimensions
                ],
                loc="lower right",
            )

        axes.set_title("Persistence barcode", fontsize=fontsize)

        # Ends plot on infinity value and starts a little bit before min_birth
        axes.axis([axis_start, infinity, 0, ind])
        return axes

    except ImportError:
        print("This function is not available, you may be missing matplotlib.")


def plot_persistence_diagram(
    persistence=[],
    persistence_file="",
    alpha=0.6,
    band=0.0,
    max_intervals=1000,
    max_plots=1000,
    inf_delta=0.1,
    legend=False,
    colormap=None,
    axes=None,
    fontsize=16,
    greyblock=True
):
    """This function plots the persistence diagram from persistence values
    list, a np.array of shape (N x 2) representing a diagram in a single
    homology dimension, or from a `persistence diagram <fileformats.html#persistence-diagram>`_ file`.

    :param persistence: Persistence intervals values list. Can be grouped by dimension or not.
    :type persistence: an array of (dimension, array of (birth, death)) or an array of (birth, death).
    :param persistence_file: A `persistence diagram <fileformats.html#persistence-diagram>`_ file style name
        (reset persistence if both are set).
    :type persistence_file: string
    :param alpha: plot transparency value (0.0 transparent through 1.0
        opaque - default is 0.6).
    :type alpha: float.
    :param band: band (not displayed if :math:`\leq` 0. - default is 0.)
    :type band: float.
    :param max_intervals: maximal number of intervals to display.
        Selected intervals are those with the longest life time. Set it
        to 0 to see all. Default value is 1000.
    :type max_intervals: int.
    :param inf_delta: Infinity is placed at :code:`((max_death - min_birth) x
        inf_delta)` above :code:`max_death` value. A reasonable value is
        between 0.05 and 0.5 - default is 0.1.
    :type inf_delta: float.
    :param legend: Display the dimension color legend (default is False).
    :type legend: boolean.
    :param colormap: A matplotlib-like qualitative colormaps. Default is None
        which means :code:`matplotlib.cm.Set1.colors`.
    :type colormap: tuple of colors (3-tuple of float between 0. and 1.).
    :param axes: A matplotlib-like subplot axes. If None, the plot is drawn on
        a new set of axes.
    :type axes: `matplotlib.axes.Axes`
    :param fontsize: Fontsize to use in axis.
    :type fontsize: int
    :param greyblock: if we want to plot a grey patch on the lower half plane for nicer rendering. Default True.
    :type greyblock: boolean
    :returns: (`matplotlib.axes.Axes`): The axes on which the plot was drawn.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib import rc
        if _gudhi_matplotlib_use_tex and _matplotlib_can_use_tex():
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')
        else:
            plt.rc('text', usetex=False)
            plt.rc('font', family='DejaVu Sans')

        if persistence_file != "":
            if path.isfile(persistence_file):
                # Reset persistence
                persistence = []
                diag = read_persistence_intervals_grouped_by_dimension(
                    persistence_file=persistence_file
                )
                for key in diag.keys():
                    for persistence_interval in diag[key]:
                        persistence.append((key, persistence_interval))
            else:
                print("file " + persistence_file + " not found.")
                return None

        persistence = _array_handler(persistence)

        if max_plots != 1000:
            print("Deprecated parameter. It has been replaced by max_intervals")
            max_intervals = max_plots

        if max_intervals > 0 and max_intervals < len(persistence):
            # Sort by life time, then takes only the max_intervals elements
            persistence = sorted(
                persistence,
                key=lambda life_time: life_time[1][1] - life_time[1][0],
                reverse=True,
            )[:max_intervals]

        if colormap == None:
            colormap = plt.cm.Set1.colors
        if axes == None:
            fig, axes = plt.subplots(1, 1)

        (min_birth, max_death) = __min_birth_max_death(persistence, band)
        delta = (max_death - min_birth) * inf_delta
        # Replace infinity values with max_death + delta for diagram to be more
        # readable
        infinity = max_death + delta
        axis_end = max_death + delta / 2
        axis_start = min_birth - delta

        # bootstrap band
        if band > 0.0:
            x = np.linspace(axis_start, infinity, 1000)
            axes.fill_between(x, x, x + band, alpha=alpha, facecolor="red")
        # lower diag patch
        if greyblock:
            axes.add_patch(mpatches.Polygon([[axis_start, axis_start], [axis_end, axis_start], [axis_end, axis_end]], fill=True, color='lightgrey'))
        # Draw points in loop
        pts_at_infty = False  # Records presence of pts at infty
        for interval in reversed(persistence):
            if float(interval[1][1]) != float("inf"):
                # Finite death case
                axes.scatter(
                    interval[1][0],
                    interval[1][1],
                    alpha=alpha,
                    color=colormap[interval[0]],
                )
            else:
                pts_at_infty = True
                # Infinite death case for diagram to be nicer
                axes.scatter(
                    interval[1][0], infinity, alpha=alpha, color=colormap[interval[0]]
                )
        if pts_at_infty:
            # infinity line and text
            axes.plot([axis_start, axis_end], [axis_start, axis_end], linewidth=1.0, color="k")
            axes.plot([axis_start, axis_end], [infinity, infinity], linewidth=1.0, color="k", alpha=alpha)
            # Infinity label
            yt = axes.get_yticks()
            yt = yt[np.where(yt < axis_end)] # to avoid ploting ticklabel higher than infinity
            yt = np.append(yt, infinity)
            ytl = ["%.3f" % e for e in yt]  # to avoid float precision error
            ytl[-1] = r'$+\infty$'
            axes.set_yticks(yt)
            axes.set_yticklabels(ytl)

        if legend:
            dimensions = list(set(item[0] for item in persistence))
            axes.legend(
                handles=[
                    mpatches.Patch(color=colormap[dim], label=str(dim))
                    for dim in dimensions
                ]
            )

        axes.set_xlabel("Birth", fontsize=fontsize)
        axes.set_ylabel("Death", fontsize=fontsize)
        axes.set_title("Persistence diagram", fontsize=fontsize)
        # Ends plot on infinity value and starts a little bit before min_birth
        axes.axis([axis_start, axis_end, axis_start, infinity + delta/2])
        return axes

    except ImportError:
        print("This function is not available, you may be missing matplotlib.")


def plot_persistence_density(
    persistence=[],
    persistence_file="",
    nbins=300,
    bw_method=None,
    max_intervals=1000,
    dimension=None,
    cmap=None,
    legend=False,
    axes=None,
    fontsize=16,
    greyblock=False
):
    """This function plots the persistence density from persistence
    values list, np.array of shape (N x 2) representing a diagram
    in a single homology dimension,
    or from a `persistence diagram <fileformats.html#persistence-diagram>`_ file.
    Be aware that this function does not distinguish the dimension, it is
    up to you to select the required one. This function also does not handle
    degenerate data set (scipy correlation matrix inversion can fail).

    :Requires: `SciPy <installation.html#scipy>`_

    :param persistence: Persistence intervals values list. 
                        Can be grouped by dimension or not.
    :type persistence: an array of (dimension, array of (birth, death)) 
                        or an array of (birth, death).
    :param persistence_file: A `persistence diagram <fileformats.html#persistence-diagram>`_
        file style name (reset persistence if both are set).
    :type persistence_file: string
    :param nbins: Evaluate a gaussian kde on a regular grid of nbins x
        nbins over data extents (default is 300)
    :type nbins: int.
    :param bw_method: The method used to calculate the estimator
        bandwidth. This can be 'scott', 'silverman', a scalar constant
        or a callable. If a scalar, this will be used directly as
        kde.factor. If a callable, it should take a gaussian_kde
        instance as only parameter and return a scalar. If None
        (default), 'scott' is used. See
        `scipy.stats.gaussian_kde documentation
        <http://scipy.github.io/devdocs/generated/scipy.stats.gaussian_kde.html>`_
        for more details.
    :type bw_method: str, scalar or callable, optional.
    :param max_intervals: maximal number of points used in the density
        estimation.
        Selected intervals are those with the longest life time. Set it
        to 0 to see all. Default value is 1000.
    :type max_intervals: int.
    :param dimension: the dimension to be selected in the intervals
        (default is None to mix all dimensions).
    :type dimension: int.
    :param cmap: A matplotlib colormap (default is
        matplotlib.pyplot.cm.hot_r).
    :type cmap: cf. matplotlib colormap.
    :param legend: Display the color bar values (default is False).
    :type legend: boolean.
    :param axes: A matplotlib-like subplot axes. If None, the plot is drawn on
        a new set of axes.
    :type axes: `matplotlib.axes.Axes`
    :param fontsize: Fontsize to use in axis.
    :type fontsize: int
    :param greyblock: if we want to plot a grey patch on the lower half plane 
                         for nicer rendering. Default False.
    :type greyblock: boolean
    :returns: (`matplotlib.axes.Axes`): The axes on which the plot was drawn.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from scipy.stats import kde
        from matplotlib import rc
        if _gudhi_matplotlib_use_tex and _matplotlib_can_use_tex():
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')
        else:
            plt.rc('text', usetex=False)
            plt.rc('font', family='DejaVu Sans')

        if persistence_file != "":
            if dimension is None:
                # All dimension case
                dimension = -1
            if path.isfile(persistence_file):
                persistence_dim = read_persistence_intervals_in_dimension(
                    persistence_file=persistence_file, only_this_dim=dimension
                )
            else:
                print("file " + persistence_file + " not found.")
                return None

        if len(persistence) > 0:
            persistence = _array_handler(persistence)
            persistence_dim = np.array(
                [
                    (dim_interval[1][0], dim_interval[1][1])
                    for dim_interval in persistence
                    if (dim_interval[0] == dimension) or (dimension is None)
                ]
            )

        persistence_dim = persistence_dim[np.isfinite(persistence_dim[:, 1])]
        if max_intervals > 0 and max_intervals < len(persistence_dim):
            # Sort by life time, then takes only the max_intervals elements
            persistence_dim = np.array(
                sorted(
                    persistence_dim,
                    key=lambda life_time: life_time[1] - life_time[0],
                    reverse=True,
                )[:max_intervals]
            )

        # Set as numpy array birth and death (remove undefined values - inf and NaN)
        birth = persistence_dim[:, 0]
        death = persistence_dim[:, 1]

        # default cmap value cannot be done at argument definition level as matplotlib is not yet defined.
        if cmap is None:
            cmap = plt.cm.hot_r
        if axes == None:
            fig, axes = plt.subplots(1, 1)

        # line display of equation : birth = death
        x = np.linspace(death.min(), birth.max(), 1000)
        axes.plot(x, x, color="k", linewidth=1.0)

        # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
        k = kde.gaussian_kde([birth, death], bw_method=bw_method)
        xi, yi = np.mgrid[
            birth.min() : birth.max() : nbins * 1j,
            death.min() : death.max() : nbins * 1j,
        ]
        zi = k(np.vstack([xi.flatten(), yi.flatten()]))

        # Make the plot
        img = axes.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=cmap)

        if greyblock:
            axes.add_patch(mpatches.Polygon([[birth.min(), birth.min()], [death.max(), birth.min()], [death.max(), death.max()]], fill=True, color='lightgrey'))

        if legend:
            plt.colorbar(img, ax=axes)

        axes.set_xlabel("Birth", fontsize=fontsize)
        axes.set_ylabel("Death", fontsize=fontsize)
        axes.set_title("Persistence density", fontsize=fontsize)

        return axes

    except ImportError:
        print(
            "This function is not available, you may be missing matplotlib and/or scipy."
        )
