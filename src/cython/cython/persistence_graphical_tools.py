"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Vincent Rouvreau, Bertrand Michel

   Copyright (C) 2016 Inria

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__author__ = "Vincent Rouvreau, Bertrand Michel"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "GPL v3"

def __min_birth_max_death(persistence, band=0.):
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
        if float(interval[1][1]) != float('inf'):
            if float(interval[1][1]) > max_death:
                max_death = float(interval[1][1])
        if float(interval[1][0]) > max_death:
            max_death = float(interval[1][0])
        if float(interval[1][0]) < min_birth:
            min_birth = float(interval[1][0])
    if band > 0.:
        max_death += band
    return (min_birth, max_death)

"""
Only 13 colors for the palette
"""
palette = ['#ff0000', '#00ff00', '#0000ff', '#00ffff', '#ff00ff', '#ffff00',
           '#000000', '#880000', '#008800', '#000088', '#888800', '#880088',
           '#008888']

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
    import os

    def plot_persistence_barcode(persistence=[], persistence_file='', alpha=0.6,
            max_barcodes=1000, inf_delta=0.1, legend=False):
        """This function plots the persistence bar code from persistence values list
        or from a :doc:`persistence file <fileformats>`.

        :param persistence: Persistence values list.
        :type persistence: list of tuples(dimension, tuple(birth, death)).
        :param persistence_file: A :doc:`persistence file <fileformats>` style name
            (reset persistence if both are set).
        :type persistence_file: string
        :param alpha: barcode transparency value (0.0 transparent through 1.0
            opaque - default is 0.6).
        :type alpha: float.
        :param max_barcodes: number of maximal barcodes to be displayed.
            Set it to 0 to see all, Default value is 1000.
            (persistence will be sorted by life time if max_barcodes is set)
        :type max_barcodes: int.
        :param inf_delta: Infinity is placed at ((max_death - min_birth) x
            inf_delta) above the highest point. A reasonable value is between
            0.05 and 0.5 - default is 0.1.
        :type inf_delta: float.
        :param legend: Display the dimension color legend (default is False).
        :type legend: boolean.
        :returns: A matplotlib object containing horizontal bar plot of persistence
            (launch `show()` method on it to display it).
        """
        if persistence_file is not '':
            if os.path.isfile(persistence_file):
                # Reset persistence
                persistence = []
                diag = read_persistence_intervals_grouped_by_dimension(persistence_file=persistence_file)
                for key in diag.keys():
                    for persistence_interval in diag[key]:
                        persistence.append((key, persistence_interval))
            else:
                print("file " + persistence_file + " not found.")
                return None

        if max_barcodes > 0 and max_barcodes < len(persistence):
            # Sort by life time, then takes only the max_plots elements
            persistence = sorted(persistence, key=lambda life_time: life_time[1][1]-life_time[1][0], reverse=True)[:max_barcodes]

        persistence = sorted(persistence, key=lambda birth: birth[1][0])

        (min_birth, max_death) = __min_birth_max_death(persistence)
        ind = 0
        delta = ((max_death - min_birth) * inf_delta)
        # Replace infinity values with max_death + delta for bar code to be more
        # readable
        infinity = max_death + delta
        axis_start = min_birth - delta
        # Draw horizontal bars in loop
        for interval in reversed(persistence):
            if float(interval[1][1]) != float('inf'):
                # Finite death case
                plt.barh(ind, (interval[1][1] - interval[1][0]), height=0.8,
                         left = interval[1][0], alpha=alpha,
                         color = palette[interval[0]],
                         linewidth=0)
            else:
                # Infinite death case for diagram to be nicer
                plt.barh(ind, (infinity - interval[1][0]), height=0.8,
                         left = interval[1][0], alpha=alpha,
                         color = palette[interval[0]],
                         linewidth=0)
            ind = ind + 1

        if legend:
            dimensions = list(set(item[0] for item in persistence))
            plt.legend(handles=[mpatches.Patch(color=palette[dim],
                                               label=str(dim)) for dim in dimensions],
                       loc='lower right')
        plt.title('Persistence barcode')
        # Ends plot on infinity value and starts a little bit before min_birth
        plt.axis([axis_start, infinity, 0, ind])
        return plt

    def plot_persistence_diagram(persistence=[], persistence_file='', alpha=0.6,
            band=0., max_plots=1000, inf_delta=0.1, legend=False):
        """This function plots the persistence diagram from persistence values
        list or from a :doc:`persistence file <fileformats>`.

        :param persistence: Persistence values list.
        :type persistence: list of tuples(dimension, tuple(birth, death)).
        :param persistence_file: A :doc:`persistence file <fileformats>` style name
            (reset persistence if both are set).
        :type persistence_file: string
        :param alpha: plot transparency value (0.0 transparent through 1.0
            opaque - default is 0.6).
        :type alpha: float.
        :param band: band (not displayed if :math:`\leq` 0. - default is 0.)
        :type band: float.
        :param max_plots: maximal number of points to display. Selected points
            are those with the longest life time. Set it to 0 to see all,
            default value is 1000.
        :type max_plots: int.
        :param inf_delta: Infinity is placed at ((max_death - min_birth) x
            inf_delta) above the highest point. A reasonable value is between
            0.05 and 0.5 - default is 0.1.
        :type inf_delta: float.
        :param legend: Display the dimension color legend (default is False).
        :type legend: boolean.
        :returns: A matplotlib object containing diagram plot of persistence
            (launch `show()` method on it to display it).
        """
        if persistence_file is not '':
            if os.path.isfile(persistence_file):
                # Reset persistence
                persistence = []
                diag = read_persistence_intervals_grouped_by_dimension(persistence_file=persistence_file)
                for key in diag.keys():
                    for persistence_interval in diag[key]:
                        persistence.append((key, persistence_interval))
            else:
                print("file " + persistence_file + " not found.")
                return None

        if max_plots > 0 and max_plots < len(persistence):
            # Sort by life time, then takes only the max_plots elements
            persistence = sorted(persistence, key=lambda life_time: life_time[1][1]-life_time[1][0], reverse=True)[:max_plots]

        (min_birth, max_death) = __min_birth_max_death(persistence, band)
        delta = ((max_death - min_birth) * inf_delta)
        # Replace infinity values with max_death + delta for diagram to be more
        # readable
        infinity = max_death + delta
        axis_start = min_birth - delta

        # line display of equation : birth = death
        x = np.linspace(axis_start, infinity, 1000)
        # infinity line and text
        plt.plot(x, x, color='k', linewidth=1.0)
        plt.plot(x, [infinity] * len(x), linewidth=1.0, color='k', alpha=alpha)
        plt.text(axis_start, infinity, r'$\infty$', color='k', alpha=alpha)
        # bootstrap band
        if band > 0.:
            plt.fill_between(x, x, x+band, alpha=alpha, facecolor='red')

        # Draw points in loop
        for interval in reversed(persistence):
            if float(interval[1][1]) != float('inf'):
                # Finite death case
                plt.scatter(interval[1][0], interval[1][1], alpha=alpha,
                            color = palette[interval[0]])
            else:
                # Infinite death case for diagram to be nicer
                plt.scatter(interval[1][0], infinity, alpha=alpha,
                            color = palette[interval[0]])

        if legend:
            dimensions = list(set(item[0] for item in persistence))
            plt.legend(handles=[mpatches.Patch(color=palette[dim], label=str(dim)) for dim in dimensions])

        plt.title('Persistence diagram')
        plt.xlabel('Birth')
        plt.ylabel('Death')
        # Ends plot on infinity value and starts a little bit before min_birth
        plt.axis([axis_start, infinity, axis_start, infinity + delta])
        return plt

    try:
        from scipy.stats import kde
        import math

        def plot_persistence_density(persistence=[], persistence_file='',
                                     nbins=300, bw_method=None,
                                     max_plots=1000, dimension=None,
                                     cmap=plt.cm.hot_r, legend=False):
            """This function plots the persistence density from persistence
            values list or from a :doc:`persistence file <fileformats>`. Be
            aware that this function does not distinguish the dimension, it is
            up to you to select the required one.

            :param persistence: Persistence values list.
            :type persistence: list of tuples(dimension, tuple(birth, death)).
            :param persistence_file: A :doc:`persistence file <fileformats>`
                style name (reset persistence if both are set).
            :type persistence_file: string
            :param nbins: Evaluate a gaussian kde on a regular grid of nbins x
                nbins over data extents (default is 300)
            :type nbins: int.
            :param bw_method: The method used to calculate the estimator
                bandwidth. This can be 'scott', 'silverman', a scalar constant
                or a callable. If a scalar, this will be used directly as
                kde.factor. If a callable, it should take a gaussian_kde
                instance as only parameter and return a scalar. If None
                (default), 'scott' is used. See scipy.stats.gaussian_kde
                documentation for more details.
            :type bw_method: str, scalar or callable, optional.
            :param max_intervals: maximal number of intervals to display.
                Selected points are those with the longest life time. Set it
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
            :returns: A matplotlib object containing diagram plot of persistence
                (launch `show()` method on it to display it).
            """
            if persistence_file is not '':
                if os.path.isfile(persistence_file):
                    # Reset persistence
                    persistence = []
                    diag = read_persistence_intervals_grouped_by_dimension(persistence_file=persistence_file)
                    for key in diag.keys():
                        for persistence_interval in diag[key]:
                            persistence.append((key, persistence_interval))
                else:
                    print("file " + persistence_file + " not found.")
                    return None

            persistence_dim = []
            if dimension is not None:
                persistence_dim = [(dim_interval) for dim_interval in persistence if (dim_interval[0] == dimension)]
            else:
                persistence_dim = persistence

            if max_intervals > 0 and max_intervals < len(persistence_dim):
                # Sort by life time, then takes only the max_intervals elements
                persistence_dim = sorted(persistence_dim,
                                     key=lambda life_time: life_time[1][1]-life_time[1][0],
                                     reverse=True)[:max_intervals]

            # Set as numpy array birth and death (remove undefined values - inf and NaN)
            birth = np.asarray([(interval[1][0]) for interval in persistence_dim if (math.isfinite(interval[1][1]) and math.isfinite(interval[1][0]))])
            death = np.asarray([(interval[1][1]) for interval in persistence_dim if (math.isfinite(interval[1][1]) and math.isfinite(interval[1][0]))])

            # line display of equation : birth = death
            x = np.linspace(death.min(), birth.max(), 1000)
            plt.plot(x, x, color='k', linewidth=1.0)

            # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
            k = kde.gaussian_kde([birth,death], bw_method=bw_method)
            xi, yi = np.mgrid[birth.min():birth.max():nbins*1j, death.min():death.max():nbins*1j]
            zi = k(np.vstack([xi.flatten(), yi.flatten()]))

            # Make the plot
            plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=cmap)

            if legend:
                plt.colorbar()

            plt.title('Persistence density')
            plt.xlabel('Birth')
            plt.ylabel('Death')
            return plt

    except ImportError:
        # Continue in case of import error, functions won't be available
        pass

except ImportError:
    # Continue in case of import error, functions won't be available
    pass
