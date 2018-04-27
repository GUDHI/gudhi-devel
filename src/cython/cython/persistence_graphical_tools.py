import matplotlib.pyplot as plt
import numpy as np
import os

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Vincent Rouvreau, Bertrand Michel

   Copyright (C) 2016 INRIA

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
__copyright__ = "Copyright (C) 2016 INRIA"
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
 
def show_palette_values(alpha=0.6):
    """This function shows palette color values in function of the dimension.

    :param alpha: alpha value in [0.0, 1.0] for horizontal bars (default is 0.6).
    :type alpha: float.
    :returns: plot the dimension palette values.
    """
    colors = []
    for color in palette:
        colors.append(color)

    y_pos = np.arange(len(palette))

    plt.barh(y_pos, y_pos + 1, align='center', alpha=alpha, color=colors)
    plt.ylabel('Dimension')
    plt.title('Dimension palette values')
    return plt

def plot_persistence_barcode(persistence=[], persistence_file='', alpha=0.6, max_barcodes=1000):
    """This function plots the persistence bar code.

    :param persistence: The persistence to plot.
    :type persistence: list of tuples(dimension, tuple(birth, death)).
    :param persistence_file: A persistence file style name (reset persistence if both are set).
    :type persistence_file: string
    :param alpha: barcode transparency value (0.0 transparent through 1.0 opaque - default is 0.6).
    :type alpha: float.
    :param max_barcodes: number of maximal barcodes to be displayed.
        Set it to 0 to see all, Default value is 1000.
        (persistence will be sorted by life time if max_barcodes is set)
    :type max_barcodes: int.
    :returns: plot -- An horizontal bar plot of persistence.
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
    delta = ((max_death - min_birth) / 10.0)
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

    plt.title('Persistence barcode')
    # Ends plot on infinity value and starts a little bit before min_birth
    plt.axis([axis_start, infinity, 0, ind])
    return plt

def plot_persistence_diagram(persistence=[], persistence_file='', alpha=0.6, band=0., max_plots=1000):
    """This function plots the persistence diagram with an optional confidence band.

    :param persistence: The persistence to plot.
    :type persistence: list of tuples(dimension, tuple(birth, death)).
    :param persistence_file: A persistence file style name (reset persistence if both are set).
    :type persistence_file: string
    :param alpha: plot transparency value (0.0 transparent through 1.0 opaque - default is 0.6).
    :type alpha: float.
    :param band: band (not displayed if :math:`\leq` 0. - default is 0.)
    :type band: float.
    :param max_plots: number of maximal plots to be displayed
        Set it to 0 to see all, Default value is 1000.
        (persistence will be sorted by life time if max_plots is set)
    :type max_plots: int.
    :returns: plot -- A diagram plot of persistence.
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
    ind = 0
    delta = ((max_death - min_birth) / 10.0)
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
        ind = ind + 1

    plt.title('Persistence diagram')
    plt.xlabel('Birth')
    plt.ylabel('Death')
    # Ends plot on infinity value and starts a little bit before min_birth
    plt.axis([axis_start, infinity, axis_start, infinity + delta])
    return plt
