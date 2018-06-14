from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair
import os

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Mathieu Carriere

   Copyright (C) 2018 INRIA

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

__author__ = "Mathieu Carriere"
__copyright__ = "Copyright (C) 2018 INRIA"
__license__ = "GPL v3"

cdef extern from "Kernels_interface.h" namespace "Gudhi::persistence_diagram":
    double                 sw            (vector[pair[double, double]],          vector[pair[double, double]],          double, int)
    vector[vector[double]] sw_matrix     (vector[vector[pair[double, double]]],  vector[vector[pair[double, double]]],  double, int)
    double                 pss           (vector[pair[double, double]],          vector[pair[double, double]],          double, int)
    vector[vector[double]] pss_matrix    (vector[vector[pair[double, double]]],  vector[vector[pair[double, double]]],  double, int)
    double                 pwg           (vector[pair[double, double]],          vector[pair[double, double]],          int, string, double, double, double)
    vector[vector[double]] pwg_matrix    (vector[vector[pair[double, double]]],  vector[vector[pair[double, double]]],  int, string, double, double, double)

def sliced_wasserstein(diagram_1, diagram_2, sigma = 1, N = 100):
    """

    :param diagram_1: The first diagram.
    :type diagram_1: vector[pair[double, double]]
    :param diagram_2: The second diagram.
    :type diagram_2: vector[pair[double, double]]
    :param sigma: bandwidth of Gaussian
    :param N: number of directions

    :returns: the sliced wasserstein kernel.
    """
    return sw(diagram_1, diagram_2, sigma, N)

def sliced_wasserstein_matrix(diagrams_1, diagrams_2, sigma = 1, N = 100):
    """

    :param diagram_1: The first set of diagrams.
    :type diagram_1: vector[vector[pair[double, double]]]
    :param diagram_2: The second set of diagrams.
    :type diagram_2: vector[vector[pair[double, double]]]
    :param sigma: bandwidth of Gaussian
    :param N: number of directions

    :returns: the sliced wasserstein kernel matrix.
    """
    return sw_matrix(diagrams_1, diagrams_2, sigma, N)

def persistence_weighted_gaussian(diagram_1, diagram_2, N = 100, weight = "arctan", sigma = 1.0, C = 1.0, p = 1.0):
    """

    :param diagram_1: The first diagram.
    :type diagram_1: vector[pair[double, double]]
    :param diagram_2: The second diagram.
    :type diagram_2: vector[pair[double, double]]
    :param N: number of Fourier features
    :param weight: weight to use for the diagram points
    :param sigma: bandwidth of Gaussian
    :param C: cost of arctan persistence weight
    :param p: power of arctan persistence weight

    :returns: the persistence weighted gaussian kernel.
    """
    return pwg(diagram_1, diagram_2, N, weight, sigma, C, p)

def persistence_weighted_gaussian_matrix(diagrams_1, diagrams_2, N = 100, weight = "arctan", sigma = 1.0, C = 1.0, p = 1.0):
    """

    :param diagram_1: The first set of diagrams.
    :type diagram_1: vector[vector[pair[double, double]]]
    :param diagram_2: The second set of diagrams.
    :type diagram_2: vector[vector[pair[double, double]]]
    :param N: number of Fourier features
    :param weight: weight to use for the diagram points
    :param sigma: bandwidth of Gaussian
    :param C: cost of arctan persistence weight
    :param p: power of arctan persistence weight

    :returns: the persistence weighted gaussian kernel matrix.
    """
    return pwg_matrix(diagrams_1, diagrams_2, N, weight, sigma, C, p)

def persistence_scale_space(diagram_1, diagram_2, sigma = 1, N = 100):
    """

    :param diagram_1: The first diagram.
    :type diagram_1: vector[pair[double, double]]
    :param diagram_2: The second diagram.
    :type diagram_2: vector[pair[double, double]]
    :param sigma: bandwidth of Gaussian
    :param N: number of Fourier features

    :returns: the persistence scale space kernel.
    """
    return pss(diagram_1, diagram_2, sigma, N)

def persistence_scale_space_matrix(diagrams_1, diagrams_2, sigma = 1, N = 100):
    """

    :param diagram_1: The first set of diagrams.
    :type diagram_1: vector[vector[pair[double, double]]]
    :param diagram_2: The second set of diagrams.
    :type diagram_2: vector[vector[pair[double, double]]]
    :param sigma: bandwidth of Gaussian
    :param N: number of Fourier features

    :returns: the persistence scale space kernel matrix.
    """
    return pss_matrix(diagrams_1, diagrams_2, sigma, N)
