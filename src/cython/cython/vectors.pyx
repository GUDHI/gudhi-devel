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

cdef extern from "Vectors_interface.h" namespace "Gudhi::persistence_diagram":
    vector[vector[double]] compute_ls    (vector[pair[double, double]], int, double, double, int)
    vector[vector[double]] compute_pim   (vector[pair[double, double]], double, double, int, double, double, int, string, double, double, double)

def landscape(diagram, nb_ls = 10, min_x = 0.0, max_x = 1.0, res_x = 100):
    """

    :param diagram: The diagram
    :type diagram: vector[pair[double, double]]
    :param nb_ls: Number of landscapes
    :param min_x: Minimum abscissa
    :param max_x: Maximum abscissa
    :param res_x: Number of samples 

    :returns: the landscape
    """
    return compute_ls(diagram, nb_ls, min_x, max_x, res_x)

def persistence_image(diagram, min_x = 0.0, max_x = 1.0, res_x = 10, min_y = 0.0, max_y = 1.0, res_y = 10, weight = "linear", sigma = 1.0, C = 1.0, p = 1.0):
    """

    :param diagram: The diagram
    :type diagram: vector[vector[pair[double, double]]]
    :param min_x: Minimum abscissa
    :param max_x: Maximum abscissa
    :param res_x: Number of abscissa pixels 
    :param min_x: Minimum ordinate
    :param max_x: Maximum ordinate
    :param res_x: Number of ordinate pixels 
    :param weight: Weight to use for the diagram points
    :param sigma: bandwidth of Gaussian
    :param C: cost of arctan persistence weight
    :param p: power of arctan persistence weight    

    :returns: the persistence image
    """
    return compute_pim(diagram, min_x, max_x, res_x, min_y, max_y, res_y, weight, sigma, C, p)
