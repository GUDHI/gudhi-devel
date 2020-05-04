# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair
import os

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "GPL v3"

cdef extern from "Bottleneck_distance_interface.h" namespace "Gudhi::persistence_diagram":
    double bottleneck(vector[pair[double, double]], vector[pair[double, double]], double) nogil
    double bottleneck(vector[pair[double, double]], vector[pair[double, double]]) nogil

def bottleneck_distance(diagram_1, diagram_2, e=None):
    """This function returns the point corresponding to a given vertex.

    :param diagram_1: The first diagram.
    :type diagram_1: vector[pair[double, double]]
    :param diagram_2: The second diagram.
    :type diagram_2: vector[pair[double, double]]
    :param e: If `e` is 0, this uses an expensive algorithm to compute the
        exact distance.
        If `e` is not 0, it asks for an additive `e`-approximation, and
        currently also allows a small multiplicative error (the last 2 or 3
        bits of the mantissa may be wrong). This version of the algorithm takes
        advantage of the limited precision of `double` and is usually a lot
        faster to compute, whatever the value of `e`.

        Thus, by default, `e` is the smallest positive double.
    :type e: float
    :rtype: float
    :returns: the bottleneck distance.
    """
    cdef vector[pair[double, double]] dgm1 = diagram_1
    cdef vector[pair[double, double]] dgm2 = diagram_2
    cdef double eps
    cdef double ret
    if e is None:
        with nogil:
            # Default value is the smallest double value (not 0, 0 is for exact version)
            ret = bottleneck(dgm1, dgm2)
    else:
        eps = e
        with nogil:
            # Can be 0 for exact version
            ret = bottleneck(dgm1, dgm2, eps)
    return ret
