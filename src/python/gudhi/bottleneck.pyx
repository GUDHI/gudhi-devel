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
cimport numpy as np
import numpy as np
import os

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "GPL v3"

cdef extern from "Bottleneck_distance_interface.h" namespace "Gudhi::persistence_diagram":
    double bottleneck(void*, int, void*, int, double) nogil
    double bottleneck(void*, int, void*, int) nogil

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
    diagram_1 = np.asarray(diagram_1, dtype=np.float64)
    diagram_2 = np.asarray(diagram_2, dtype=np.float64)
    if diagram_1.size == 0:
        diagram_1, diagram_2 = diagram_2, diagram_1
    if diagram_1.size == 0:
        return 0.
    if diagram_2.size == 0:
        return (diagram_1[:,1] - diagram_1[:,0]).max() / 2
    return _bottleneck_distance(diagram_1, diagram_2, e)

def _bottleneck_distance(np.ndarray[double, ndim=2, mode="c"] diagram_1, np.ndarray[double, ndim=2, mode="c"] diagram_2, e):
    assert diagram_1.shape[1] == 2
    assert diagram_2.shape[1] == 2
    cdef long n1 = diagram_1.shape[0]
    cdef long n2 = diagram_2.shape[0]
    cdef double eps
    cdef double ret
    if e is None:
        with nogil:
            # Default value is the smallest double value (not 0, 0 is for exact version)
            ret = bottleneck(diagram_1.data, n1, diagram_2.data, n2)
    else:
        eps = e
        with nogil:
            # Can be 0 for exact version
            ret = bottleneck(diagram_1.data, n1, diagram_2.data, n2, eps)
    return ret
