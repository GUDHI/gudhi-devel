from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.string cimport string
import os

# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

cdef extern from "Off_reader_interface.h" namespace "Gudhi":
    vector[vector[double]] read_points_from_OFF_file(string off_file)

def read_off(off_file=''):
    """Read points from OFF file.

    :param off_file: An OFF file style name.
    :type off_file: string

    :returns:  The point set.
    :rtype: vector[vector[double]]
    """
    if off_file:
        if os.path.isfile(off_file):
            return read_points_from_OFF_file(str.encode(off_file))
        else:
            print("file " + off_file + " not found.")
            return []

