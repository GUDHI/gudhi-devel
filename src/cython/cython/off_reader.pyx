from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.string cimport string
import os

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Vincent Rouvreau

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

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 INRIA"
__license__ = "GPL v3"

cdef extern from "Off_reader_interface.h" namespace "Gudhi":
    vector[vector[double]] read_points_from_OFF_file(string off_file)

def read_off(off_file=''):
    """Read points from OFF file.

    :param off_file: An OFF file style name.
    :type off_file: string

    :returns:  The point set.
    :rtype: vector[vector[double]]
    """
    if off_file is not '':
        if os.path.isfile(off_file):
            return read_points_from_OFF_file(off_file)
        else:
            print("file " + off_file + " not found.")

