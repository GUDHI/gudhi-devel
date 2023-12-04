# This file is part of the Gudhi Library - https://gudhi.inria.fr/ -
# which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full
# license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - 2023/11 Vincent Rouvreau: numpy interface for read_points_from_off_file
#   - YYYY/MM Author: Description of the modification

from __future__ import print_function
from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.string cimport string
cimport cython
import errno
import os
import numpy as np

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

cdef extern from "Off_reader_interface.h" namespace "Gudhi":
    vector[vector[double]] read_points_from_OFF_file(string off_file)

def read_points_from_off_file(off_file=''):
    """Read points from an `OFF file <fileformats.html#off-file-format>`_.

    :param off_file: An OFF file style name.
    :type off_file: string

    :returns:  The point set.
    :rtype: numpy.ndarray
    """
    if off_file:
        with open(off_file) as input_file:
            line = next(input_file).lower()
            # First line should be "OFF" (3d case) or "nOFF" (dD case)
            if line.startswith("off"):
                index = 0
            elif line.startswith("noff"):
                index = 1
            else:
                raise ValueError(f"Inconsistent OFF header for {off_file}, got '{line.rstrip()}', should be 'OFF' or 'nOFF'")
            # Number of points is the first number ("OFF" case) or the second one ("nOFF" case) of the second line
            nb_points = int(next(input_file).split()[index])
        
        return np.loadtxt(off_file, skiprows=2, max_rows=nb_points, dtype=np.floating)

@cython.embedsignature(True)
def write_points_to_off_file(fname, points):
    """Write points to an `OFF file <fileformats.html#off-file-format>`_.

    A simple wrapper for `numpy.savetxt`.

    :param fname: Name of the OFF file.
    :type fname: str or file handle
    :param points: Point coordinates.
    :type points: numpy array of shape (n, dim)
    """
    points = np.array(points, copy=False)
    assert len(points.shape) == 2
    dim = points.shape[1]
    if dim == 3:
        head = 'OFF\n{} 0 0'.format(points.shape[0])
    else:
        head = 'nOFF\n{} {} 0 0'.format(dim, points.shape[0])
    np.savetxt(fname, points, header=head, comments='')
