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

def _get_next_line(skiprows, file_desc, comment='#'):
    while True:
        line = next(file_desc)
        skiprows += 1
        if not line.startswith(comment):
            break
    return skiprows, line

def _read_off_file_header(off_file=''):
    # Must start at -1 as first line (aka. 0), must be at least read
    skiprows = -1
    nb_vertices = -1
    with open(off_file) as input_file:
        skiprows, line = _get_next_line(skiprows, input_file)
        # First line should be "OFF" (3d case) or "nOFF" (dD case)
        if line.lower().startswith("off"):
            dim = 3
        elif line.lower().startswith("4off"):
            dim = 4
        elif line.lower().startswith("noff"):
            # "nOFF" case, next line is the dimension
            # can also contain nb_vertices, nb_faces nb_edges (can also be on the next line)
            skiprows, line = _get_next_line(skiprows, input_file)
            digits = [int(s) for s in line.split() if s.isdigit()]
            dim = digits[0]
            if len(digits) > 1:
                nb_vertices = digits[1]
                # nb_faces =  digits[2]
                # nb_edges =  digits[3] # not used - can be ignored
                # nb_cells =  digits[4]
        else:
            raise ValueError(f"Inconsistent OFF header for {off_file}, got '{line.rstrip()}', should be 'OFF', '4OFF' or 'nOFF'")
        
        # nb_vertices can be already set by "nOFF" case, when 'dim nb_vertices nb_faces nb_edges' on the same line
        if nb_vertices < 0:
            # Number of points is the first number ("OFF" case) or the second one ("nOFF" case) of the second line
            skiprows, line = _get_next_line(skiprows, input_file)
            digits = [int(s) for s in line.split() if s.isdigit()]
            nb_vertices = digits[0]
            # nb_faces =  digits[1]
            # nb_edges =  digits[2] # not used - can be ignored
            # nb_cells =  digits[3]
        # Increment skiprows if comments between header and the first point
        skiprows, line = _get_next_line(skiprows, input_file)
        return skiprows, dim, nb_vertices

def read_points_from_off_file(off_file=''):
    """Read points from an `OFF file <fileformats.html#off-file-format>`_.

    :param off_file: An OFF file style name.
    :type off_file: string

    :returns:  The point set.
    :rtype: numpy.ndarray
    """
    skiprows, dim, nb_points = _read_off_file_header(off_file)
    print(skiprows, dim, nb_points)
    # usecols=list(range(dim)) stands here to avoid comments at the end of line
    # or colors that can be added in RGB format after the points, the faces, ...
    return np.loadtxt(off_file, dtype=np.floating, comments='#', skiprows=skiprows,
                      usecols=list(range(dim)), max_rows=nb_points)

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
