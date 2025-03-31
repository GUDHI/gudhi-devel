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
import errno
import os
import numpy as np

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

def _get_next_line(file_desc, comment='#'):
    """Return the next line that is not a comment.

    :param file_desc: An open file in read mode.
    :type file_desc: file
    :param comment: The characters or list of characters used to indicate the start of a comment.
    :type comment: string

    :returns: The next line.
    :rtype: string
    """
    while True:
        # file_desc.readline() is preferred to next(file_desc), as the second option is not compatible with
        # seek and tell methods
        line = file_desc.readline()
        if (not line.startswith(comment)) and len(line.split()) > 0:
            break
    return line

def _read_off_file_header(file_desc):
    """Return the information contained in the header of an OFF file.

    :param file_desc: An open file in read mode.
    :type file_desc: file

    :returns: The point cloud dimension, and the number of points.
    :rtype: tuple(int, int)

    Raises:
        ValueError: If the file does not respect the OFF file format.
    """
    nb_vertices = -1

    line = _get_next_line(file_desc)
    # First line should be "OFF" (3d case with some variants) "4OFF" (4d case) or "nOFF" (dD case)
    if line.lower().startswith("noff"):
        # "nOFF" case, next line is the dimension
        # can also contain nb_vertices, nb_faces nb_edges (can also be on the next line)
        line = _get_next_line(file_desc)
        digits = [int(s) for s in line.split() if s.isdigit()]
        dim = digits[0]
        if len(digits) > 1:
            nb_vertices = digits[1]
            # nb_faces =  digits[2]
            # nb_edges =  digits[3] # not used - can be ignored
            # nb_cells =  digits[4]
    elif line.lower().startswith("4off"):
        dim = 4
    # "OFF", "COFF" and "STOFF" are 3d cases - let's stick with the C++ interface
    elif line.lower().find("off") >= 0:
        dim = 3
    else:
        raise ValueError(f"Inconsistent OFF header, got '{line.rstrip()}', should be 'OFF', '4OFF' or 'nOFF'")
        
    # nb_vertices can be already set by "nOFF" case, when 'dim nb_vertices nb_faces nb_edges' on the same line
    if nb_vertices < 0:
        # Number of points is the first number ("OFF" case) or the second one ("nOFF" case) of the second line
        line = _get_next_line(file_desc)
        digits = [int(s) for s in line.split() if s.isdigit()]
        nb_vertices = digits[0]
        # nb_faces =  digits[1]
        # nb_edges =  digits[2] # not used - can be ignored
        # nb_cells =  digits[3]
    # "_get_next_line + go back to the previous line" is just a hack for comments in the most likely places
    # TODO: remove "_get_next_line + go back to the previous line" when numpy ≥ 1.23.0 will be the standard
    line = _get_next_line(file_desc)
    # Here the first line without comment is read - let's go back to the beginning of this line
    file_desc.seek(file_desc.tell() - len(line))
    return dim, nb_vertices

def read_points_from_off_file(off_file=''):
    """Read points from an `OFF file <fileformats.html#off-file-format>`_.

    :param off_file: An OFF file style name.
    :type off_file: string

    :returns:  The point set.
    :rtype: numpy.ndarray

    .. warning::
        This function is using `numpy.loadtxt <https://numpy.org/doc/stable/reference/generated/numpy.loadtxt.html>`_
        with `comments='#'` as an argument. Empty and or comment lines between the points are only supported with numpy
        &ge; 1.23.0.
    """
    # newline='' is required for Windows, otherwise end of line with '\r\n' are only detected as '\n'
    # This is required by _read_off_file_header that needs the exact length of the line (to go backward in the file reading)
    with open(off_file, newline='') as input_file:
        dim, nb_points = _read_off_file_header(input_file)
        # usecols=list(range(dim)) stands here to avoid comments at the end of line
        # or colors that can be added in RGB format after the points, the faces, ...
        points = np.loadtxt(input_file, dtype=np.float64, comments='#',
                            usecols=range(dim), max_rows=nb_points)
        assert points.shape == (nb_points, dim), f"{points.shape} is different from expected ({nb_points}, {dim})"
        return points

def write_points_to_off_file(fname, points):
    """Write points to an `OFF file <fileformats.html#off-file-format>`_.

    A simple wrapper for `numpy.savetxt`.

    :param fname: Name of the OFF file.
    :type fname: str or file handle
    :param points: Point coordinates.
    :type points: numpy array of shape (n, dim)
    """
    points = np.asarray(points)
    assert len(points.shape) == 2
    dim = points.shape[1]
    if dim == 3:
        head = 'OFF\n{} 0 0'.format(points.shape[0])
    else:
        head = 'nOFF\n{} {} 0 0'.format(dim, points.shape[0])
    np.savetxt(fname, points, header=head, comments='')
