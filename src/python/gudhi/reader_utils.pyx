from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.pair cimport pair

from os import path
from numpy import array as np_array

# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2017 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2017 Inria"
__license__ = "MIT"

cdef extern from "Reader_utils_interface.h" namespace "Gudhi":
    vector[vector[double]] read_matrix_from_csv_file(string off_file, char separator)
    map[int, vector[pair[double, double]]] read_pers_intervals_grouped_by_dimension(string filename)
    vector[pair[double, double]] read_pers_intervals_in_dimension(string filename, int only_this_dim)

def read_lower_triangular_matrix_from_csv_file(csv_file='', separator=';'):
    """Read lower triangular matrix from a CSV style file.

    :param csv_file: A CSV file style name.
    :type csv_file: string
    :param separator: The value separator in the CSV file. Default value is ';'
    :type separator: char

    :returns:  The lower triangular matrix.
    :rtype: List[List[float]]
    """
    if csv_file:
        if path.isfile(csv_file):
            return read_matrix_from_csv_file(str.encode(csv_file), ord(separator[0]))
    print("file " + csv_file + " not set or not found.")
    return []

def read_persistence_intervals_grouped_by_dimension(persistence_file=''):
    """Reads a file containing persistence intervals.
    Each line might contain 2, 3 or 4 values: [[field] dimension] birth death
    The return value is a `dict(dim, list(tuple(birth, death)))`
    where `dim` is an `int`, `birth` a `float`, and `death` a `float`.
    Note: the function does not check that birth <= death.

    :param persistence_file: A persistence file style name.
    :type persistence_file: string

    :returns:  The persistence pairs grouped by dimension.
    :rtype: Dict[int, List[Tuple[float, float]]]
    """
    if persistence_file:
        if path.isfile(persistence_file):
            return read_pers_intervals_grouped_by_dimension(str.encode(persistence_file))
    print("file " + persistence_file + " not set or not found.")
    return []

def read_persistence_intervals_in_dimension(persistence_file='', only_this_dim=-1):
    """Reads a file containing persistence intervals.
    Each line of persistence_file might contain 2, 3 or 4 values:
    [[field] dimension] birth death
    Note: the function does not check that birth <= death.

    :param persistence_file: A persistence file style name.
    :type persistence_file: string
    :param only_this_dim: The specific dimension. Default value is -1.
        If `only_this_dim` = -1, dimension is ignored and all lines are returned.
        If `only_this_dim` is >= 0, only the lines where dimension =
        `only_this_dim` (or where dimension is not specified) are returned.
    :type only_this_dim: int.

    :returns:  The persistence intervals.
    :rtype: numpy array of dimension 2
    """
    if persistence_file:
        if path.isfile(persistence_file):
            return np_array(read_pers_intervals_in_dimension(str.encode(
                persistence_file), only_this_dim))
    print("file " + persistence_file + " not set or not found.")
    return []
