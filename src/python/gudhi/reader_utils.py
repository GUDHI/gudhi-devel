# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2017 Inria
#
# Modification(s):
#   - 2025/03 Thibaud Kloczko: Use nanobind instead of Cython for python bindings.
#   - YYYY/MM Author: Description of the modification

__license__ = "MIT"


from os import path
from numpy import array as np_array

from gudhi import _reader_utils_ext as t


def read_lower_triangular_matrix_from_csv_file(csv_file="", separator=";"):
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
            return t.read_matrix_from_csv_file(csv_file, separator[0])
    print("file " + csv_file + " not set or not found.")
    return []


def read_persistence_intervals_grouped_by_dimension(persistence_file=""):
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
            return t.read_pers_intervals_grouped_by_dimension(persistence_file)
    print("file " + persistence_file + " not set or not found.")
    return []


def read_persistence_intervals_in_dimension(persistence_file="", only_this_dim=-1):
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
            return np_array(
                t.read_pers_intervals_in_dimension(persistence_file, only_this_dim)
            )
    print("file " + persistence_file + " not set or not found.")
    return []
