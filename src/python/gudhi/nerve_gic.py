# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2018 Inria
#
# Modification(s):
#   - 2025/03 Hannah Schreiber: Use nanobind instead of Cython for python bindings.
#   - YYYY/MM Author: Description of the modification

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2018 Inria"
__license__ = "GPL v3"

import errno
import os

from gudhi._nerve_gic_ext import _Nerve_gic_interface
from gudhi.simplex_tree import SimplexTree

# CoverComplex python interface
class CoverComplex(_Nerve_gic_interface):
    """Cover complex data structure.

    The data structure is a simplicial complex, representing a Graph Induced
    simplicial Complex (GIC) or a Nerve, and whose simplices are computed with
    a cover C of a point cloud P, which often comes from the preimages of
    intervals covering the image of a function f defined on P. These intervals
    are parameterized by their resolution (either their length or their number)
    and their gain (percentage of overlap). To compute a GIC, one also needs a
    graph G built on top of P, whose cliques with vertices belonging to
    different elements of C correspond to the simplices of the GIC.
    """

    def __init__(self):
        """CoverComplex constructor.
        """
        super().__init__()

    def create_simplex_tree(self):
        """
        :returns: A simplex tree created from the Cover complex.
        :rtype: SimplexTree
        """
        simplex_tree = SimplexTree()
        super().create_simplex_tree(simplex_tree)
        return simplex_tree

    def read_point_cloud(self, off_file):
        """Reads and stores the input point cloud from .(n)OFF file.

        :param off_file: Name of the input .OFF or .nOFF file.
        :type off_file: string

        :rtype: bool
        :returns: Read file status.
        """
        if os.path.isfile(off_file):
            return super().read_point_cloud(off_file)
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), off_file)

    def set_color_from_file(self, color_file_name):
        """Computes the function used to color the nodes of the simplicial
        complex from a file containing the function values.

        :param color_file_name: Name of the input color file.
        :type color_file_name: string
        """
        if os.path.isfile(color_file_name):
            super().set_color_from_file(color_file_name)
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), color_file_name)

    def set_cover_from_file(self, cover_file_name):
        """Creates the cover C from a file containing the cover elements of
        each point (the order has to be the same as in the input file!).

        :param cover_file_name: Name of the input cover file.
        :type cover_file_name: string
        """
        if os.path.isfile(cover_file_name):
            super().set_cover_from_file(cover_file_name)
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), cover_file_name)

    def set_function_from_file(self, func_file_name):
        """Creates the function f from a file containing the function values.

        :param func_file_name: Name of the input function file.
        :type func_file_name: string
        """
        if os.path.isfile(func_file_name):
            super().set_function_from_file(func_file_name)
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), func_file_name)

    def set_graph_from_file(self, graph_file_name):
        """Creates a graph G from a file containing the edges.

        :param graph_file_name: Name of the input graph file. The graph file
            contains one edge per line, each edge being represented by the IDs
            of its two nodes.
        :type graph_file_name: string
        """
        if os.path.isfile(graph_file_name):
            super().set_graph_from_file(graph_file_name)
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), graph_file_name)

