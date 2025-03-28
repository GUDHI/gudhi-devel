# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

# from libcpp.vector cimport vector
# from libcpp.utility cimport pair
# from libc.stdint cimport intptr_t
#
# from gudhi.simplex_tree cimport *
# from gudhi.simplex_tree import SimplexTree

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

# StrongWitnessComplex python interface
class StrongWitnessComplex:
    """Constructs (strong) witness complex for a given table of nearest
    landmarks with respect to witnesses.
    """

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self, nearest_landmark_table=None):
        """StrongWitnessComplex constructor.

        :param nearest_landmark_table: A list of lists of nearest landmarks and their distances.
            `nearest_landmark_table[w][k]==(l,d)` means that l is the k-th nearest landmark to
            witness w, and d is the (squared) distance between l and w.
        :type nearest_landmark_table: list of list of pair of int and float
        """

    def create_simplex_tree(self, max_alpha_square = float('inf'), limit_dimension = -1):
        """
        :param max_alpha_square: The maximum relaxation parameter.
            Default is set to infinity.
        :type max_alpha_square: float
        :returns: A simplex tree created from the Delaunay Triangulation.
        :rtype: SimplexTree
        """
        stree = SimplexTree()
        return stree
