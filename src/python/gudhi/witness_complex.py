# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - 2025/03 Vincent Rouvreau: Use nanobind instead of Cython for python bindings.
#   - 2025/04 Hannah Schreiber: Re-add possibility of tensors (numpy, torch etc.) as input.
#   - YYYY/MM Author: Description of the modification

__license__ = "MIT"


from gudhi import _witness_complex_ext as t
from gudhi.simplex_tree import SimplexTree


# WitnessComplex python interface
class WitnessComplex(t.Witness_complex_interface):
    """Constructs (weak) witness complex for a given table of nearest landmarks
    with respect to witnesses.
    """

    def __init__(self, nearest_landmark_table=None):
        """WitnessComplex constructor.

        :param nearest_landmark_table: A list of lists of nearest landmarks and their distances.
            `nearest_landmark_table[w][k]==(l,d)` means that l is the k-th nearest landmark to
            witness w, and d is the (squared) distance between l and w.
        :type nearest_landmark_table: list of list of pair of int and float
        """
        if nearest_landmark_table is None:
            super().__init__()
        else:
            super().__init__(nearest_landmark_table)

    def create_simplex_tree(self, max_alpha_square=float("inf"), limit_dimension=-1):
        """
        :param max_alpha_square: The maximum relaxation parameter.
            Default is set to infinity.
        :type max_alpha_square: float
        :param limit_dimension: Represents the maximal dimension of the simplicial complex.
            Default value (-1) means no limit.
        :type limit_dimension: int
        :returns: A simplex tree created from the Delaunay Triangulation.
        :rtype: SimplexTree
        """
        stree = SimplexTree()
        if limit_dimension != -1:
            super().create_simplex_tree(stree, max_alpha_square, limit_dimension)
        else:
            super().create_simplex_tree(stree, max_alpha_square)
        return stree
