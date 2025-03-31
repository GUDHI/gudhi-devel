# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2025 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from gudhi import _strong_witness_complex_ext as t
from gudhi.simplex_tree import SimplexTree

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2025 Inria"
__license__ = "MIT"

# StrongWitnessComplex python interface
class StrongWitnessComplex(t.Strong_witness_complex_interface):
    """Constructs (strong) witness complex for a given table of nearest
    landmarks with respect to witnesses.
    """

    def __init__(self, nearest_landmark_table=None):
        """StrongWitnessComplex constructor.
        Args:
            param nearest_landmark_table (Iterable[Iterable[Pair[float]]): A list of lists of nearest landmarks and their distances.
                                                                          `nearest_landmark_table[w][k]==(l,d)` means that l is the k-th nearest landmark to
                                                                          witness w, and d is the (squared) distance between l and w.
        """
        if nearest_landmark_table is not None:
            super().__init__(nearest_landmark_table)

    def create_simplex_tree(self, max_alpha_square: float = float('inf'), limit_dimension = -1) -> SimplexTree:
        """
        Args:
            max_alpha_square (float): The maximum relaxation parameter. Default is set to infinity.
            limit_dimension (int):
        Returns:
            SimplexTree: A simplex tree created from the Delaunay Triangulation.
        """
        stree = SimplexTree()

        if limit_dimension != -1:
            super().create_simplex_tree(stree, max_alpha_square, limit_dimension)
        else:
            super().create_simplex_tree(stree, max_alpha_square)
        return stree
