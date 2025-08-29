# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - 2025/03 Thibaud Kloczko: Use nanobind instead of Cython for python bindings.
#   - 2025/04 Hannah Schreiber: Re-add possibility of tensors (numpy, torch etc.) as input.
#   - YYYY/MM Author: Description of the modification

__license__ = "GPL v3"


from gudhi import _euclidean_strong_witness_complex_ext as t
from gudhi.simplex_tree import SimplexTree


# EuclideanStrongWitnessComplex python interface
class EuclideanStrongWitnessComplex(t.Euclidean_strong_witness_complex_interface):
    """Constructs strong witness complex for given sets of witnesses and
    landmarks in Euclidean space.
    """

    def __init__(self, landmarks=None, witnesses=None):
        """WitnessComplex constructor.

        :param landmarks: A list of landmarks (in the point cloud).
        :type landmarks: list of list of double

        :param witnesses: The point cloud.
        :type witnesses: list of list of double
        """

        if landmarks is None or witnesses is None:
            super().__init__()
        else:
            super().__init__(landmarks, witnesses)

    def create_simplex_tree(
        self, max_alpha_square: float = float("inf"), limit_dimension: int = -1
    ) -> SimplexTree:
        """
        :param max_alpha_square: The maximum alpha square threshold the
            simplices shall not exceed. Default is set to infinity.
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


#
# euclidean_strong_witness_complex.py ends here
