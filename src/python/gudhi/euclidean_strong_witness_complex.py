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

__author__ = "Vincent Rouvreau"
__maintainer__ = "Thibaud Kloczko, Hannah Schreiber"
__copyright__ = "Copyright (C) 2016 Inria"
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
        Args:
            landmarks (Sequence[Sequence[float]]): A list of landmarks (in the point cloud).
            witnesses (Sequence[Sequence[float]]): The point cloud (list of list of double).
        """

        if landmarks is not None and witnesses is not None:
            super().__init__(landmarks, witnesses)

    def create_simplex_tree(
        self, max_alpha_square: float = float("inf"), limit_dimension: int = -1
    ) -> SimplexTree:
        """
        Args:
            max_alpha_square (float): The maximum alpha square threshold the simplices shall not exceed.
                                      Default is set to infinity.
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


#
# euclidean_strong_witness_complex.py ends here
