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
__license__ = "MIT"


from typing import Literal, Optional
from collections.abc import Sequence
from numpy.typing import ArrayLike

from gudhi import _rips_complex_ext as t
from gudhi.simplex_tree import SimplexTree


# RipsComplex python interface
class RipsComplex(t.Rips_complex_interface):
    """The data structure is a one skeleton graph, or Rips graph, containing edges when the edge length is less or
    equal to a given threshold. Edge length is computed from a user given point cloud with a given distance function,
    or a distance matrix.
    """

    def __init__(
        self,
        *,
        points: ArrayLike = [],
        distance_matrix: ArrayLike = [],
        max_edge_length: float = float("inf"),
        sparse: Optional[float] = None
    ):
        """RipsComplex constructor.

        :param points: A list of points in d-Dimension.
        :type points: Sequence[Sequence[float]] or any array like object of ndim 2 and dtype convertible to float.

        Or

        :param distance_matrix: A distance matrix (full square or lower triangular).
        :type distance_matrix: Sequence[Sequence[float]] (square or just the lower triangle) or any square array like
            object of dtype convertible to float.

        And in both cases

        :param max_edge_length: Maximal edge length. All edges of the graph strictly greater than `threshold` are not
            inserted in the graph.
        :type max_edge_length: float
        :param sparse: If this is not None, it switches to building a sparse Rips and represents the approximation
            parameter epsilon.
        :type sparse: float
        """
        if sparse is not None:
            if len(distance_matrix) == 0:
                super().__init__(points, max_edge_length, sparse, True)
            else:
                super().__init__(distance_matrix, max_edge_length, sparse, False)
        else:
            if len(distance_matrix) == 0:
                super().__init__(points, max_edge_length, True)
            else:
                super().__init__(distance_matrix, max_edge_length, False)

    def create_simplex_tree(self, max_dimension: int = 1):
        """
        :param max_dimension: graph expansion for Rips until this given maximal dimension.
        :type max_dimension: int
        :returns: A simplex tree encoding the Vietoris–Rips filtration.
        :rtype: SimplexTree
        """
        stree = SimplexTree()
        super().create_simplex_tree(stree, max_dimension)
        return stree
