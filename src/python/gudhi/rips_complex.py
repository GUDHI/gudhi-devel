# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - 2025/03 Thibaud Kloczko: Use nanobind instead of Cython for python bindings.
#   - YYYY/MM Author: Description of the modification

__author__ = "Vincent Rouvreau"
__maintainer__ = "Thibaud Kloczko"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

from typing import Literal, Optional
from collections.abc import Iterable

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
        points: Iterable[Iterable[float]] = [],
        distance_matrix: Iterable[Iterable[float]] = [],
        max_edge_length: float = float("inf"),
        sparse: Optional[float] = None
    ):
        """RipsComplex constructor.

        Args:
            points (Iterable[Iterable(float)]): A list of points in d-Dimension.
        Or
            distance_matrix (Iterable[Iterable(float)]: A distance matrix (full square or lower triangular).

        And in both cases

            max_edge_length (float, optional): Maximal edge length. All edges of the graph strictly greater than `threshold` are not
                                               inserted in the graph.
            sparse (float, optional): If this is not None, it switches to building a sparse Rips and represents the approximation
                                      parameter epsilon.
        """
        super().__init__()
        if sparse is not None:
            if len(distance_matrix) == 0:
                super().init_points_sparse(points, max_edge_length, sparse)
            else:
                super().init_matrix_sparse(distance_matrix, max_edge_length, sparse)
        else:
            if len(distance_matrix) == 0:
                super().init_points(points, max_edge_length)
            else:
                super().init_matrix(distance_matrix, max_edge_length)

    def create_simplex_tree(self, max_dimension: int = 1):
        """
        Args:
            max_dimension (int): graph expansion for Rips until this given maximal dimension.
        Returns:
            SimplexTree: A simplex tree encoding the Vietoris–Rips filtration.
        """
        stree = SimplexTree()
        super().create_simplex_tree(stree, max_dimension)
        return stree
