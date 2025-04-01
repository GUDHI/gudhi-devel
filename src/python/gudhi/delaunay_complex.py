# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):  Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

__author__ = "Vincent Rouvreau"
__maintainer__ = "Thibaud Kloczko"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

import warnings

from gudhi import _delaunay_complex_ext as t

from typing import Literal, Optional
from collections.abc import Iterable

from gudhi.simplex_tree import SimplexTree


# DelaunayComplex python interface
class DelaunayComplex(t.Delaunay_complex_interface):
    """DelaunayComplex is a simplicial complex constructed from the finite cells of a Delaunay Triangulation.

    When :paramref:`~gudhi.DelaunayComplex.create_simplex_tree.filtration` is:

    * `None` (default value) - The filtration value of each simplex is not computed (set to `NaN`)
    * `'alpha'`              - The filtration value of each simplex is computed as an :class:`~gudhi.AlphaComplex`
    * `'cech'`               - The filtration value of each simplex is computed as a :class:`~gudhi.DelaunayCechComplex`
    """

    def __init__(
        self,
        points: Iterable[Iterable[float]] = [],
        weights: Optional[Iterable[float]] = None,
        precision: Literal["fast", "safe", "exact"] = "safe",
    ):
        """Args:
              points (Iterable[Iterable[float]]): A list of points in d-Dimension.
              weights (Optional[Iterable[float]]): A list of weights. If set, the number of weights must correspond to
              the number of points.
              precision (str): Complex precision can be `'fast'`, `'safe'` or `'exact'`. Default is `'safe'`.

        :raises ValueError: In case of inconsistency between the number of points and weights.
        """
        assert precision in [
            "fast",
            "safe",
            "exact",
        ], "Delaunay complex precision can only be 'fast', 'safe' or 'exact'"
        fast = precision == "fast"
        exact = precision == "exact"

        # weights are set but is inconsistent with the number of points
        if weights is not None and len(weights) != len(points):
            raise ValueError("Inconsistency between the number of points and weights")
        else:
            if weights is None:
                weights = []

        super().__init__(points, weights, fast, exact)

    def create_simplex_tree(
        self,
        max_alpha_square: float = float("inf"),
        filtration: Optional[Literal["alpha", "cech"]] = None,
        output_squared_values: bool = True,
    ) -> SimplexTree:
        """
        Args:
            max_alpha_square (float): The maximum alpha square threshold the simplices shall not exceed. Default is set to
                                      infinity, and there is very little point using anything else since it does not save time.
            filtration (Literal): Set this value to `None` (default value) if filtration values are not needed to be computed
                                  (will be set to `NaN`). Set it to `alpha` to compute the filtration values with the Alpha complex, or
                                  to `cech` to compute the Delaunay Cech complex.
            output_squared_values (bool): Square filtration values when `True`. Default is `True`.
        Returns:
            SimplexTree: A simplex tree created from the Delaunay Triangulation. The vertex `k` corresponds to the k-th
                         input point. The vertices may not be numbered contiguously as some points may be discarded in the
                         triangulation (duplicate points, weighted hidden point, ...).
        """
        if not filtration in [None, "alpha", "cech"]:
            raise ValueError(
                f"'{filtration}' is not a valid filtration value. Must be None, 'alpha' or 'cech'"
            )

        filt = t.Filtration.NONE
        if filtration == "cech":
            filt = t.Filtration.CECH
        elif filtration == "alpha":
            filt = t.Filtration.ALPHA

        stree = SimplexTree()
        super().create_simplex_tree(stree, max_alpha_square, filt, output_squared_values)
        return stree


class AlphaComplex(DelaunayComplex):
    """AlphaComplex is a simplicial complex constructed from the finite cells of a Delaunay Triangulation.

    The filtration value of each simplex is computed as the squared radius of the smallest empty sphere passing through
    all of its vertices.

    All simplices that have a filtration value strictly greater than a given alpha squared value are not inserted into
    the complex.

    For more details about the algorithm, please refer to the
    `Alpha complex C++ documentation <https://gudhi.inria.fr/doc/latest/group__alpha__complex.html>`_

    .. note::

        When DelaunayComplex is constructed with an infinite value of alpha, the complex is a Delaunay complex.
    """

    def create_simplex_tree(
        self,
        max_alpha_square: float = float("inf"),
        default_filtration_value: bool = False,
        output_squared_values: bool = True,
    ) -> SimplexTree:
        """
        Args:
            max_alpha_square (float): The maximum alpha square threshold the simplices shall not exceed. Default is set to
                                      infinity, and there is very little point using anything else since it does not save time.
            default_filtration_value (bool): Default value is `False` (which means compute the filtration values). Set this
                                             value to `True` if filtration values are not needed to be computed (will be set to `NaN`), but please
                                             consider constructing a :class:`~gudhi.DelaunayComplex` instead.
            output_squared_values (bool): Square filtration values when `True`. Default is `True` to keep backward
                                          compatibility.
        Returns:
            SimplexTree: A simplex tree created from the Delaunay Triangulation. The vertex `k` corresponds to the k-th
                         input point. The vertices may not be numbered contiguously as some points may be discarded in the
                         triangulation (duplicate points, weighted hidden point, ...).
        """
        filtration = "alpha"
        if default_filtration_value:
            filtration = None
            warnings.warn(
                """Since Gudhi 3.10, creating an AlphaComplex with default_filtration_value=True is deprecated.
                          Please consider constructing a DelaunayComplex instead.
                          """,
                DeprecationWarning,
            )
        return super().create_simplex_tree(max_alpha_square, filtration, output_squared_values)


class DelaunayCechComplex(DelaunayComplex):
    """DelaunayCechComplex is a simplicial complex constructed from the finite cells of a Delaunay Triangulation.

    The filtration value of each simplex is equal to the squared radius of its minimal enclosing ball (MEB).

    All simplices that have a filtration value strictly greater than a given alpha squared value are not inserted into
    the complex.

    .. note::

        When DelaunayCechComplex is constructed with an infinite value of alpha, the complex is a Delaunay complex.
    """

    def __init__(self, points=[], precision="safe"):
        """
        Args:
            points (Iterable[Iterable[float]]): A list of points in d-Dimension.
            precision (str): Complex precision can be `'fast'`, `'safe'` or `'exact'`. Default is `'safe'`.
        """
        super().__init__(points=points, weights=None, precision=precision)

    def create_simplex_tree(
        self, max_alpha_square: float = float("inf"), output_squared_values: bool = True
    ) -> SimplexTree:
        """
        Args:
            max_alpha_square (float): The maximum alpha square threshold the simplices shall not exceed. Default is set to
                                      infinity, and there is very little point using anything else since it does not save time.
            output_squared_values (bool): Square filtration values when `True`. Default is `True`.
        Returns:
            SimplexTree: A simplex tree created from the Delaunay Čech Triangulation. The vertex `k` corresponds to the
                         k-th input point. The vertices may not be numbered contiguously as some points may be discarded in the
                         triangulation(duplicate points, ...).
        """
        return super().create_simplex_tree(max_alpha_square, "cech", output_squared_values)


#
# delaunay_complex.py ends here
