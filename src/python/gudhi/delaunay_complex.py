# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):  Thibaud Kloczko
#
# Copyright (C) 2025 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from gudhi import delaunay_complex_ext as t

class DelaunayComplex(t.Delaunay_complex_interface):
    """DelaunayComplex is a simplicial complex constructed from the finite cells of a Delaunay Triangulation.

    When :paramref:`~gudhi.DelaunayComplex.create_simplex_tree.filtration` is:

    * `None` (default value) - The filtration value of each simplex is not computed (set to `NaN`)
    * `'alpha'`              - The filtration value of each simplex is computed as an :class:`~gudhi.AlphaComplex`
    * `'cech'`               - The filtration value of each simplex is computed as a :class:`~gudhi.DelaunayCechComplex`
    """

    def __init__(self, points=[], weights=None, precision='safe'):
        """DelaunayComplex constructor.

        :param points: A list of points in d-Dimension.
        :type points: Iterable[Iterable[float]]

        :param weights: A list of weights. If set, the number of weights must correspond to the number of points.
        :type weights: Iterable[float]

        :param precision: Delaunay complex precision can be 'fast', 'safe' or 'exact'. Default is 'safe'.
        :type precision: string

        :raises ValueError: In case of inconsistency between the number of points and weights.
        """
        assert precision in ['fast', 'safe', 'exact'], "Delaunay complex precision can only be 'fast', 'safe' or 'exact'"
        fast = precision == 'fast'
        exact = precision == 'exact'

        # weights are set but is inconsistent with the number of points
        if weights is not None and len(weights) != len(points):
            raise ValueError("Inconsistency between the number of points and weights")
        else:
            weights = []

        super().__init__(points, weights, fast, exact)


#
# delaunay_complex.py ends here
