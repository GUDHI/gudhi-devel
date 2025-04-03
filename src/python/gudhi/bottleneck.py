# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Hannah Schreiber
#
# Copyright (C) 2025 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

__author__ = "Hannah Schreiber"
__maintainer__ = ""
__copyright__ = "Copyright (C) 2025 Inria"
__license__ = "MIT"

from collections.abc import Sequence
import warnings

from gudhi import _bottleneck_ext as t

def bottleneck_distance(diagram_1, diagram_2, e = None):
    """Compute the Bottleneck distance between two diagrams.
    Points at infinity and on the diagonal are supported.

    :param diagram_1: The first diagram.
    :type diagram_1: numpy array of shape (m,2)
    :param diagram_2: The second diagram.
    :type diagram_2: numpy array of shape (n,2)
    :param e: If `e` is 0, this uses an expensive algorithm to compute the
        exact distance.
        If `e` is not 0, it asks for an additive `e`-approximation, and
        currently also allows a small multiplicative error (the last 2 or 3
        bits of the mantissa may be wrong). This version of the algorithm takes
        advantage of the limited precision of `double` and is usually a lot
        faster to compute, whatever the value of `e`.
        Thus, by default (`e=None`), `e` is the smallest positive double.
    :type e: float
    :rtype: float
    :returns: the bottleneck distance.
    """
    if isinstance(diagram_1, list):
      warnings.warn("Using a list will produce a copy. We recommend using a numpy array instead", RuntimeWarning)
      return t._bottleneck_distance_list(diagram_1, diagram_2, e)
    if isinstance(diagram_1, Sequence):
      warnings.warn("Using a sequence will produce a copy. We recommend using a numpy array instead", RuntimeWarning)
      return t._bottleneck_distance_sequence(diagram_1, diagram_2, e)
    return t._bottleneck_distance_tensor(diagram_1, diagram_2, e)

