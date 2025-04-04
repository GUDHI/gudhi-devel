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
__license__ = "BSD 3-Clause"

from collections.abc import Sequence
import warnings

from gudhi import _hera_ext as t

def bottleneck_distance(X, Y, delta = .01):
    """Compute the Bottleneck distance between two diagrams.
    Points at infinity are supported.

    .. note::
       Points on the diagonal are not supported and must be filtered out before calling this function.

    Parameters:
        X (n x 2 numpy array): First diagram
        Y (n x 2 numpy array): Second diagram
        delta (float): Relative error 1+delta

    Returns:
        float: (approximate) bottleneck distance d_B(X,Y)
    """
    if isinstance(X, list):
      warnings.warn("Using a list will produce a copy. We recommend using a numpy array instead", RuntimeWarning)
      return t._bottleneck_distance_list(X, Y, delta)
    if isinstance(X, Sequence):
      warnings.warn("Using a sequence will produce a copy. We recommend using a numpy array instead", RuntimeWarning)
      return t._bottleneck_distance_sequence(X, Y, delta)
    return t._bottleneck_distance_tensor(X, Y, delta)

def wasserstein_distance(X, Y, order = 1, internal_p = float('Inf'), delta = .01, matching = False):
    """Compute the Wasserstein distance between two diagrams.
    Points at infinity are supported.

    Parameters:
        X (n x 2 numpy array): First diagram
        Y (n x 2 numpy array): Second diagram
        order (float): Wasserstein exponent W_q
        internal_p (float): Internal Minkowski norm L^p in R^2
        delta (float): Relative error 1+delta
        matching (bool): if ``True``, computes and returns the optimal matching between X and Y, encoded as a
            (n x 2) np.array [...[i,j]...], meaning the i-th point in X is matched to the j-th point in Y, with the
            convention that (-1) represents the diagonal. If the distance between two diagrams is +inf (which happens
            if the cardinalities of essential parts differ) and the matching is requested, it will be set to ``None``
            (any matching is optimal).

        Returns:
            float|Tuple[float,numpy.array|None]: Approximate Wasserstein distance W_q(X,Y), and optionally the
                corresponding matching
    """
    if isinstance(X, list):
      warnings.warn("Using a list will produce a copy. We recommend using a numpy array instead", RuntimeWarning)
      return t._wasserstein_distance_list(X, Y, order, internal_p, delta, matching)
    if isinstance(X, Sequence):
      warnings.warn("Using a sequence will produce a copy. We recommend using a numpy array instead", RuntimeWarning)
      return t._wasserstein_distance_sequence(X, Y, order, internal_p, delta, matching)
    return t._wasserstein_distance_tensor(X, Y, order, internal_p, delta, matching)



