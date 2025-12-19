# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Hannah Schreiber
#
# Copyright (C) 2025 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

__license__ = "BSD 3-Clause"

import numpy as np
from numpy.typing import ArrayLike
import warnings

from gudhi import _hera_ext as t


def _diagram_as_numpy_array(diagram: ArrayLike) -> np.ndarray:
    dgm = np.asarray(diagram, dtype=np.double)
    # to allow empty diagrams (invalid shapes will be deled with by nanobind later)
    if dgm.size == 0 and dgm.ndim == 1:
        dgm = np.empty((0, 2))
    return dgm


def bottleneck_distance(X: ArrayLike, Y: ArrayLike, delta: float = 0.01) -> float:
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
    return t._bottleneck_distance(
        _diagram_as_numpy_array(X), _diagram_as_numpy_array(Y), delta
    )


def wasserstein_distance(
    X: ArrayLike,
    Y: ArrayLike,
    order: float = 1,
    internal_p: float = float("Inf"),
    delta: float = 0.01,
    matching: bool = False,
):
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
            
            .. warning::
            
                For matching request, please consider using :func:`~gudhi.wasserstein.wasserstein_distance` (POT
                version) instead. This version with ``matching=True`` is known to have bugs.

    Returns:
        float|Tuple[float,numpy.array|None]: Approximate Wasserstein distance W_q(X,Y), and optionally the 
            corresponding matching
    
    """
    if matching:
        warnings.warn(
            """
            There is known bug (https://github.com/GUDHI/gudhi-devel/issues/1158) when `matching` is set to `True` with
            the Hera backend. For the moment, we recommend using `gudhi.wasserstein.wasserstein_distance` instead.
            """, UserWarning)
        
    return t._wasserstein_distance(
        _diagram_as_numpy_array(X),
        _diagram_as_numpy_array(Y),
        order,
        internal_p,
        delta,
        matching,
    )
