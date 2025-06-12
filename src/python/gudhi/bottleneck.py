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
__license__ = "GPL v3"

from collections.abc import Sequence
import warnings
import numpy as np
from numpy.typing import ArrayLike

from gudhi import _bottleneck_ext as t


def _diagram_as_numpy_array(diagram: ArrayLike) -> np.ndarray:
    dgm = np.asarray(diagram, dtype=np.double)
    # to allow empty diagrams (invalid shapes will be deled with by nanobind later)
    if dgm.size == 0 and dgm.ndim == 1:
        dgm = np.empty((0, 2))
    return dgm


def bottleneck_distance(diagram_1: ArrayLike, diagram_2: ArrayLike, e: float = None) -> float:
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
    return t._bottleneck_distance(
        _diagram_as_numpy_array(diagram_1), _diagram_as_numpy_array(diagram_2), e
    )
