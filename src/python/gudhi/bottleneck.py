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

def bottleneck_distance(diagram_1:ArrayLike, diagram_2:ArrayLike, e:float = None) -> float:
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
    if len(diagram_1) == 0: # to allow empty diagrams, but I am not sure if every ArrayLike object works with `len`?
      dgm1 = np.empty((0,2))
    else:
      # delegates some format errors to numpy and assures single C++ format
      # a copy is unavoidable for sequences like lists etc. anyway with the C++ bindings
      dgm1 = np.asarray(diagram_1)
    if len(diagram_2) == 0:
      dgm2 = np.empty((0,2))
    else:
      dgm2 = np.asarray(diagram_2)
    
    return t._bottleneck_distance(dgm1, dgm2, e)

