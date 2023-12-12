# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Theo Lacombe
#
# Copyright (C) 2021 Université Gustave Eiffel
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

import numpy as np
import scipy.spatial.distance as sc


def _dist_to_diag(X, internal_p):
    """
    :param X: (n x 2) array encoding the points of a persistent diagram.
    :param internal_p: Ground metric (i.e. norm L^p).
    :returns: (n) array encoding the (respective orthogonal) distances of the points to the diagonal

    .. note::
        Assumes that the points are above the diagonal.
    """
    return (X[:, 1] - X[:, 0]) * 2 ** (1.0 / internal_p - 1)


def _build_dist_matrix(X, Y, order, internal_p):
    """
    :param X: (n x 2) numpy.array encoding the (points of the) first diagram.
    :param Y: (m x 2) numpy.array encoding the second diagram.
    :param order: exponent for the Wasserstein metric.
    :param internal_p: Ground metric (i.e. norm L^p).
    :returns: (n+1) x (m+1) np.array encoding the cost matrix C.
                For 0 <= i < n, 0 <= j < m, C[i,j] encodes the distance between X[i] and Y[j],
                while C[i, m] (resp. C[n, j]) encodes the distance (to the p) between X[i] (resp Y[j])
                and its orthogonal projection onto the diagonal.
                note also that C[n, m] = 0  (it costs nothing to move from the diagonal to the diagonal).
    """
    Cxd = _dist_to_diag(X, internal_p)**order
    Cdy = _dist_to_diag(Y, internal_p)**order
    if np.isinf(internal_p):
        C = sc.cdist(X,Y, metric='chebyshev')**order
    else:
        C = sc.cdist(X,Y, metric='minkowski', p=internal_p)**order
    Cf = np.hstack((C, Cxd[:,None]))
    Cdy = np.append(Cdy, 0)

    Cf = np.vstack((Cf, Cdy[None,:]))

    return Cf