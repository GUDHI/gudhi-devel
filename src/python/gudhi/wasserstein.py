# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Theo Lacombe
#
# Copyright (C) 2019 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

import numpy as np
import scipy.spatial.distance as sc
try:
    import ot
except ImportError:
    print("POT (Python Optimal Transport) package is not installed. Try to run $ conda install -c conda-forge pot ; or $ pip install POT")

def _proj_on_diag(X):
    '''
    :param X: (n x 2) array encoding the points of a persistent diagram.
    :returns: (n x 2) array encoding the (respective orthogonal) projections of the points onto the diagonal
    '''
    Z = (X[:,0] + X[:,1]) / 2.
    return np.array([Z , Z]).T


def _build_dist_matrix(X, Y, order=2., internal_p=2.):
    '''
    :param X: (n x 2) numpy.array encoding the (points of the) first diagram.
    :param Y: (m x 2) numpy.array encoding the second diagram.
    :param order: exponent for the Wasserstein metric.
    :param internal_p: Ground metric (i.e. norm L^p).
    :returns: (n+1) x (m+1) np.array encoding the cost matrix C. 
                For 1 <= i <= n, 1 <= j <= m, C[i,j] encodes the distance between X[i] and Y[j], while C[i, m+1] (resp. C[n+1, j]) encodes the distance (to the p) between X[i] (resp Y[j]) and its orthogonal proj onto the diagonal.
                note also that C[n+1, m+1] = 0  (it costs nothing to move from the diagonal to the diagonal).
    '''
    Xdiag = _proj_on_diag(X)
    Ydiag = _proj_on_diag(Y)
    if np.isinf(internal_p):
        C = sc.cdist(X,Y, metric='chebyshev')**order
        Cxd = np.linalg.norm(X - Xdiag, ord=internal_p, axis=1)**order
        Cdy = np.linalg.norm(Y - Ydiag, ord=internal_p, axis=1)**order
    else:
        C = sc.cdist(X,Y, metric='minkowski', p=internal_p)**order
        Cxd = np.linalg.norm(X - Xdiag, ord=internal_p, axis=1)**order
        Cdy = np.linalg.norm(Y - Ydiag, ord=internal_p, axis=1)**order
    Cf = np.hstack((C, Cxd[:,None]))
    Cdy = np.append(Cdy, 0)

    Cf = np.vstack((Cf, Cdy[None,:]))

    return Cf


def _perstot(X, order, internal_p):
    '''
    :param X: (n x 2) numpy.array (points of a given diagram).
    :param order: exponent for Wasserstein. Default value is 2.
    :param internal_p: Ground metric on the (upper-half) plane (i.e. norm L^p in R^2); Default value is 2 (Euclidean norm).
    :returns: float, the total persistence of the diagram (that is, its distance to the empty diagram).    
    '''
    Xdiag = _proj_on_diag(X)
    return (np.sum(np.linalg.norm(X - Xdiag, ord=internal_p, axis=1)**order))**(1./order)


def wasserstein_distance(X, Y, order=2., internal_p=2.):
    '''
    :param X: (n x 2) numpy.array encoding the (finite points of the) first diagram. Must not contain essential points (i.e. with infinite coordinate).
    :param Y: (m x 2) numpy.array encoding the second diagram.
    :param order: exponent for Wasserstein; Default value is 2.
    :param internal_p: Ground metric on the (upper-half) plane (i.e. norm L^p in R^2); Default value is 2 (Euclidean norm).
    :returns: the Wasserstein distance of order q (1 <= q < infinity) between persistence diagrams with respect to the internal_p-norm as ground metric.
    :rtype: float
    '''
    n = len(X)
    m = len(Y)

    # handle empty diagrams
    if X.size == 0:
        if Y.size == 0:
            return 0.
        else:
            return _perstot(Y, order, internal_p)
    elif Y.size == 0:
        return _perstot(X, order, internal_p)

    M = _build_dist_matrix(X, Y, order=order, internal_p=internal_p)
    a = np.full(n+1, 1. / (n + m) )  # weight vector of the input diagram. Uniform here.
    a[-1] = a[-1] * m                # normalized so that we have a probability measure, required by POT
    b = np.full(m+1, 1. / (n + m) )  # weight vector of the input diagram. Uniform here.
    b[-1] = b[-1] * n                # so that we have a probability measure, required by POT

    # Comptuation of the otcost using the ot.emd2 library.
    # Note: it is the Wasserstein distance to the power q.
    # The default numItermax=100000 is not sufficient for some examples with 5000 points, what is a good value?
    ot_cost = (n+m) * ot.emd2(a, b, M, numItermax=2000000)

    return ot_cost ** (1./order)
