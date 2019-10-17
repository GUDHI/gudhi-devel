import numpy as np
import scipy.spatial.distance as sc
try:
    import ot
except ImportError:
    print("POT (Python Optimal Transport) package is not installed. Try to run $ conda install -c conda-forge pot ; or $ pip install POT")

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Theo Lacombe

    Copyright (C) 2019 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

def _proj_on_diag(X):
    '''
    :param X: (n x 2) array encoding the points of a persistent diagram.
    :returns: (n x 2) arary encoding the (respective orthogonal) projections of the points onto the diagonal
    '''
    Z = (X[:,0] + X[:,1]) / 2.
    return np.array([Z , Z]).T


def _build_dist_matrix(X, Y, p=2., q=2.):
    '''
    :param X: (n x 2) numpy.array encoding the (points of the) first diagram.
    :param Y: (m x 2) numpy.array encoding the second diagram.
    :param q: Ground metric (i.e. norm l_q).
    :param p: exponent for the Wasserstein metric.
    :returns: (n+1) x (m+1) np.array encoding the cost matrix C. 
                For 1 <= i <= n, 1 <= j <= m, C[i,j] encodes the distance between X[i] and Y[j], while C[i, m+1] (resp. C[n+1, j]) encodes the distance (to the p) between X[i] (resp Y[j]) and its orthogonal proj onto the diagonal.
                note also that C[n+1, m+1] = 0  (it costs nothing to move from the diagonal to the diagonal).
    '''
    Xdiag = _proj_on_diag(X)
    Ydiag = _proj_on_diag(Y)
    if np.isinf(q):
        C = sc.cdist(X,Y, metric='chebyshev')**p
        Cxd = np.linalg.norm(X - Xdiag, ord=q, axis=1)**p
        Cdy = np.linalg.norm(Y - Ydiag, ord=q, axis=1)**p
    else:
        C = sc.cdist(X,Y, metric='minkowski', p=q)**p
        Cxd = np.linalg.norm(X - Xdiag, ord=q, axis=1)**p
        Cdy = np.linalg.norm(Y - Ydiag, ord=q, axis=1)**p
    Cf = np.hstack((C, Cxd[:,None]))
    Cdy = np.append(Cdy, 0)

    Cf = np.vstack((Cf, Cdy[None,:]))

    return Cf


def _perstot(X, p, q):
    '''
    :param X: (n x 2) numpy.array (points of a given diagram).
    :param q: Ground metric on the (upper-half) plane (i.e. norm l_q in R^2); Default value is 2 (Euclidean norm).
    :param p: exponent for Wasserstein; Default value is 2.
    :returns: float, the total persistence of the diagram (that is, its distance to the empty diagram).    
    '''
    Xdiag = _proj_on_diag(X)
    return (np.sum(np.linalg.norm(X - Xdiag, ord=q, axis=1)**p))**(1./p)


def wasserstein_distance(X, Y, p=2., q=2.):
    '''
    :param X: (n x 2) numpy.array encoding the (finite points of the) first diagram. Must not contain essential points (i.e. with infinite coordinate).
    :param Y: (m x 2) numpy.array encoding the second diagram.
    :param q: Ground metric on the (upper-half) plane (i.e. norm l_q in R^2); Default value is 2 (euclidean norm).
    :param p: exponent for Wasserstein; Default value is 2.
    :returns: the p-Wasserstein distance (1 <= p < infty) with respect to the q-norm as ground metric.
    :rtype: float
    '''
    n = len(X)
    m = len(Y)

    # handle empty diagrams
    if X.size == 0:
        if Y.size == 0:
            return 0.
        else:
            return _perstot(Y, p, q)
    elif Y.size == 0:
        return _perstot(X, p, q)

    M = _build_dist_matrix(X, Y, p=p, q=q)
    a = np.full(n+1, 1. / (n + m) )  # weight vector of the input diagram. Uniform here.
    a[-1] = a[-1] * m                # normalized so that we have a probability measure, required by POT
    b = np.full(m+1, 1. / (n + m) )  # weight vector of the input diagram. Uniform here.
    b[-1] = b[-1] * n                # so that we have a probability measure, required by POT

    # Comptuation of the otcost using the ot.emd2 library.
    # Note: it is the squared Wasserstein distance.
    ot_cost = (n+m) * ot.emd2(a, b, M)

    return ot_cost ** (1./p)

