import numpy as np
import scipy.spatial.distance as sc
try:
    import ot
except ImportError:
    print("POT (Python Optimal Transport) package is not installed. Try to run $ pip install POT")

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Theo Lacombe

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

def proj_on_diag(X):
    '''
    param X: (n x 2) array encoding the points of a persistent diagram.
    return: (n x 2) arary encoding the (respective orthogonal) projections of the points onto the diagonal
    '''
    Z = (X[:,0] + X[:,1]) / 2.
    return np.array([Z , Z]).T


def build_dist_matrix(X, Y, p=2., q=2.):
    '''
    param X: (n x 2) np.array encoding the (points of the) first diagram.
    param Y: (m x 2) np.array encoding the second diagram.
    param q: Ground metric (i.e. norm l_q).
    param p: exponent for the Wasserstein metric.
    return: (n+1) x (m+1) np.array encoding the cost matrix C. 
                For 1 <= i <= n, 1 <= j <= m, C[i,j] encodes the distance between X[i] and Y[j], while C[i, m+1] (resp. C[n+1, j]) encodes the distance (to the p) between X[i] (resp Y[j]) and its orthogonal proj onto the diagonal.
                note also that C[n+1, m+1] = 0  (it costs nothing to move from the diagonal to the diagonal).
    '''
    Xdiag = proj_on_diag(X)
    Ydiag = proj_on_diag(Y)
    if np.isinf(p):
        C = sc.cdist(X,Y, metric='chebyshev', p=q)**p
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


def wasserstein_distance(X, Y, p=2., q=2.):
    '''
    param X, Y: (n x 2) and (m x 2) numpy array (points of persistence diagrams)
    param q: Ground metric (i.e. norm l_q); Default value is 2 (euclidean norm).
    param p: exponent for Wasserstein; Default value is 2.
    return: float, the p-Wasserstein distance (1 <= p < infty) with respect to the q-norm as ground metric.
    '''
    M = build_dist_matrix(X, Y, p=p, q=q)
    n = len(X)
    m = len(Y)
    a = 1.0 / (n + m) * np.ones(n)  # weight vector of the input diagram. Uniform here.
    hat_a = np.append(a, m/(n+m))  # so that we have a probability measure, required by POT
    b = 1.0 / (n + m) * np.ones(m)  # weight vector of the input diagram. Uniform here.
    hat_b = np.append(b, n/(m+n))  # so that we have a probability measure, required by POT

    # Comptuation of the otcost using the ot.emd2 library.
    # Note: it is the squared Wasserstein distance.
    ot_cost = (n+m) * ot.emd2(hat_a, hat_b, M)

    return np.power(ot_cost, 1./p)

