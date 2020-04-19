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

def _dist_to_diag(X, internal_p):
    '''
    :param X: (n x 2) array encoding the points of a persistent diagram.
    :param internal_p: Ground metric (i.e. norm L^p).
    :returns: (n) array encoding the (respective orthogonal) distances of the points to the diagonal

    .. note::
        Assumes that the points are above the diagonal.
    '''
    return (X[:, 1] - X[:, 0]) * 2 ** (1.0 / internal_p - 1)


def _build_dist_matrix(X, Y, order, internal_p):
    '''
    :param X: (n x 2) numpy.array encoding the (points of the) first diagram.
    :param Y: (m x 2) numpy.array encoding the second diagram.
    :param order: exponent for the Wasserstein metric.
    :param internal_p: Ground metric (i.e. norm L^p).
    :returns: (n+1) x (m+1) np.array encoding the cost matrix C.
                For 0 <= i < n, 0 <= j < m, C[i,j] encodes the distance between X[i] and Y[j],
                while C[i, m] (resp. C[n, j]) encodes the distance (to the p) between X[i] (resp Y[j])
                and its orthogonal projection onto the diagonal.
                note also that C[n, m] = 0  (it costs nothing to move from the diagonal to the diagonal).
    '''
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


def _perstot_autodiff(X, order, internal_p):
    '''
    Version of _perstot that works on eagerpy tensors.
    '''
    return _dist_to_diag(X, internal_p).norms.lp(order)

def _perstot(X, order, internal_p, enable_autodiff):
    '''
    :param X: (n x 2) numpy.array (points of a given diagram).
    :param order: exponent for Wasserstein. Default value is 2.
    :param internal_p: Ground metric on the (upper-half) plane (i.e. norm L^p in R^2); Default value is 2 (Euclidean norm).
    :param enable_autodiff: If X is torch.tensor, tensorflow.Tensor or jax.numpy.ndarray, make the computation
        transparent to automatic differentiation.
    :type enable_autodiff: bool
    :returns: float, the total persistence of the diagram (that is, its distance to the empty diagram).
    '''
    if enable_autodiff:
        import eagerpy as ep

        return _perstot_autodiff(ep.astensor(X), order, internal_p).raw
    else:
        return np.linalg.norm(_dist_to_diag(X, internal_p), ord=order)


def wasserstein_distance(X, Y, matching=False, order=2., internal_p=2., enable_autodiff=False):
    '''
    :param X: (n x 2) numpy.array encoding the (finite points of the) first diagram. Must not contain essential points
                (i.e. with infinite coordinate).
    :param Y: (m x 2) numpy.array encoding the second diagram.
    :param matching: if True, computes and returns the optimal matching between X and Y, encoded as
                     a (n x 2) np.array  [...[i,j]...], meaning the i-th point in X is matched to
                     the j-th point in Y, with the convention (-1) represents the diagonal.
    :param order: exponent for Wasserstein; Default value is 2.
    :param internal_p: Ground metric on the (upper-half) plane (i.e. norm L^p in R^2);
                       Default value is 2 (Euclidean norm).
    :param enable_autodiff: If X and Y are torch.tensor, tensorflow.Tensor or jax.numpy.ndarray, make the computation
        transparent to automatic differentiation.
    :type enable_autodiff: bool
    :returns: the Wasserstein distance of order q (1 <= q < infinity) between persistence diagrams with
              respect to the internal_p-norm as ground metric.
              If matching is set to True, also returns the optimal matching between X and Y.
    '''
    n = len(X)
    m = len(Y)

    # handle empty diagrams
    if n == 0:
        if m == 0:
            if not matching:
                # What if enable_autodiff?
                return 0.
            else:
                return 0., np.array([])
        else:
            if not matching:
                return _perstot(Y, order, internal_p, enable_autodiff)
            else:
                return _perstot(Y, order, internal_p, enable_autodiff), np.array([[-1, j] for j in range(m)])
    elif m == 0:
        if not matching:
            return _perstot(X, order, internal_p, enable_autodiff)
        else:
            return _perstot(X, order, internal_p, enable_autodiff), np.array([[i, -1] for i in range(n)])

    if enable_autodiff:
        import eagerpy as ep

        X_orig = ep.astensor(X)
        Y_orig = ep.astensor(Y)
        X = X_orig.numpy()
        Y = Y_orig.numpy()
    M = _build_dist_matrix(X, Y, order=order, internal_p=internal_p)
    a = np.ones(n+1) # weight vector of the input diagram. Uniform here.
    a[-1] = m
    b = np.ones(m+1) # weight vector of the input diagram. Uniform here.
    b[-1] = n

    if matching:
        assert not enable_autodiff, "matching and enable_autodiff are currently incompatible"
        P = ot.emd(a=a,b=b,M=M, numItermax=2000000)
        ot_cost = np.sum(np.multiply(P,M))
        P[-1, -1] = 0  # Remove matching corresponding to the diagonal
        match = np.argwhere(P)
        # Now we turn to -1 points encoding the diagonal
        match[:,0][match[:,0] >= n] = -1
        match[:,1][match[:,1] >= m] = -1
        return ot_cost ** (1./order) , match

    if enable_autodiff:
        P = ot.emd(a=a, b=b, M=M, numItermax=2000000)
        pairs = np.argwhere(P[:-1, :-1])
        diag1 = np.nonzero(P[:-1, -1])
        diag2 = np.nonzero(P[-1, :-1])
        dists = []
        # empty arrays are not handled properly by the helpers, so we avoid calling them
        if len(pairs):
            dists.append((Y_orig[pairs[:, 1]] - X_orig[pairs[:, 0]]).norms.lp(internal_p, axis=-1).norms.lp(order))
        if len(diag1):
            dists.append(_perstot_autodiff(X_orig[diag1], order, internal_p))
        if len(diag2):
            dists.append(_perstot_autodiff(Y_orig[diag2], order, internal_p))
        dists = [dist.reshape(1) for dist in dists]
        return ep.concatenate(dists).norms.lp(order).raw
        # Should just compute the L^order norm manually?
        # We can also concatenate the 3 vectors to compute just one norm.

    # Comptuation of the otcost using the ot.emd2 library.
    # Note: it is the Wasserstein distance to the power q.
    # The default numItermax=100000 is not sufficient for some examples with 5000 points, what is a good value?
    ot_cost = ot.emd2(a, b, M, numItermax=2000000)

    return ot_cost ** (1./order)
