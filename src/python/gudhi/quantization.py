""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Theo Lacombe

    Copyright (C) 2019 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""
import numpy as np
try:
    import ot
except ImportError:
    print("POT (Python Optimal Transport) package is not installed. Try to run $ conda install -c conda-forge pot ; or $ pip install POT")
import scipy.spatial.distance as sc

from wasserstein import _proj_on_diag, _build_dist_matrix


def _dist_to_diag(X, order, internal_p):
    Xdiag = _proj_on_diag(X)
    return np.linalg.norm(X - Xdiag, axis=1)**order


def _greed_init(Y,k, nb_trial=10):
    """
    :param Y: Input diagram
    :type Y: ``np.array``
    :param k: size of the target dgm
    :type k: ``np.array``
    :param nb_trial: number of repetition of this greedy initialization
    :type nb_trial: ``np.int``
    :returns: An initialization with points "nicely spread".
    """
    n = Y.shape[0]
    X = Y[np.random.choice(n, k)]  # randomly initialize positions
    d = np.mean(sc.pdist(X))
    for _ in range(nb_trial):
        X2 = Y[np.random.choice(n,k)]
        d2 = np.mean(sc.pdist(X2))
        if d2 > d:
            d = d2
            X = X2.copy()
    return X


def _balanced_kmeans(Y, w, k, t, gamma, nb_max_iter, stopping_criterion, init=None):
    '''
    :param Y: encoding support of the input (weighted) measure (dgm)
    :type Y: (n x 2) ``np.array``
    :param w: weights on each point of the input (eg distance to diag). Must have shape (n x 1).
    :type w: ``numpy.array``
    :param k: number of centroid in the output.
    :type k: ``int``
    :param t: Learning rate, must be in (0,1).
    :type t: ``float``
    :param gamma: Parameter to use in Sinkhorn approximation of opt transport plan. If 0, exact transport is computed. 
                  Using Sinkhorn approximation faster the computations for large values of ``n``.
    :type gamma: ``float``
    :param nb_max_iter: Maximum number of interation in the process.
    :type nb_max_iter: ``int``
    :param stopping_criterion: stopping criterion of the process, measured in W2.
    :type stopping_criterion: ``float``
    :param init: If ``None``, init is randomly made on points ``k`` points of the input diagram ``Y``. 
                    Otherwise, ``init`` must be a (k x 2) ``numpy.array``
    :returns: (quantized_diagram, transport_plan), `(k x 2)` ``numpy.array`` 

    Inspired from [Fast Wasserstein barycenters ; Cuturi, Doucet, 2014]. Only true in euclidean.

    TODO: line search for t ?
    '''
    n = Y.shape[0]
    assert (Y.shape[1] == 2)
    assert (w.shape[0] == n)
    w = w/np.sum(w)  # Normalization of input mesure (to 1). This has *no* impact on quantization
                     #    since computed cells are mass-scale invariant.
    if init is None:
        X = _greed_init(Y, k)
    else:
        X = init.copy()
    a = 1./k * np.ones(k)  # Uniform weights enforced

    for i in range(nb_max_iter):
        C = sc.cdist(X,Y)**2
        if gamma > 0:  # Apply sinkhorn reg
            P = ot.bregman.sinkhorn(a, w, C, gamma)  # size k x n
        else:  # Exact computation of OT
            P = ot.emd(a, w, C)
        new_X = (1 - t) * X + t * k * np.dot(P,Y)
        # We compute the evolution of X wrt the W2 metric to check for the stopping criterion.
        # L2 distance could also be used (faster to compute, slower to converge).
        e = ot.emd2([],[], sc.cdist(X, new_X)**2)**(1./2)
        if e < stopping_criterion:
            break
        else:
            X = new_X

    return X, P


def _weight_optim(hat_a, hat_b, C, t, gamma, nb_max_iter, stopping_criterion):
    '''
    :param a: `k+1` weights on each point of the (current estimate of) quantization
    :type a: ``numpy.array``
    :param hat_b: n+1  ; weights on each point (eg distance to diag)
    :type hat_b: ``numpy.array``
    :param C: `(k+1, n+1)` cost matrix (previously computed, include the diagonal)
    :type C: ``numpy.array``
    :param t: Learning rate.
    :type t: ``float``
    :param gamma: Parameter to use in Sinkhorn approximation of optimal transport plan. If 0, exact transport is computed.
    :type gamma: ``float``
    :param nb_max_iter: nombre max iter in the process
    :type nb_max_iter: ``int``
    :param stopping_criterion: stopping criterion of the process (measured in L1)
    :type stopping_criterion: ``float``
    :returns: new weights for the current estimate.

    Inspired from [Fast Wasserstein barycenters ; Cuturi, Doucet, 2014].
    	Adapted to persistence diagrams to take the diagonal into account.
    '''
    k = len(hat_a) - 1

    for i in range(nb_max_iter):
        if gamma > 0:
            dico = ot.bregman.sinkhorn(hat_a, hat_b, C, gamma, log=True)[1]
        else:
            dico = ot.emd(hat_a, hat_b, C, log=True)[1]
        
        alpha = dico['u'][:k]
        a = hat_a[:k]
        new_a = np.multiply(a, np.exp(- t * alpha))
        new_a = new_a/(2 * np.sum(new_a))

        e = np.mean(np.abs(new_a - a))

        if e < stopping_criterion:
            break
        else:
            hat_a = np.append(new_a, 0.5)

    return hat_a


def _loc_update(hat_a, hat_b, Y, C, gamma):
    '''
    :param a: Weights on X locations (size k).
    :type a: ``numpy.array``
    :param Y: location of point in attach data (size (n x 2))
    :type Y: ``numpy.array``
    :param P: , optimal (partial) transport plan between X and Y (and the diag) (size (k+1) x (n+1))
    :type P: ``numpy.array``
    :return: updated position of the points.
    '''
    k = len(hat_a) - 1
    a = hat_a[:k]

    if gamma > 0:  # We apply sinkhorn reg
        P = ot.bregman.sinkhorn(hat_a, hat_b, C, gamma)  # size k x n
    else:  # exact computation of OT (ok for n not too large)
        P = ot.emd(hat_a, hat_b, C)

    k= P.shape[0] -1
    n = P.shape[1] -1
    Pxy = P[:k,:n]  # size k x n
    new_X = np.divide(Pxy.dot(Y), a[:,None])  # size k

    Y_mean = (Y[:,0] + Y[:,1]) / 2  # size n

    t = 1/np.sum(Pxy, axis=1) - 1/a  # size k
    new_X = new_X + np.multiply(t, Pxy.dot(Y_mean))[:,None]

    return new_X, P


def balanced_quantization(Y, k,
                          weight_function = lambda X: _dist_to_diag(X, order=2., internal_p=2.),
                          t=0.5,
                          gamma=0., nb_max_iter=100, stopping_criterion=0.001):
    """
    Compute a quantization of a diagram ``Y`` with ``k`` centroids, where each centroid defines a cluster of points
    in the input diagram with same total mass. Here, the mass is given by ``weight_function``. For instance, 
    weighting by the distance to the diagonal (default) will define clusters of same `total persistence`.

    :param Y: encoding support of the input (weighted) measure (dgm)
    :type Y: (n x 2) ``np.array``
    :param weight_function: weight function on each point of the input (default: squared distance to the diagonal).
                            Must return a ``numpy.array`` of shape (n) with non-negative entries.
    :param k: number of centroids in the output.
    :type k: ``int``
    :param t: Learning rate, must be in (0,1).
    :type t: ``float``
    :param gamma: Parameter to use in Sinkhorn approximation of opt transport plan. If 0, exact transport is computed. 
                  Using Sinkhorn approximation faster the computations for large values of ``n``.
    :type gamma: ``float``
    :param nb_max_iter: Maximum number of interation in the process.
    :type nb_max_iter: ``int``
    :param stopping_criterion: stopping criterion of the process, measured in W2.
    :type stopping_criterion: ``float``
    :returns: (quantized_diagram, transport_plan), `(k x 2)` ``np.array`` 
    """
    b = weight_function(Y)
    return _balanced_kmeans(Y, b, k, t, gamma, nb_max_iter, stopping_criterion)


def kmeans_quantization(Y, k, weight_update=False, 
                        gamma=0., 
                        nb_max_iter=100, stopping_criterion=0.001, t=10):
    """
    Compute a quantization of a diagram ``Y`` (of ``n`` points) with ``k`` centroids, by solving the following problem:
    find a diagram with ``k`` points, with a mass ``n/k`` on each point that would minimize the Wasserstein distance 
    to the first diagram (taking the diagonal into account).

    If weight_update is True, optimization is performed in an EM-style: optimizing locations, then weights,
    then locations...

    :param Y: Encoding support of the input diagram.
    :type Y: `(n x 2)` ``numpy.array``
    :param k: number of centroids
    :type k: ``int``
    :param weight_update: If true, centroids weights are optimized. Otherwise stays uniform. Default ``False``.
    :type weight_update: ``boolean``
    :param gamma: Parameter to use in Sinkhorn approximation of opt transport plan. If 0, exact transport is computed.
    :type gamma: ``float``
    :param nb_max_iter: maximum iter in the process
    :type nb_max_iter: ``int``
    :param stopping_criterion: stopping criterion of the process (measured in L2)
    :type stopping_criterion: ``float``
    :returns: (centroids_positions, centroids_weights, transport_plan)
    """
    n = Y.shape[0]
    assert (Y.shape[1] == 2)

    b = np.full(n, 1./(2 * n))  # weight vector of the input diagram. Uniform here.
    hat_b = np.append(b, 0.5)   # so that we have a probability measure

    X = _greed_init(Y, k)

    a = np.full(k, 1./(2 * k))  # Uniform initialization of weights
    hat_a = np.append(a , 0.5)  # so that we have a probability measure

    for i in range(nb_max_iter):

        C = _build_dist_matrix(X, Y)

        if weight_update:
            hat_a = _weight_optim(hat_a, hat_b, C, t, gamma, nb_max_iter, stopping_criterion)

        new_X, P = _loc_update(hat_a, hat_b, Y, C, gamma)

        ### Compute error update
        diff = np.linalg.norm(new_X - X, axis=1)
        diff = diff[~np.isnan(diff)]
        e = np.linalg.norm(diff)
        if e < stopping_criterion:
            break
        else:
            X = new_X
   
    if weight_update:
        return new_X, P, 2*n * hat_a[:k]

    return new_X, P

