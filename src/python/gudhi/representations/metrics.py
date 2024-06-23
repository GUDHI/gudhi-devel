# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Mathieu Carrière
#
# Copyright (C) 2018-2019 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.metrics import pairwise_distances
from gudhi.hera import wasserstein_distance as hera_wasserstein_distance
from .preprocessing import Padding
from joblib import Parallel, delayed

#############################################
# Metrics ###################################
#############################################

def _sliced_wasserstein_distance(D1, D2, num_directions):
    """
    This is a function for computing the sliced Wasserstein distance from two persistence diagrams. The Sliced Wasserstein distance is computed by projecting the persistence diagrams onto lines, comparing the projections with the 1-norm, and finally averaging over the lines. See http://proceedings.mlr.press/v70/carriere17a.html for more details.
    
    Parameters:
        D1: (n x 2) numpy.array encoding the (finite points of the) first diagram. Must not contain essential points (i.e. with infinite coordinate).
        D2: (m x 2) numpy.array encoding the second diagram.
        num_directions (int): number of lines evenly sampled from [-pi/2,pi/2] in order to approximate and speed up the distance computation.

    Returns: 
        float: the sliced Wasserstein distance between persistence diagrams. 
    """
    thetas = np.linspace(-np.pi/2, np.pi/2, num=num_directions, endpoint=False)[np.newaxis,:-1]
    lines = np.concatenate([np.cos(thetas), np.sin(thetas)], axis=0)
    approx1 = np.matmul(D1, lines)
    approx_diag1 = np.matmul(np.broadcast_to(D1.sum(-1,keepdims=True)/2,(len(D1),2)), lines)
    approx2 = np.matmul(D2, lines)
    approx_diag2 = np.matmul(np.broadcast_to(D2.sum(-1,keepdims=True)/2,(len(D2),2)), lines)
    A = np.sort(np.concatenate([approx1, approx_diag2], axis=0), axis=0)
    B = np.sort(np.concatenate([approx2, approx_diag1], axis=0), axis=0)
    L1 = np.sum(np.abs(A-B), axis=0)
    return np.mean(L1)

def _compute_persistence_diagram_projections(X, num_directions):
    """
    This is a function for projecting the points of a list of persistence diagrams (as well as their diagonal projections) onto a fixed number of lines sampled uniformly on [-pi/2, pi/2]. This function can be used as a preprocessing step in order to speed up the running time for computing all pairwise sliced Wasserstein distances / kernel values on a list of persistence diagrams. 

    Parameters:
        X (list of n numpy arrays of shape (num,2)): list of persistence diagrams.
        num_directions (int): number of lines evenly sampled from [-pi/2,pi/2] in order to approximate and speed up the distance computation.

    Returns: 
        list of n numpy arrays of shape (2*num,num_directions): list of projected persistence diagrams.
    """
    thetas = np.linspace(-np.pi/2, np.pi/2, num=num_directions, endpoint=False)[np.newaxis,:-1]
    lines = np.concatenate([np.cos(thetas), np.sin(thetas)], axis=0)
    XX = [np.vstack([np.matmul(D, lines), np.matmul(np.matmul(D, .5 * np.ones((2,2))), lines)]) for D in X]
    return XX

def _sliced_wasserstein_distance_on_projections(D1, D2):
    """
    This is a function for computing the sliced Wasserstein distance between two persistence diagrams that have already been projected onto some lines. It simply amounts to comparing the sorted projections with the 1-norm, and averaging over the lines. See http://proceedings.mlr.press/v70/carriere17a.html for more details.

    Parameters: 
        D1: (2n x number_of_lines) numpy.array containing the n projected points of the first diagram, and the n projections of their diagonal projections.
        D2: (2m x number_of_lines) numpy.array containing the m projected points of the second diagram, and the m projections of their diagonal projections.

    Returns: 
        float: the sliced Wasserstein distance between the projected persistence diagrams. 
    """
    lim1, lim2 = int(len(D1)/2), int(len(D2)/2)
    approx1, approx_diag1, approx2, approx_diag2 = D1[:lim1], D1[lim1:], D2[:lim2], D2[lim2:]
    A = np.sort(np.concatenate([approx1, approx_diag2], axis=0), axis=0)
    B = np.sort(np.concatenate([approx2, approx_diag1], axis=0), axis=0)
    L1 = np.sum(np.abs(A-B), axis=0)
    return np.mean(L1)

def _persistence_fisher_distance(D1, D2, kernel_approx=None, bandwidth=1.):
    """
    This is a function for computing the persistence Fisher distance from two persistence diagrams. The persistence Fisher distance is obtained by computing the original Fisher distance between the probability distributions associated to the persistence diagrams given by convolving them with a Gaussian kernel. See http://papers.nips.cc/paper/8205-persistence-fisher-kernel-a-riemannian-manifold-kernel-for-persistence-diagrams for more details.

    Parameters: 
        D1: (n x 2) numpy.array encoding the (finite points of the) first diagram). Must not contain essential points (i.e. with infinite coordinate).
        D2: (m x 2) numpy.array encoding the second diagram.
        bandwidth (float): bandwidth of the Gaussian kernel used to turn persistence diagrams into probability distributions.
        kernel_approx: kernel approximation class used to speed up computation. Common kernel approximations classes can be found in the scikit-learn library (such as RBFSampler for instance).   

    Returns: 
        float: the persistence Fisher distance between persistence diagrams. 
    """
    projection = (1./2) * np.ones((2,2))
    diagonal_projections1 = np.matmul(D1, projection)
    diagonal_projections2 = np.matmul(D2, projection)
    if kernel_approx is not None:
        approx1 = kernel_approx.transform(D1)
        approx_diagonal1 = kernel_approx.transform(diagonal_projections1)
        approx2 = kernel_approx.transform(D2)
        approx_diagonal2 = kernel_approx.transform(diagonal_projections2)
        Z = np.concatenate([approx1, approx_diagonal1, approx2, approx_diagonal2], axis=0)
        U, V = np.sum(np.concatenate([approx1, approx_diagonal2], axis=0), axis=0), np.sum(np.concatenate([approx2, approx_diagonal1], axis=0), axis=0) 
        vectori, vectorj = np.abs(np.matmul(Z, U.T)), np.abs(np.matmul(Z, V.T))
        vectori_sum, vectorj_sum = np.sum(vectori), np.sum(vectorj)
        if vectori_sum != 0:
            vectori = vectori/vectori_sum
        if vectorj_sum != 0:
            vectorj = vectorj/vectorj_sum
        return np.arccos(  min(np.dot(np.sqrt(vectori), np.sqrt(vectorj)), 1.)  )
    else:
        Z = np.concatenate([D1, diagonal_projections1, D2, diagonal_projections2], axis=0)
        U, V = np.concatenate([D1, diagonal_projections2], axis=0), np.concatenate([D2, diagonal_projections1], axis=0) 
        vectori = np.sum(np.exp(-np.square(pairwise_distances(Z,U))/(2 * np.square(bandwidth)))/(bandwidth * np.sqrt(2*np.pi)), axis=1)
        vectorj = np.sum(np.exp(-np.square(pairwise_distances(Z,V))/(2 * np.square(bandwidth)))/(bandwidth * np.sqrt(2*np.pi)), axis=1)
        vectori_sum, vectorj_sum = np.sum(vectori), np.sum(vectorj)
        if vectori_sum != 0:
            vectori = vectori/vectori_sum
        if vectorj_sum != 0:
            vectorj = vectorj/vectorj_sum
        return np.arccos(  min(np.dot(np.sqrt(vectori), np.sqrt(vectorj)), 1.)  )

def _pairwise(fallback, skipdiag, X, Y, metric, n_jobs):
    if Y is not None:
        return fallback(X, Y, metric=metric, n_jobs=n_jobs)
    triu = np.triu_indices(len(X), k=skipdiag)
    tril = (triu[1], triu[0])
    par = Parallel(n_jobs=n_jobs, prefer="threads")
    d = par(delayed(metric)([triu[0][i]], [triu[1][i]]) for i in range(len(triu[0])))
    m = np.empty((len(X), len(X)))
    m[triu] = d
    m[tril] = d
    if skipdiag:
        np.fill_diagonal(m, 0)
    return m

def _sklearn_wrapper(metric, X, Y, **kwargs):
    """
    This function is a wrapper for any metric between two persistence diagrams that takes two numpy arrays of shapes (nx2) and (mx2) as arguments.
    """
    if Y is None:
        def flat_metric(a, b):
            return metric(X[int(a[0])], X[int(b[0])], **kwargs)
    else:
        def flat_metric(a, b):
            return metric(X[int(a[0])], Y[int(b[0])], **kwargs)
    return flat_metric

PAIRWISE_DISTANCE_FUNCTIONS = {
    "wasserstein": hera_wasserstein_distance,
    "hera_wasserstein": hera_wasserstein_distance,
    "persistence_fisher": _persistence_fisher_distance,
}

def pairwise_persistence_diagram_distances(X, Y=None, metric="bottleneck", n_jobs=None, **kwargs):
    """
    This function computes the distance matrix between two lists of persistence diagrams given as numpy arrays of shape (nx2).

    Parameters:
        X (list of n numpy arrays of shape (numx2)): first list of persistence diagrams. 
        Y (list of m numpy arrays of shape (numx2)): second list of persistence diagrams (optional). If None, pairwise distances are computed from the first list only.
        metric: distance to use. It can be either a string ("sliced_wasserstein", "wasserstein", "hera_wasserstein" (Wasserstein distance computed with Hera---note that Hera is also used for the default option "wasserstein"), "pot_wasserstein" (Wasserstein distance computed with POT), "bottleneck", "persistence_fisher") or a function taking two numpy arrays of shape (nx2) and (mx2) as inputs. If it is a function, make sure that it is symmetric and that it outputs 0 if called on the same two arrays. 
        n_jobs (int): number of jobs to use for the computation. This uses joblib.Parallel(prefer="threads"), so metrics that do not release the GIL may not scale unless run inside a `joblib.parallel_backend <https://joblib.readthedocs.io/en/latest/parallel.html#joblib.parallel_backend>`_ block.
        **kwargs: optional keyword parameters. Any further parameters are passed directly to the distance function. See the docs of the various distance classes in this module.

    Returns: 
        numpy array of shape (nxm): distance matrix
    """
    XX = np.reshape(np.arange(len(X)), [-1,1])
    YY = None if Y is None or Y is X else np.reshape(np.arange(len(Y)), [-1,1])
    if metric == "bottleneck":
        try: 
            from .. import bottleneck_distance
            return _pairwise(pairwise_distances, True, XX, YY, metric=_sklearn_wrapper(bottleneck_distance, X, Y, **kwargs), n_jobs=n_jobs)
        except ImportError:
            print("Gudhi built without CGAL")
            raise
    elif metric == "pot_wasserstein":
        try:
            from gudhi.wasserstein import wasserstein_distance as pot_wasserstein_distance
            return _pairwise(pairwise_distances, True, XX, YY, metric=_sklearn_wrapper(pot_wasserstein_distance,  X, Y, **kwargs), n_jobs=n_jobs)
        except ImportError:
            print("POT (Python Optimal Transport) is not installed. Please install POT or use metric='wasserstein' or metric='hera_wasserstein'")
            raise
    elif metric == "sliced_wasserstein":
        Xproj = _compute_persistence_diagram_projections(X, **kwargs)
        Yproj = None if Y is None else _compute_persistence_diagram_projections(Y, **kwargs)
        return _pairwise(pairwise_distances, True, XX, YY, metric=_sklearn_wrapper(_sliced_wasserstein_distance_on_projections, Xproj, Yproj), n_jobs=n_jobs)
    elif type(metric) == str:
        return _pairwise(pairwise_distances, True, XX, YY, metric=_sklearn_wrapper(PAIRWISE_DISTANCE_FUNCTIONS[metric], X, Y, **kwargs), n_jobs=n_jobs)
    else:
        return _pairwise(pairwise_distances, True, XX, YY, metric=_sklearn_wrapper(metric, X, Y, **kwargs), n_jobs=n_jobs)

class SlicedWassersteinDistance(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the sliced Wasserstein distance matrix from a list of persistence diagrams. The Sliced Wasserstein distance is computed by projecting the persistence diagrams onto lines, comparing the projections with the 1-norm, and finally integrating over all possible lines. See http://proceedings.mlr.press/v70/carriere17a.html for more details. 
    """
    def __init__(self, num_directions=10, n_jobs=None):
        """
        Constructor for the SlicedWassersteinDistance class.

        Parameters:
            num_directions (int): number of lines evenly sampled from [-pi/2,pi/2] in order to approximate and speed up the distance computation (default 10). 
            n_jobs (int): number of jobs to use for the computation. See :func:`pairwise_persistence_diagram_distances` for details.
        """
        self.num_directions = num_directions
        self.n_jobs = n_jobs

    def fit(self, X, y=None):
        """
        Fit the SlicedWassersteinDistance class on a list of persistence diagrams: persistence diagrams are projected onto the different lines. The diagrams themselves and their projections are then stored in numpy arrays, called **diagrams_** and **approx_diag_**.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        self.diagrams_ = X
        return self

    def transform(self, X):
        """
        Compute all sliced Wasserstein distances between the persistence diagrams that were stored after calling the fit() method, and a given list of (possibly different) persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.

        Returns:
            numpy array of shape (number of diagrams in **diagrams**) x (number of diagrams in X): matrix of pairwise sliced Wasserstein distances.
        """
        return pairwise_persistence_diagram_distances(X, self.diagrams_, metric="sliced_wasserstein", num_directions=self.num_directions, n_jobs=self.n_jobs)

    def __call__(self, diag1, diag2):
        """
        Apply SlicedWassersteinDistance on a single pair of persistence diagrams and outputs the result.

        Parameters:
            diag1 (n x 2 numpy array): first input persistence diagram.
            diag2 (n x 2 numpy array): second input persistence diagram.

        Returns:
            float: sliced Wasserstein distance.
        """
        return _sliced_wasserstein_distance(diag1, diag2, num_directions=self.num_directions)

class BottleneckDistance(BaseEstimator, TransformerMixin):
    r"""
    This is a class for computing the bottleneck distance matrix from a list of persistence diagrams.

    :Requires: `CGAL <installation.html#cgal>`_
    """
    def __init__(self, epsilon=None, n_jobs=None):
        """
        Constructor for the BottleneckDistance class.

        Parameters:
            epsilon (double): absolute (additive) error tolerated on the distance (default is the smallest positive float), see :func:`gudhi.bottleneck_distance`.
            n_jobs (int): number of jobs to use for the computation. See :func:`pairwise_persistence_diagram_distances` for details.
        """
        self.epsilon = epsilon
        self.n_jobs = n_jobs

    def fit(self, X, y=None):
        """
        Fit the BottleneckDistance class on a list of persistence diagrams: persistence diagrams are stored in a numpy array called **diagrams**.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        self.diagrams_ = X
        return self

    def transform(self, X):
        """
        Compute all bottleneck distances between the persistence diagrams that were stored after calling the fit() method, and a given list of (possibly different) persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.

        Returns:
            numpy array of shape (number of diagrams in **diagrams**) x (number of diagrams in X): matrix of pairwise bottleneck distances.
        """
        Xfit = pairwise_persistence_diagram_distances(X, self.diagrams_, metric="bottleneck", e=self.epsilon, n_jobs=self.n_jobs)
        return Xfit

    def __call__(self, diag1, diag2):
        """
        Apply BottleneckDistance on a single pair of persistence diagrams and outputs the result.

        Parameters:
            diag1 (n x 2 numpy array): first input persistence diagram.
            diag2 (n x 2 numpy array): second input persistence diagram.

        Returns:
            float: bottleneck distance.
        """
        try: 
            from .. import bottleneck_distance
            return bottleneck_distance(diag1, diag2, e=self.epsilon)
        except ImportError:
            print("Gudhi built without CGAL")
            raise

class PersistenceFisherDistance(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the persistence Fisher distance matrix from a list of persistence diagrams. The persistence Fisher distance is obtained by computing the original Fisher distance between the probability distributions associated to the persistence diagrams given by convolving them with a Gaussian kernel. See http://papers.nips.cc/paper/8205-persistence-fisher-kernel-a-riemannian-manifold-kernel-for-persistence-diagrams for more details. 
    """
    def __init__(self, bandwidth=1., kernel_approx=None, n_jobs=None):
        """
        Constructor for the PersistenceFisherDistance class.

        Parameters:
            bandwidth (double): bandwidth of the Gaussian kernel used to turn persistence diagrams into probability distributions (default 1.).
            kernel_approx (class): kernel approximation class used to speed up computation (default None). Common kernel approximations classes can be found in the scikit-learn library (such as RBFSampler for instance).   
            n_jobs (int): number of jobs to use for the computation. See :func:`pairwise_persistence_diagram_distances` for details.
        """
        self.bandwidth, self.kernel_approx = bandwidth, kernel_approx
        self.n_jobs = n_jobs

    def fit(self, X, y=None):
        """
        Fit the PersistenceFisherDistance class on a list of persistence diagrams: persistence diagrams are stored in a numpy array called **diagrams** and the kernel approximation class (if not None) is applied on them.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        self.diagrams_ = X
        return self

    def transform(self, X):
        """
        Compute all persistence Fisher distances between the persistence diagrams that were stored after calling the fit() method, and a given list of (possibly different) persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.

        Returns:
            numpy array of shape (number of diagrams in **diagrams**) x (number of diagrams in X): matrix of pairwise persistence Fisher distances.
        """
        return pairwise_persistence_diagram_distances(X, self.diagrams_, metric="persistence_fisher", bandwidth=self.bandwidth, kernel_approx=self.kernel_approx, n_jobs=self.n_jobs)

    def __call__(self, diag1, diag2):
        """
        Apply PersistenceFisherDistance on a single pair of persistence diagrams and outputs the result.

        Parameters:
            diag1 (n x 2 numpy array): first input persistence diagram.
            diag2 (n x 2 numpy array): second input persistence diagram.

        Returns:
            float: persistence Fisher distance.
        """
        return _persistence_fisher_distance(diag1, diag2, bandwidth=self.bandwidth, kernel_approx=self.kernel_approx)


class WassersteinDistance(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the Wasserstein distance matrix from a list of persistence diagrams. 
    """

    def __init__(self, order=1, internal_p=np.inf, mode="hera", delta=0.01, n_jobs=None):
        """
        Constructor for the WassersteinDistance class.

        Parameters:
            order (int): exponent for Wasserstein, default value is 1., see :func:`gudhi.wasserstein.wasserstein_distance`.
            internal_p (int): ground metric on the (upper-half) plane (i.e. norm l_p in R^2), default value is `np.inf`, see :func:`gudhi.wasserstein.wasserstein_distance`.
            mode (str): method for computing Wasserstein distance. Either "pot" or "hera". Default set to "hera".
            delta (float): relative error 1+delta. Used only if mode == "hera".
            n_jobs (int): number of jobs to use for the computation. See :func:`pairwise_persistence_diagram_distances` for details.
        """
        self.order, self.internal_p, self.mode = order, internal_p, mode
        self.delta = delta
        self.n_jobs = n_jobs

    def fit(self, X, y=None):
        """
        Fit the WassersteinDistance class on a list of persistence diagrams: persistence diagrams are stored in a numpy array called **diagrams**.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        if self.mode not in ("pot", "hera"):
            raise NameError("Unknown mode. Current available values for mode are 'hera' and 'pot'")
        self.diagrams_ = X
        return self

    def transform(self, X):
        """
        Compute all Wasserstein distances between the persistence diagrams that were stored after calling the fit() method, and a given list of (possibly different) persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.

        Returns:
            numpy array of shape (number of diagrams in **diagrams**) x (number of diagrams in X): matrix of pairwise Wasserstein distances.
        """
        if self.mode == "hera":
            Xfit = pairwise_persistence_diagram_distances(X, self.diagrams_, metric="hera_wasserstein", order=self.order, internal_p=self.internal_p, delta=self.delta, n_jobs=self.n_jobs)
        else:
            Xfit = pairwise_persistence_diagram_distances(X, self.diagrams_, metric="pot_wasserstein", order=self.order, internal_p=self.internal_p, matching=False, n_jobs=self.n_jobs)
        return Xfit

    def __call__(self, diag1, diag2):
        """
        Apply WassersteinDistance on a single pair of persistence diagrams and outputs the result.

        Parameters:
            diag1 (n x 2 numpy array): first input persistence diagram.
            diag2 (n x 2 numpy array): second input persistence diagram.

        Returns:
            float: Wasserstein distance.
        """
        if self.mode == "hera":
            return hera_wasserstein_distance(diag1, diag2, order=self.order, internal_p=self.internal_p, delta=self.delta)
        elif self.mode == "pot":
            try:
                from gudhi.wasserstein import wasserstein_distance as pot_wasserstein_distance
                return pot_wasserstein_distance(diag1, diag2, order=self.order, internal_p=self.internal_p, matching=False)
            except ImportError:
                print("POT (Python Optimal Transport) is not installed. Please install POT or use mode='hera'")
                raise
        else:
            raise NameError("Unknown mode. Current available values for mode are 'hera' and 'pot'")
