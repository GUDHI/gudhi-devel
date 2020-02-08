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
from gudhi.wasserstein import wasserstein_distance
from .preprocessing import Padding

try:
    from .. import bottleneck_distance
    USE_GUDHI = True
except ImportError:
    USE_GUDHI = False
    print("Gudhi built without CGAL: BottleneckDistance will return a null matrix")

#############################################
# Metrics ###################################
#############################################

def sliced_wasserstein_distance(D1, D2, num_directions):
    """
    This is a function for computing the sliced Wasserstein distance from two persistence diagrams. The Sliced Wasserstein distance is computed by projecting the persistence diagrams onto lines, comparing the projections with the 1-norm, and finally integrating over all possible lines. See http://proceedings.mlr.press/v70/carriere17a.html for more details.
    :param D1: (n x 2) numpy.array encoding the (finite points of the) first diagram. Must not contain essential points (i.e. with infinite coordinate).
    :param D2: (m x 2) numpy.array encoding the second diagram.
    :param num_directions: number of lines evenly sampled from [-pi/2,pi/2] in order to approximate and speed up the distance computation.
    :returns: the sliced Wasserstein distance between persistence diagrams.
    :rtype: float 
    """
    thetas = np.linspace(-np.pi/2, np.pi/2, num=num_directions+1)[np.newaxis,:-1]
    lines = np.concatenate([np.cos(thetas), np.sin(thetas)], axis=0)
    approx1 = np.matmul(D1, lines)
    diag_proj1 = (1./2) * np.ones((2,2))
    approx_diag1 = np.matmul(np.matmul(D1, diag_proj1), lines)
    approx2 = np.matmul(D2, lines)
    diag_proj2 = (1./2) * np.ones((2,2))
    approx_diag2 = np.matmul(np.matmul(D2, diag_proj2), lines)
    A = np.sort(np.concatenate([approx1, approx_diag2], axis=0), axis=0)
    B = np.sort(np.concatenate([approx2, approx_diag1], axis=0), axis=0)
    L1 = np.sum(np.abs(A-B), axis=0)
    return np.mean(L1)

def persistence_fisher_distance(D1, D2, kernel_approx=None, bandwidth=1.):
    """
    This is a function for computing the persistence Fisher distance from two persistence diagrams. The persistence Fisher distance is obtained by computing the original Fisher distance between the probability distributions associated to the persistence diagrams given by convolving them with a Gaussian kernel. See http://papers.nips.cc/paper/8205-persistence-fisher-kernel-a-riemannian-manifold-kernel-for-persistence-diagrams for more details.
    :param D1: (n x 2) numpy.array encoding the (finite points of the) first diagram. Must not contain essential points (i.e. with infinite coordinate).
    :param D2: (m x 2) numpy.array encoding the second diagram.
    :param bandwidth: bandwidth of the Gaussian kernel used to turn persistence diagrams into probability distributions.
    :param kernel_approx: kernel approximation class used to speed up computation. Common kernel approximations classes can be found in the scikit-learn library (such as RBFSampler for instance).   
    :returns: the persistence Fisher distance between persistence diagrams.
    :rtype: float 
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

def sklearn_wrapper(metric, X, Y, **kwargs):
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

def pairwise_persistence_diagram_distances(X, Y=None, metric="bottleneck", **kwargs):
    """
    This function computes the distance matrix between two lists of persistence diagrams given as numpy arrays of shape (nx2).
    :param X: first list of persistence diagrams. 
    :param Y: second list of persistence diagrams (optional). If None, pairwise distances are computed from the first list only.
    :param metric: distance to use. It can be either a string ("sliced_wasserstein", "wasserstein", "bottleneck", "persistence_fisher") or a function taking two numpy arrays of shape (nx2) and (mx2) as inputs.
    :returns: distance matrix, i.e., numpy array of shape (num diagrams 1 x num diagrams 2)
    :rtype: float
    """
    XX = np.reshape(np.arange(len(X)), [-1,1])
    YY = None if Y is None else np.reshape(np.arange(len(Y)), [-1,1]) 
    if metric == "bottleneck":
        return pairwise_distances(XX, YY, metric=sklearn_wrapper(bottleneck_distance, X, Y, **kwargs))
    elif metric == "wasserstein":
        return pairwise_distances(XX, YY, metric=sklearn_wrapper(wasserstein_distance, X, Y, **kwargs))
    elif metric == "sliced_wasserstein":
        return pairwise_distances(XX, YY, metric=sklearn_wrapper(sliced_wasserstein_distance, X, Y, **kwargs))
    elif metric == "persistence_fisher":
        return pairwise_distances(XX, YY, metric=sklearn_wrapper(persistence_fisher_distance, X, Y, **kwargs))
    else:
        return pairwise_distances(XX, YY, metric=sklearn_wrapper(metric, X, Y, **kwargs))

class SlicedWassersteinDistance(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the sliced Wasserstein distance matrix from a list of persistence diagrams. The Sliced Wasserstein distance is computed by projecting the persistence diagrams onto lines, comparing the projections with the 1-norm, and finally integrating over all possible lines. See http://proceedings.mlr.press/v70/carriere17a.html for more details. 
    """
    def __init__(self, num_directions=10):
        """
        Constructor for the SlicedWassersteinDistance class.

        Parameters:
            num_directions (int): number of lines evenly sampled from [-pi/2,pi/2] in order to approximate and speed up the distance computation (default 10). 
        """
        self.num_directions = num_directions

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
        return pairwise_persistence_diagram_distances(X, self.diagrams_, metric="sliced_wasserstein", num_directions=self.num_directions)

class BottleneckDistance(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the bottleneck distance matrix from a list of persistence diagrams. 
    """
    def __init__(self, epsilon=None):
        """
        Constructor for the BottleneckDistance class.

        Parameters:
            epsilon (double): absolute (additive) error tolerated on the distance (default is the smallest positive float), see :func:`gudhi.bottleneck_distance`.
        """
        self.epsilon = epsilon

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
        if not USE_GUDHI:
            print("Gudhi built without CGAL: returning a null matrix")
        Xfit = pairwise_persistence_diagram_distances(X, self.diagrams_, metric="bottleneck", e=self.epsilon) if USE_GUDHI else np.zeros((len(X), len(self.diagrams_)))
        return Xfit

class WassersteinDistance(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the Wasserstein distance matrix from a list of persistence diagrams. 
    """
    def __init__(self, order=2, internal_p=2):
        """
        Constructor for the WassersteinDistance class.

        Parameters:
            order (int): exponent for Wasserstein, default value is 2., see :func:`gudhi.wasserstein.wasserstein_distance`.
            internal_p (int): ground metric on the (upper-half) plane (i.e. norm l_p in R^2), default value is 2 (euclidean norm), see :func:`gudhi.wasserstein.wasserstein_distance`.
        """
        self.order, self.internal_p = order, internal_p

    def fit(self, X, y=None):
        """
        Fit the WassersteinDistance class on a list of persistence diagrams: persistence diagrams are stored in a numpy array called **diagrams**.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
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
        return pairwise_persistence_diagram_distances(X, self.diagrams_, metric="wasserstein", order=self.order, internal_p=self.internal_p)

class PersistenceFisherDistance(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the persistence Fisher distance matrix from a list of persistence diagrams. The persistence Fisher distance is obtained by computing the original Fisher distance between the probability distributions associated to the persistence diagrams given by convolving them with a Gaussian kernel. See http://papers.nips.cc/paper/8205-persistence-fisher-kernel-a-riemannian-manifold-kernel-for-persistence-diagrams for more details. 
    """
    def __init__(self, bandwidth=1., kernel_approx=None):
        """
        Constructor for the PersistenceFisherDistance class.

        Parameters:
            bandwidth (double): bandwidth of the Gaussian kernel used to turn persistence diagrams into probability distributions (default 1.).
            kernel_approx (class): kernel approximation class used to speed up computation (default None). Common kernel approximations classes can be found in the scikit-learn library (such as RBFSampler for instance).   
        """
        self.bandwidth, self.kernel_approx = bandwidth, kernel_approx

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
        return pairwise_persistence_diagram_distances(X, self.diagrams_, metric="persistence_fisher", bandwidth=self.bandwidth, kernel_approx=self.kernel_approx)
