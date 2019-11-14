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
try:
    from .. import bottleneck_distance
    USE_GUDHI = True
except ImportError:
    USE_GUDHI = False
    print("Gudhi built without CGAL: BottleneckDistance will return a null matrix")

#############################################
# Metrics ###################################
#############################################

class SlicedWassersteinDistance(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the sliced Wasserstein distance matrix from a list of persistence diagrams. The Sliced Wasserstein distance is computed by projecting the persistence diagrams onto lines, comparing the projections with the 1-norm, and finally integrating over all possible lines. See http://proceedings.mlr.press/v70/carriere17a.html for more details. 
    """
    def __init__(self, num_directions=10):
        """
        Constructor for the SlicedWassersteinDistance class.

        Attributes:
            num_directions (int): number of lines evenly sampled from [-pi/2,pi/2] in order to approximate and speed up the distance computation (default 10). 
        """
        self.num_directions = num_directions
        thetas = np.linspace(-np.pi/2, np.pi/2, num=self.num_directions+1)[np.newaxis,:-1]
        self.lines_ = np.concatenate([np.cos(thetas), np.sin(thetas)], axis=0)

    def fit(self, X, y=None):
        """
        Fit the SlicedWassersteinDistance class on a list of persistence diagrams: persistence diagrams are projected onto the different lines. The diagrams themselves and their projections are then stored in numpy arrays, called **diagrams_** and **approx_diag_**.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.
            y (n x 1 array): persistence diagram labels (unused).
        """
        self.diagrams_ = X
        self.approx_ = [np.matmul(X[i], self.lines_) for i in range(len(X))]
        diag_proj = (1./2) * np.ones((2,2))
        self.approx_diag_ = [np.matmul(np.matmul(X[i], diag_proj), self.lines_) for i in range(len(X))]
        return self

    def transform(self, X):
        """
        Compute all sliced Wasserstein distances between the persistence diagrams that were stored after calling the fit() method, and a given list of (possibly different) persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.

        Returns:
            Xfit (numpy array of shape (number of diagrams in **diagrams**) x (number of diagrams in X)): matrix of pairwise sliced Wasserstein distances.
        """
        Xfit = np.zeros((len(X), len(self.approx_)))
        if len(self.diagrams_) == len(X) and np.all([np.array_equal(self.diagrams_[i], X[i]) for i in range(len(X))]):
            for i in range(len(self.approx_)):
                for j in range(i+1, len(self.approx_)):
                    A = np.sort(np.concatenate([self.approx_[i], self.approx_diag_[j]], axis=0), axis=0)
                    B = np.sort(np.concatenate([self.approx_[j], self.approx_diag_[i]], axis=0), axis=0)
                    L1 = np.sum(np.abs(A-B), axis=0)
                    Xfit[i,j] = np.mean(L1)
                    Xfit[j,i] = Xfit[i,j]
        else:
            diag_proj = (1./2) * np.ones((2,2))
            approx = [np.matmul(X[i], self.lines_) for i in range(len(X))]
            approx_diag = [np.matmul(np.matmul(X[i], diag_proj), self.lines_) for i in range(len(X))]
            for i in range(len(approx)):
                for j in range(len(self.approx_)):
                    A = np.sort(np.concatenate([approx[i], self.approx_diag_[j]], axis=0), axis=0)
                    B = np.sort(np.concatenate([self.approx_[j], approx_diag[i]], axis=0), axis=0)
                    L1 = np.sum(np.abs(A-B), axis=0)
                    Xfit[i,j] = np.mean(L1)

        return Xfit

class BottleneckDistance(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the bottleneck distance matrix from a list of persistence diagrams. 
    """
    def __init__(self, epsilon=1e-3):
        """
        Constructor for the BottleneckDistance class.

        Attributes:
            epsilon (double): approximation quality (default 1e-4).
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
            Xfit (numpy array of shape (number of diagrams in **diagrams**) x (number of diagrams in X)): matrix of pairwise bottleneck distances.
        """
        num_diag1 = len(X)

        if len(self.diagrams_) == len(X) and np.all([np.array_equal(self.diagrams_[i], X[i]) for i in range(len(X))]):
            matrix = np.zeros((num_diag1, num_diag1))

            if USE_GUDHI:
                for i in range(num_diag1):
                    for j in range(i+1, num_diag1):
                        matrix[i,j] = bottleneck_distance(X[i], X[j], self.epsilon)
                        matrix[j,i] = matrix[i,j]
            else:
                print("Gudhi required---returning null matrix")

        else:
            num_diag2 = len(self.diagrams_)
            matrix = np.zeros((num_diag1, num_diag2))

            if USE_GUDHI:
                for i in range(num_diag1):
                    for j in range(num_diag2):
                        matrix[i,j] = bottleneck_distance(X[i], self.diagrams_[j], self.epsilon)
            else:
                print("Gudhi required---returning null matrix")

        Xfit = matrix

        return Xfit

class PersistenceFisherDistance(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the persistence Fisher distance matrix from a list of persistence diagrams. The persistence Fisher distance is obtained by computing the original Fisher distance between the probability distributions associated to the persistence diagrams given by convolving them with a Gaussian kernel. See http://papers.nips.cc/paper/8205-persistence-fisher-kernel-a-riemannian-manifold-kernel-for-persistence-diagrams for more details. 
    """
    def __init__(self, bandwidth=1., kernel_approx=None):
        """
        Constructor for the PersistenceFisherDistance class.

        Attributes:
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
        projection = (1./2) * np.ones((2,2))
        self.diagonal_projections_ = [np.matmul(X[i], projection) for i in range(len(X))]
        if self.kernel_approx is not None:
            self.approx_ = [self.kernel_approx.transform(X[i]) for i in range(len(X))]
            self.approx_diagonal_ = [self.kernel_approx.transform(self.diagonal_projections_[i]) for i in range(len(X))]
        return self

    def transform(self, X):
        """
        Compute all persistence Fisher distances between the persistence diagrams that were stored after calling the fit() method, and a given list of (possibly different) persistence diagrams.

        Parameters:
            X (list of n x 2 numpy arrays): input persistence diagrams.

        Returns:
            Xfit (numpy array of shape (number of diagrams in **diagrams**) x (number of diagrams in X)): matrix of pairwise persistence Fisher distances.
        """
        Xfit = np.zeros((len(X), len(self.diagrams_)))
        if len(self.diagrams_) == len(X) and np.all([np.array_equal(self.diagrams_[i], X[i]) for i in range(len(X))]):
            for i in range(len(self.diagrams_)):
                for j in range(i+1, len(self.diagrams_)):
                    if self.kernel_approx is not None:
                        Z = np.concatenate([self.approx_[i], self.approx_diagonal_[i], self.approx_[j], self.approx_diagonal_[j]], axis=0)
                        U, V = np.sum(np.concatenate([self.approx_[i], self.approx_diagonal_[j]], axis=0), axis=0), np.sum(np.concatenate([self.approx_[j], self.approx_diagonal_[i]], axis=0), axis=0) 
                        vectori, vectorj = np.abs(np.matmul(Z, U.T)), np.abs(np.matmul(Z, V.T))
                        vectori_sum, vectorj_sum = np.sum(vectori), np.sum(vectorj)
                        if vectori_sum != 0:
                            vectori = vectori/vectori_sum
                        if vectorj_sum != 0:
                            vectorj = vectorj/vectorj_sum
                        Xfit[i,j] = np.arccos(  min(np.dot(np.sqrt(vectori), np.sqrt(vectorj)), 1.)  )
                        Xfit[j,i] = Xfit[i,j]
                    else:
                        Z = np.concatenate([self.diagrams_[i], self.diagonal_projections_[i], self.diagrams_[j], self.diagonal_projections_[j]], axis=0)
                        U, V = np.concatenate([self.diagrams_[i], self.diagonal_projections_[j]], axis=0), np.concatenate([self.diagrams_[j], self.diagonal_projections_[i]], axis=0) 
                        vectori = np.sum(np.exp(-np.square(pairwise_distances(Z,U))/(2 * np.square(self.bandwidth)))/(self.bandwidth * np.sqrt(2*np.pi)), axis=1)
                        vectorj = np.sum(np.exp(-np.square(pairwise_distances(Z,V))/(2 * np.square(self.bandwidth)))/(self.bandwidth * np.sqrt(2*np.pi)), axis=1)
                        vectori_sum, vectorj_sum = np.sum(vectori), np.sum(vectorj)
                        if vectori_sum != 0:
                            vectori = vectori/vectori_sum
                        if vectorj_sum != 0:
                            vectorj = vectorj/vectorj_sum
                        Xfit[i,j] = np.arccos(  min(np.dot(np.sqrt(vectori), np.sqrt(vectorj)), 1.)  )
                        Xfit[j,i] = Xfit[i,j]
        else:
            projection = (1./2) * np.ones((2,2))
            diagonal_projections = [np.matmul(X[i], projection) for i in range(len(X))]
            if self.kernel_approx is not None:
                approx = [self.kernel_approx.transform(X[i]) for i in range(len(X))]
                approx_diagonal = [self.kernel_approx.transform(diagonal_projections[i]) for i in range(len(X))]
            for i in range(len(X)):
                for j in range(len(self.diagrams_)):
                    if self.kernel_approx is not None:
                        Z = np.concatenate([approx[i], approx_diagonal[i], self.approx_[j], self.approx_diagonal_[j]], axis=0)
                        U, V = np.sum(np.concatenate([approx[i], self.approx_diagonal_[j]], axis=0), axis=0), np.sum(np.concatenate([self.approx_[j], approx_diagonal[i]], axis=0), axis=0) 
                        vectori, vectorj = np.abs(np.matmul(Z, U.T)), np.abs(np.matmul(Z, V.T))
                        vectori_sum, vectorj_sum = np.sum(vectori), np.sum(vectorj)
                        if vectori_sum != 0:
                            vectori = vectori/vectori_sum
                        if vectorj_sum != 0:
                            vectorj = vectorj/vectorj_sum
                        Xfit[i,j] = np.arccos(  min(np.dot(np.sqrt(vectori), np.sqrt(vectorj)), 1.)  )
                    else:
                        Z = np.concatenate([X[i], diagonal_projections[i], self.diagrams_[j], self.diagonal_projections_[j]], axis=0)
                        U, V = np.concatenate([X[i], self.diagonal_projections_[j]], axis=0), np.concatenate([self.diagrams_[j], diagonal_projections[i]], axis=0) 
                        vectori = np.sum(np.exp(-np.square(pairwise_distances(Z,U))/(2 * np.square(self.bandwidth)))/(self.bandwidth * np.sqrt(2*np.pi)), axis=1)
                        vectorj = np.sum(np.exp(-np.square(pairwise_distances(Z,V))/(2 * np.square(self.bandwidth)))/(self.bandwidth * np.sqrt(2*np.pi)), axis=1)
                        vectori_sum, vectorj_sum = np.sum(vectori), np.sum(vectorj)
                        if vectori_sum != 0:
                            vectori = vectori/vectori_sum
                        if vectorj_sum != 0:
                            vectorj = vectorj/vectorj_sum
                        Xfit[i,j] = np.arccos(  min(np.dot(np.sqrt(vectori), np.sqrt(vectorj)), 1.)  )
        return Xfit
