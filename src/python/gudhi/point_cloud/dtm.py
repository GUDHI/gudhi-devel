# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Marc Glisse
#
# Copyright (C) 2020 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from .knn import KNN
import numpy as np


class DTM:
    """
    Class to compute the distance to the empirical measure defined by a point set, as introduced in :cite:`dtm`.
    """

    def __init__(self, k, q=2, **kwargs):
        """
        Args:
            k (int): number of neighbors (possibly including the point itself).
            q (float): order used to compute the distance to measure. Defaults to 2.
            kwargs: same parameters as :class:`~gudhi.point_cloud.knn.KNN`, except that metric="neighbors" means that :func:`transform` expects an array with the distances to the k nearest neighbors.
        """
        self.k = k
        self.q = q
        self.params = kwargs

    def fit_transform(self, X, y=None):
        return self.fit(X).transform(X)

    def fit(self, X, y=None):
        """
        Args:
            X (numpy.array): coordinates for mass points.
        """
        if self.params.setdefault("metric", "euclidean") != "neighbors":
            self.knn = KNN(self.k, return_index=False, return_distance=True, sort_results=False, **self.params)
            self.knn.fit(X)
        return self

    def transform(self, X):
        """
        Args:
            X (numpy.array): coordinates for query points, or distance matrix if metric is "precomputed", or distances to the k nearest neighbors if metric is "neighbors" (if the array has more than k columns, the remaining ones are ignored).
        """
        if self.params["metric"] == "neighbors":
            distances = X[:, : self.k]
        else:
            distances = self.knn.transform(X)
        distances = distances ** self.q
        dtm = distances.sum(-1) / self.k
        dtm = dtm ** (1.0 / self.q)
        # We compute too many powers, 1/p in knn then q in dtm, 1/q in dtm then q or some log in the caller.
        # Add option to skip the final root?
        return dtm


class DTMDensity:
    """
    Density estimator based on the distance to the empirical measure defined by a point set, as defined in :cite:`dtmdensity`. Note that this implementation does not renormalize so the total measure is not 1, see the reference for suitable normalization factors in the Euclidean case.
    """

    def __init__(self, k=None, weights=None, q=None, dim=None, **kwargs):
        """
        Args:
            k (int): number of neighbors (possibly including the point itself).
            weights (numpy.array): weights of each of the k neighbors, optional.
            q (float): order used to compute the distance to measure. Defaults to dim.
            dim (float): final exponent representing the dimension. Defaults to the dimension, and must be specified when the dimension cannot be read from the input (metric="neighbors" or metric="precomputed").
            kwargs: same parameters as :class:`~gudhi.point_cloud.knn.KNN`, except that metric="neighbors" means that :func:`transform` expects an array with the distances to the k nearest neighbors.
        """
        if weights is None:
            assert k is not None, "Must specify k or weights"
            self.k = k
            self.weights = np.full(k, 1.0 / k)
        else:
            self.weights = weights
            self.k = len(weights)
            assert k is None or k == self.k, "k differs from the length of weights"
        self.q = q
        self.dim = dim
        self.params = kwargs

    def fit_transform(self, X, y=None):
        return self.fit(X).transform(X)

    def fit(self, X, y=None):
        """
        Args:
            X (numpy.array): coordinates for mass points.
        """
        if self.params.setdefault("metric", "euclidean") != "neighbors":
            self.knn = KNN(self.k, return_index=False, return_distance=True, sort_results=False, **self.params)
            self.knn.fit(X)
        return self

    def transform(self, X):
        """
        Args:
            X (numpy.array): coordinates for query points, or distance matrix if metric is "precomputed", or distances to the k nearest neighbors if metric is "neighbors" (if the array has more than k columns, the remaining ones are ignored).
        """
        q = self.q
        dim = self.dim
        if dim is None:
            assert self.params["metric"] not in {
                "neighbors",
                "precomputed",
            }, "dim not specified and cannot guess the dimension"
            dim = len(X[0])
        if q is None:
            q = dim
        if self.params["metric"] == "neighbors":
            distances = X[:, : self.k]
        else:
            distances = self.knn.transform(X)
        distances = distances ** q
        dtm = (distances * weights).sum(-1)
        return dtm ** (-dim / q)
        # We compute too many powers, 1/p in knn then q in dtm, d/q in dtm then whatever in the caller.
        # Add option to skip the final root?
