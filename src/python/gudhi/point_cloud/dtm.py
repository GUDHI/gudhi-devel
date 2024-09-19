# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Marc Glisse
#
# Copyright (C) 2020 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from .knn import KNearestNeighbors
import numpy as np

__author__ = "Marc Glisse"
__copyright__ = "Copyright (C) 2020 Inria"
__license__ = "MIT"


class DistanceToMeasure:
    """
    Class to compute the distance to the empirical measure defined by a point set, 
    as introduced in :cite:`dtmgeoinference2011`.
    """

    def __init__(self, k, q=2, **kwargs):
        """
        Args:
            k (int): number of neighbors (possibly including the point itself).
            q (float): order used to compute the distance to measure. Defaults to 2.
            kwargs: same parameters as :class:`~gudhi.point_cloud.knn.KNearestNeighbors`, except that
                metric="neighbors" means that :func:`transform` expects an array with the distances
                to the k nearest neighbors.
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
            self.knn = KNearestNeighbors(
                self.k, return_index=False, return_distance=True, sort_results=False, **self.params
            )
            self.knn.fit(X)
        return self

    def transform(self, X):
        """
        Args:
            X (numpy.array): coordinates for query points, or distance matrix if metric is "precomputed",
                or distances to the k nearest neighbors if metric is "neighbors" (if the array has more
                than k columns, the remaining ones are ignored).

        Returns:
            numpy.array: a 1-d array with, for each point of X, its distance to the measure defined
            by the argument of :func:`fit`.
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
    Density estimator based on the distance to the empirical measure defined by a point set, as defined
    in :cite:`dtmdensity`. Note that this implementation only renormalizes when asked, and the renormalization
    only works for a Euclidean metric, so in other cases the total measure may not be 1.

    .. note:: When the dimension is high, using it as an exponent can quickly lead to under- or overflows.
        We recommend using a small fixed value instead (for both dim and q) in those cases, even if it won't
        have the same nice theoretical properties as the dimension.
    """

    def __init__(self, k=None, weights=None, q=None, dim=None, normalize=False, n_samples=None, **kwargs):
        """
        Args:
            k (int): number of neighbors (possibly including the point itself). Optional if it can be guessed
                from weights or metric="neighbors".
            weights (numpy.array): weights of each of the k neighbors, optional. They are supposed to sum to 1.
            q (float): order used to compute the distance to measure. Defaults to dim.
            dim (float): final exponent representing the dimension. Defaults to the dimension, and must be specified
                when the dimension cannot be read from the input (metric is "neighbors" or "precomputed").
            normalize (bool): normalize the density so it corresponds to a probability measure on ℝᵈ.
                Only available for the Euclidean metric, defaults to False.
            n_samples (int): number of sample points used for fitting. Only needed if `normalize` is True and
                metric is "neighbors".
            kwargs: same parameters as :class:`~gudhi.point_cloud.knn.KNearestNeighbors`, except that
                metric="neighbors" means that :func:`transform` expects an array with the distances to
                the k nearest neighbors.
        """
        if weights is None:
            self.k = k
            if k is None:
                assert kwargs.get("metric") == "neighbors", 'Must specify k or weights, unless metric is "neighbors"'
                self.weights = None
            else:
                self.weights = np.full(k, 1.0 / k)
        else:
            self.weights = weights
            self.k = len(weights)
            assert k is None or k == self.k, "k differs from the length of weights"
        self.q = q
        self.dim = dim
        self.params = kwargs
        self.normalize = normalize
        self.n_samples = n_samples

    def fit_transform(self, X, y=None):
        return self.fit(X).transform(X)

    def fit(self, X, y=None):
        """
        Args:
            X (numpy.array): coordinates for mass points.
        """
        if self.params.setdefault("metric", "euclidean") != "neighbors":
            self.knn = KNearestNeighbors(
                self.k, return_index=False, return_distance=True, sort_results=False, **self.params
            )
            self.knn.fit(X)
            if self.params["metric"] != "precomputed":
                self.n_samples = len(X)
        return self

    def transform(self, X):
        """
        Args:
            X (numpy.array): coordinates for query points, or distance matrix if metric is "precomputed",
                or distances to the k nearest neighbors if metric is "neighbors" (if the array has more
                than k columns, the remaining ones are ignored).
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
        k = self.k
        weights = self.weights
        if self.params["metric"] == "neighbors":
            distances = np.asarray(X)
            if weights is None:
                k = distances.shape[1]
                weights = np.full(k, 1.0 / k)
            else:
                distances = distances[:, :k]
        else:
            distances = self.knn.transform(X)
        distances = distances ** q
        dtm = (distances * weights).sum(-1)
        if self.normalize:
            dtm /= (np.arange(1, k + 1) ** (q / dim) * weights).sum()
        density = dtm ** (-dim / q)
        if self.normalize:
            import math

            if self.params["metric"] == "precomputed":
                self.n_samples = len(X[0])
            # Volume of d-ball
            Vd = math.pi ** (dim / 2) / math.gamma(dim / 2 + 1)
            density /= self.n_samples * Vd
        return density
        # We compute too many powers, 1/p in knn then q in dtm, d/q in dtm then whatever in the caller.
        # Add option to skip the final root?
