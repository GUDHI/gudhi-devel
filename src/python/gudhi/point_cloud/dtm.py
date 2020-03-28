# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Marc Glisse
#
# Copyright (C) 2020 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from .knn import KNN


class DTM:
    """
    Class to compute the distance to the empirical measure defined by a point set.
    """

    def __init__(self, k, q=None, **kwargs):
        """
        Args:
            k (int): number of neighbors (possibly including the point itself).
            q (float): order used to compute the distance to measure. Defaults to the dimension, or 2 if metric is "neighbors" or "distance_matrix".
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
            # KNN gives sorted distances, which is unnecessary here.
            # Maybe add a parameter to say we don't need sorting?
            self.knn = KNN(self.k, return_index=False, return_distance=True, **self.params)
            self.knn.fit(X)
        return self

    def transform(self, X):
        """
        Args:
            X (numpy.array): coordinates for query points, or distance matrix if metric is "precomputed", or distances to the k nearest neighbors if metric is "neighbors" (if the array has more than k columns, the remaining ones are ignored).
        """
        q = self.q
        if q is None:
            if self.params["metric"] in {"neighbors", "precomputed"}:
                q = 2
            else:
                q = len(X[0])
        if self.params["metric"] == "neighbors":
            distances = X[:, : self.k]
        else:
            distances = self.knn.transform(X)
        distances = distances ** q
        dtm = distances.sum(-1) / self.k
        dtm = dtm ** (1.0 / q)
        # We compute too many powers, 1/p in knn then q in dtm, 1/q in dtm then q or some log in the caller.
        # Add option to skip the final root?
        return dtm
