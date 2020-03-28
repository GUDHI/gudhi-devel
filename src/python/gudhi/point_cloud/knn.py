# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Marc Glisse
#
# Copyright (C) 2020 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

import numpy


class KNN:
    """
    Class wrapping several implementations for computing the k nearest neighbors in a point set.
    """

    def __init__(self, k, return_index=True, return_distance=False, metric="euclidean", **kwargs):
        """
        Args:
            k (int): number of neighbors (possibly including the point itself).
            return_index (bool): if True, return the index of each neighbor.
            return_distance (bool): if True, return the distance to each neighbor.
            implementation (str): choice of the library that does the real work.

                * 'keops' for a brute-force, CUDA implementation through pykeops. Useful when the dimension becomes large (10+) but the number of points remains low (less than a million). Only "minkowski" and its aliases are supported.
                * 'ckdtree' for scipy's cKDTree. Only "minkowski" and its aliases are supported.
                * 'sklearn' for scikit-learn's NearestNeighbors. Note that this provides in particular an option algorithm="brute".
                * 'hnsw' for hnswlib.Index. It can be very fast but does not provide guarantees. Only supports "euclidean" for now.
                * None will try to select a sensible one (scipy if possible, scikit-learn otherwise).
            metric (str): see `sklearn.neighbors.NearestNeighbors`.
            eps (float): relative error when computing nearest neighbors with the cKDTree.
            p (float): norm L^p on input points (including numpy.inf) if metric is "minkowski". Defaults to 2.
            n_jobs (int): number of jobs to schedule for parallel processing of nearest neighbors on the CPU.
                If -1 is given all processors are used. Default: 1.
            kwargs: additional parameters are forwarded to the backends.
        """
        self.k = k
        self.return_index = return_index
        self.return_distance = return_distance
        self.metric = metric
        self.params = kwargs
        # canonicalize
        if metric == "euclidean":
            self.params["p"] = 2
            self.metric = "minkowski"
        elif metric == "manhattan":
            self.params["p"] = 1
            self.metric = "minkowski"
        elif metric == "chebyshev":
            self.params["p"] = numpy.inf
            self.metric = "minkowski"
        elif metric == "minkowski":
            self.params["p"] = kwargs.get("p", 2)
        if self.params.get("implementation") in {"keops", "ckdtree"}:
            assert self.metric == "minkowski"
        if self.params.get("implementation") == "hnsw":
            assert self.metric == "minkowski" and self.params["p"] == 2
        if not self.params.get("implementation"):
            if self.metric == "minkowski":
                self.params["implementation"] = "ckdtree"
            else:
                self.params["implementation"] = "sklearn"

    def fit_transform(self, X, y=None):
        return self.fit(X).transform(X)

    def fit(self, X, y=None):
        """
        Args:
            X (numpy.array): coordinates for reference points.
        """
        self.ref_points = X
        if self.params["implementation"] == "ckdtree":
            # sklearn could handle this, but it is much slower
            from scipy.spatial import cKDTree

            self.kdtree = cKDTree(X)

        if self.params["implementation"] == "sklearn" and self.metric != "precomputed":
            # FIXME: sklearn badly handles "precomputed"
            from sklearn.neighbors import NearestNeighbors

            nargs = {
                k: v for k, v in self.params.items() if k in {"p", "n_jobs", "metric_params", "algorithm", "leaf_size"}
            }
            self.nn = NearestNeighbors(self.k, metric=self.metric, **nargs)
            self.nn.fit(X)

        if self.params["implementation"] == "hnsw":
            import hnswlib

            self.graph = hnswlib.Index("l2", len(X[0]))  # Actually returns squared distances
            self.graph.init_index(
                len(X), **{k: v for k, v in self.params.items() if k in {"ef_construction", "M", "random_seed"}}
            )
            n = self.params.get("num_threads")
            if n is None:
                n = self.params.get("n_jobs", 1)
                self.params["num_threads"] = n
            self.graph.add_items(X, num_threads=n)

        return self

    def transform(self, X):
        """
        Args:
            X (numpy.array): coordinates for query points, or distance matrix if metric is "precomputed".
        """
        metric = self.metric
        k = self.k

        if metric == "precomputed":
            # scikit-learn could handle that, but they insist on calling fit() with an unused square array, which is too unnatural.
            X = numpy.array(X)
            if self.return_index:
                neighbors = numpy.argpartition(X, k - 1)[:, 0:k]
                distances = numpy.take_along_axis(X, neighbors, axis=-1)
                ngb_order = numpy.argsort(distances, axis=-1)
                neighbors = numpy.take_along_axis(neighbors, ngb_order, axis=-1)
                if self.return_distance:
                    distances = numpy.take_along_axis(distances, ngb_order, axis=-1)
                    return neighbors, distances
                else:
                    return neighbors
            if self.return_distance:
                distances = numpy.partition(X, k - 1)[:, 0:k]
                # partition is not guaranteed to sort the lower half, although it often does
                distances.sort(axis=-1)
                return distances
            return None

        if self.params["implementation"] == "hnsw":
            ef = self.params.get("ef")
            if ef is not None:
                self.graph.set_ef(ef)
            neighbors, distances = self.graph.knn_query(X, k, num_threads=self.params["num_threads"])
            # The k nearest neighbors are always sorted. I couldn't find it in the doc, but the code calls searchKnn,
            # which returns a priority_queue, and then fills the return array backwards with top/pop on the queue.
            if self.return_index:
                if self.return_distance:
                    return neighbors, numpy.sqrt(distances)
                else:
                    return neighbors
            if self.return_distance:
                return numpy.sqrt(distances)
            return None

        if self.params["implementation"] == "keops":
            import torch
            from pykeops.torch import LazyTensor

            # 'float64' is slow except on super expensive GPUs. Allow it with some param?
            XX = torch.tensor(X, dtype=torch.float32)
            if X is self.ref_points:
                YY = XX
            else:
                YY = torch.tensor(self.ref_points, dtype=torch.float32)

            p = self.params["p"]
            if p == numpy.inf:
                # Requires pykeops 1.4 or later
                mat = (LazyTensor(XX[:, None, :]) - LazyTensor(YY[None, :, :])).abs().max(-1)
            elif p == 2:  # Any even integer?
                mat = ((LazyTensor(XX[:, None, :]) - LazyTensor(YY[None, :, :])) ** p).sum(-1)
            else:
                mat = ((LazyTensor(XX[:, None, :]) - LazyTensor(YY[None, :, :])).abs() ** p).sum(-1)

            if self.return_index:
                if self.return_distance:
                    distances, neighbors = mat.Kmin_argKmin(k, dim=1)
                    if p != numpy.inf:
                        distances = distances ** (1.0 / p)
                    return neighbors, distances
                else:
                    neighbors = mat.argKmin(k, dim=1)
                    return neighbors
            if self.return_distance:
                distances = mat.Kmin(k, dim=1)
                if p != numpy.inf:
                    distances = distances ** (1.0 / p)
                return distances
            return None
        # FIXME: convert everything back to numpy arrays or not?

        if self.params["implementation"] == "ckdtree":
            qargs = {key: val for key, val in self.params.items() if key in {"p", "eps", "n_jobs"}}
            distances, neighbors = self.kdtree.query(X, k=self.k, **qargs)
            if self.return_index:
                if self.return_distance:
                    return neighbors, distances
                else:
                    return neighbors
            if self.return_distance:
                return distances
            return None

        assert self.params["implementation"] == "sklearn"
        if self.return_distance:
            distances, neighbors = self.nn.kneighbors(X, return_distance=True)
            if self.return_index:
                return neighbors, distances
            else:
                return distances
        if self.return_index:
            neighbors = self.nn.kneighbors(X, return_distance=False)
            return neighbors
        return None
