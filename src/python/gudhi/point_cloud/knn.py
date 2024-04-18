# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Marc Glisse
#
# Copyright (C) 2020 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

import numpy
import warnings

# TODO: https://github.com/facebookresearch/faiss

__author__ = "Marc Glisse"
__copyright__ = "Copyright (C) 2020 Inria"
__license__ = "MIT"


class KNearestNeighbors:
    """
    Class wrapping several implementations for computing the k nearest neighbors in a point set.

    :Requires: `PyKeOps <installation.html#pykeops>`_, `SciPy <installation.html#scipy>`_,
        `Scikit-learn <installation.html#scikit-learn>`_, and/or `Hnswlib <installation.html#hnswlib>`_
        in function of the selected `implementation`.
    """

    def __init__(self, k, return_index=True, return_distance=False, metric="euclidean", **kwargs):
        """
        Args:
            k (int): number of neighbors (possibly including the point itself).
            return_index (bool): if True, return the index of each neighbor.
            return_distance (bool): if True, return the distance to each neighbor.
            implementation (str): choice of the library that does the real work.

                * 'keops' for a brute-force, CUDA implementation through pykeops. Useful when the dimension becomes
                  large (10+) but the number of points remains low (less than a million). Only "minkowski" and its
                  aliases are supported.
                * 'ckdtree' for scipy's cKDTree. Only "minkowski" and its aliases are supported.
                * 'sklearn' for scikit-learn's NearestNeighbors. Note that this provides in particular an option
                  algorithm="brute".
                * 'hnsw' for hnswlib.Index. It can be very fast but does not provide guarantees. Only supports
                  "euclidean" for now.
                * None will try to select a sensible one (scipy if possible, scikit-learn otherwise).
            metric (str): see `sklearn.neighbors.NearestNeighbors`.
            eps (float): relative error when computing nearest neighbors with the cKDTree.
            p (float): norm L^p on input points (including numpy.inf) if metric is "minkowski". Defaults to 2.
            n_jobs (int): number of jobs to schedule for parallel processing of nearest neighbors on the CPU.
                If -1 is given all processors are used. Default: 1.
            sort_results (bool): if True, then distances and indices of each point are
                sorted on return, so that the first column contains the closest points.
                Otherwise, neighbors are returned in an arbitrary order. Defaults to True.
            enable_autodiff (bool): if the input is a torch.tensor or tensorflow.Tensor, this
                instructs the function to compute distances in a way that works with automatic differentiation.
                This is experimental, not supported for all metrics, and requires the package EagerPy.
                Defaults to False.
            kwargs: additional parameters are forwarded to the backends.
        """
        if k < 1:
            raise ValueError(f"Expected number of neighbors (aka. 'k') > 0. Got {k}")
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
        if not return_distance:
            self.params["enable_autodiff"] = False

    def fit_transform(self, X, y=None):
        return self.fit(X).transform(X)

    def fit(self, X, y=None):
        """
        Args:
            X (numpy.array): coordinates for reference points.
        """
        if self.k > len(X):
            raise ValueError(
                f"Expected number of neighbors (aka. 'k') <= number of samples, but k={self.k} and number of samples={len(X)}"
            )
        self.ref_points = X
        if self.params.get("enable_autodiff", False):
            import eagerpy as ep

            X = ep.astensor(X)
            if self.params["implementation"] != "keops" or not isinstance(X, ep.PyTorchTensor):
                # I don't know a clever way to reuse a GPU tensor from tensorflow in pytorch
                # without copying to/from the CPU.
                X = X.numpy()
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
            self.nn = NearestNeighbors(n_neighbors=self.k, metric=self.metric, **nargs)
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

        Returns:
            numpy.array: if return_index, an array of shape (len(X), k) with the indices (in the argument
            of :func:`fit`) of the k nearest neighbors to the points of X. If return_distance, an array of the
            same shape with the distances to those neighbors. If both, a tuple with the two arrays, in this order.
        """
        if self.params.get("enable_autodiff", False):
            # pykeops does not support autodiff for kmin yet, but when it does in the future,
            # we may want a special path.
            import eagerpy as ep

            save_return_index = self.return_index
            self.return_index = True
            self.return_distance = False
            self.params["enable_autodiff"] = False
            try:
                newX = ep.astensor(X)
                if self.params["implementation"] != "keops" or (
                    not isinstance(newX, ep.PyTorchTensor) and not isinstance(newX, ep.NumPyTensor)
                ):
                    newX = newX.numpy()
                else:
                    newX = newX.raw
                neighbors = self.transform(newX)
            finally:
                self.return_index = save_return_index
                self.return_distance = True
                self.params["enable_autodiff"] = True
            # We can implement more later as needed
            assert self.metric == "minkowski"
            p = self.params["p"]
            Y = ep.astensor(self.ref_points)
            neighbor_pts = Y[neighbors,]
            diff = neighbor_pts - X[:, None, :]
            if isinstance(diff, ep.PyTorchTensor):
                # https://github.com/jonasrauber/eagerpy/issues/6
                distances = ep.astensor(diff.raw.norm(p, -1))
            else:
                distances = diff.norms.lp(p, -1)
            if self.return_index:
                return neighbors, distances.raw
            else:
                return distances.raw

        metric = self.metric
        k = self.k

        if metric == "precomputed":
            # scikit-learn could handle that, but they insist on calling fit() with an unused square array
            # which is too unnatural.
            if self.return_index:
                n_jobs = self.params.get("n_jobs", 1)
                # Supposedly numpy can be compiled with OpenMP and handle this, but nobody does that?!
                if n_jobs == 1:
                    neighbors = numpy.argpartition(X, k - 1)[:, 0:k]
                    if self.params.get("sort_results", True):
                        X = numpy.take_along_axis(X, neighbors, axis=-1)
                        ngb_order = numpy.argsort(X, axis=-1)
                        neighbors = numpy.take_along_axis(neighbors, ngb_order, axis=-1)
                    else:
                        ngb_order = neighbors
                    if self.return_distance:
                        distances = numpy.take_along_axis(X, ngb_order, axis=-1)
                        return neighbors, distances
                    else:
                        return neighbors
                else:
                    from joblib import Parallel, delayed, effective_n_jobs
                    from sklearn.utils import gen_even_slices

                    slices = gen_even_slices(len(X), effective_n_jobs(n_jobs))
                    parallel = Parallel(prefer="threads", n_jobs=n_jobs)
                    if self.params.get("sort_results", True):

                        def func(M):
                            neighbors = numpy.argpartition(M, k - 1)[:, 0:k]
                            Y = numpy.take_along_axis(M, neighbors, axis=-1)
                            ngb_order = numpy.argsort(Y, axis=-1)
                            return numpy.take_along_axis(neighbors, ngb_order, axis=-1)

                    else:

                        def func(M):
                            return numpy.argpartition(M, k - 1)[:, 0:k]

                    neighbors = numpy.concatenate(parallel(delayed(func)(X[s]) for s in slices))
                    if self.return_distance:
                        distances = numpy.take_along_axis(X, neighbors, axis=-1)
                        return neighbors, distances
                    else:
                        return neighbors
            if self.return_distance:
                n_jobs = self.params.get("n_jobs", 1)
                if n_jobs == 1:
                    distances = numpy.partition(X, k - 1)[:, 0:k]
                    if self.params.get("sort_results"):
                        # partition is not guaranteed to sort the lower half, although it often does
                        distances.sort(axis=-1)
                else:
                    from joblib import Parallel, delayed, effective_n_jobs
                    from sklearn.utils import gen_even_slices

                    if self.params.get("sort_results"):

                        def func(M):
                            # Not partitioning in place, because we should not modify the user's array?
                            r = numpy.partition(M, k - 1)[:, 0:k]
                            r.sort(axis=-1)
                            return r

                    else:
                        func = lambda M: numpy.partition(M, k - 1)[:, 0:k]
                    slices = gen_even_slices(len(X), effective_n_jobs(n_jobs))
                    parallel = Parallel(prefer="threads", n_jobs=n_jobs)
                    distances = numpy.concatenate(parallel(delayed(func)(X[s]) for s in slices))
                return distances
            return None

        if self.params["implementation"] == "hnsw":
            ef = self.params.get("ef")
            if ef is not None:
                self.graph.set_ef(ef)
            neighbors, distances = self.graph.knn_query(X, k, num_threads=self.params["num_threads"])
            with warnings.catch_warnings():
                if not (numpy.all(numpy.isfinite(distances))):
                    warnings.warn("Overflow/infinite value encountered while computing 'distances'", RuntimeWarning)
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
            XX = torch.as_tensor(X, dtype=torch.float32)
            if X is self.ref_points:
                YY = XX
            else:
                YY = torch.as_tensor(self.ref_points, dtype=torch.float32)
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
                    with warnings.catch_warnings():
                        if not (torch.isfinite(distances).all()):
                            warnings.warn("Overflow/infinite value encountered while computing 'distances'",
                                          RuntimeWarning)
                    if p != numpy.inf:
                        distances = distances ** (1.0 / p)
                    return neighbors, distances
                else:
                    neighbors = mat.argKmin(k, dim=1)
                    return neighbors
            if self.return_distance:
                distances = mat.Kmin(k, dim=1)
                with warnings.catch_warnings():
                    if not (torch.isfinite(distances).all()):
                        warnings.warn("Overflow/infinite value encountered while computing 'distances'",
                                      RuntimeWarning)
                if p != numpy.inf:
                    distances = distances ** (1.0 / p)
                return distances
            return None

        if self.params["implementation"] == "ckdtree":
            qargs = {key: val for key, val in self.params.items() if key in {"p", "eps"}}
            # SciPy renamed n_jobs to workers
            qargs["workers"] = self.params.get("workers") or self.params.get("n_jobs") or 1
            distances, neighbors = self.kdtree.query(X, k=self.k, **qargs)
            if k == 1:
                # SciPy decided to squeeze the last dimension for k=1
                distances = distances[:, None]
                neighbors = neighbors[:, None]
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
