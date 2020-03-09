import numpy
from ._tomato import *

# The fit/predict interface is not so well suited...
# TODO: option for a faster, weaker (probabilistic) knn


class Tomato:
    """
    Clustering

    This clustering algorithm needs a neighborhood graph on the points, and an estimation of the density at each point. A few possible graph constructions and density estimators are provided for convenience, but it is perfectly natural to provide your own. In particular, we do not provide anything specific to cluster pixels on images yet.

    Attributes
    ----------
    n_clusters_: int
        The number of clusters. Writing to it automatically adjusts labels_.
    merge_threshold_: float
        minimum prominence of a cluster so it doesn't get merged. Writing to it automatically adjusts labels_.
    n_leaves_: int
        number of leaves (unstable clusters) in the hierarchical tree
    leaf_labels_: ndarray of shape (n_samples)
        cluster labels for each point, at the very bottom of the hierarchy
    labels_: ndarray of shape (n_samples)
        cluster labels for each point, after merging
    diagram_: ndarray of shape (n_leaves_,2)
        persistence diagram (only the finite points)
    children_: ndarray of shape (n_leaves_-1,2)
        The children of each non-leaf node. Values less than n_leaves_ correspond to leaves of the tree. A node i greater than or equal to n_leaves_ is a non-leaf node and has children children_[i - n_leaves_]. Alternatively at the i-th iteration, children[i][0] and children[i][1] are merged to form node n_leaves_ + i
    params_: dict
        Parameters like input_type, metric, etc
    """

    # Not documented for now, because I am not sure how useful it is.
    #    max_density_per_cc_: ndarray of shape (n_connected_components)
    #        maximum of the density function on each connected component

    def __init__(
        self,
        input_type="points",
        metric=None,
        graph_type="knn",
        density_type="logDTM",
        n_clusters=None,
        merge_threshold=None,
        #       eliminate_threshold=None,
        #           eliminate_threshold (float): minimum max weight of a cluster so it doesn't get eliminated
        **params
    ):
        """
        Each parameter has a corresponding attribute, like self.merge_threshold_, that can be changed later.

        Args:
            input_type (str): 'points', 'distance_matrix' or 'neighbors'.
            metric (None|Callable): If None, use Minkowski of parameter p.
            graph_type (str): 'manual', 'knn' or 'radius'. Ignored if input_type is 'neighbors'.
            density_type (str): 'manual', 'DTM', 'logDTM' or 'kde'.
            kde_params (dict): if density_type is 'kde', additional parameters passed directly to sklearn.neighbors.KernelDensity.
            k (int): number of neighbors for a knn graph (including the vertex itself). Defaults to 10.
            k_DTM (int): number of neighbors for the DTM density estimation (including the vertex itself). Defaults to k.
            r (float): size of a neighborhood if graph_type is 'radius'
            eps (float): approximation factor when computing nearest neighbors (currently ignored with a GPU)
            gpu (bool): enable use of CUDA (through pykeops) to compute k nearest neighbors. This is useful when the dimension becomes large (10+) but the number of points remains low (less than a million).
            n_clusters (int): number of clusters requested. Defaults to ???
            merge_threshold (float): minimum prominence of a cluster so it doesn't get merged.
            symmetrize_graph (bool): whether we should add edges to make the neighborhood graph symmetric. This can be useful with k-NN for small k. Defaults to false.
            p (float): norm L^p on input points (numpy.inf is supported without gpu). Defaults to 2.
            p_DTM (float): order used to compute the distance to measure. Defaults to the dimension, or 2 if input_type is 'distance_matrix'.
            n_jobs (int): Number of jobs to schedule for parallel processing of nearest neighbors on the CPU. If -1 is given all processors are used. Default: 1.
        """
        # Should metric='precomputed' mean input_type='distance_matrix'?
        # Should we be able to pass metric='minkowski' (what None does currently)?
        self.input_type_ = input_type
        self.metric_ = metric
        self.graph_type_ = graph_type
        self.density_type_ = density_type
        self.params_ = params
        self.__n_clusters = n_clusters
        self.__merge_threshold = merge_threshold
        # self.eliminate_threshold_ = eliminate_threshold
        if n_clusters and merge_threshold:
            raise ValueError("Cannot specify both a merge threshold and a number of clusters")

    def fit(self, X, y=None, weights=None):
        """
        Args:
            X ((n,d)-array of float|(n,n)-array of float|Iterable[Iterable[int]]): coordinates of the points, or distance_matrix (full, not just a triangle), or list of neighbors for each point (points are represented by their index, starting from 0), according to input_type.
            weights (ndarray of shape (n_samples)): if density_type is 'manual', a density estimate at each point
        """
        # TODO: First detect if this is a new call with the same data (only threshold changed?)
        # TODO: less code duplication (subroutines?), less spaghetti, but don't compute neighbors twice if not needed. Clear error message for missing or contradictory parameters.
        if weights is not None:
            density_type = "manual"
        else:
            density_type = self.density_type_
            assert density_type != "manual"
            if density_type == "manual":
                raise ValueError("If density_type is 'manual', you must provide weights to fit()")

        input_type = self.input_type_
        if input_type == "points" and self.metric_:
            from sklearn.metrics import pairwise_distances

            X = pairwise_distances(X, metric=self.metric_, n_jobs=self.params_.get("n_jobs"))
            input_type = "distance_matrix"

        if input_type == "distance_matrix" and self.graph_type_ == "radius":
            X = numpy.array(X)
            r = self.params_["r"]
            self.neighbors_ = [numpy.nonzero(l <= r) for l in X]

        if input_type == "distance_matrix" and self.graph_type_ == "knn":
            k = self.params_["k"]
            self.neighbors_ = numpy.argpartition(X, k - 1)[:, 0:k]

        if input_type == "neighbors":
            self.neighbors_ = X
            assert density_type == "manual"

        if input_type == "points" and self.graph_type_ == "knn" and self.density_type_ in {"DTM", "logDTM"}:
            self.points_ = X
            q = self.params_.get("p_DTM", len(X[0]))
            p = self.params_.get("p", 2)
            k = self.params_.get("k", 10)
            k_DTM = self.params_.get("k_DTM", k)
            kmax = max(k, k_DTM)
            if self.params_.get("gpu"):
                import torch
                from pykeops.torch import LazyTensor

                # 'float64' is slow except on super expensive GPUs. Allow it with some param?
                XX = torch.tensor(self.points_, dtype=torch.float32)
                if p == numpy.inf:
                    # Requires a version of pykeops strictly more recent than 1.3
                    dd, nn = (
                        (LazyTensor(XX[:, None, :]) - LazyTensor(XX[None, :, :]))
                        .abs()
                        .max(-1)
                        .Kmin_argKmin(kmax, dim=1)
                    )
                elif p == 2:  # Any even integer?
                    dd, nn = (
                        ((LazyTensor(XX[:, None, :]) - LazyTensor(XX[None, :, :])) ** p)
                        .sum(-1)
                        .Kmin_argKmin(kmax, dim=1)
                    )
                else:
                    dd, nn = (
                        ((LazyTensor(XX[:, None, :]) - LazyTensor(XX[None, :, :])).abs() ** p)
                        .sum(-1)
                        .Kmin_argKmin(kmax, dim=1)
                    )
                if k < kmax:
                    nn = nn[:, 0:k]
                if k_DTM < kmax:
                    dd = dd[:, 0:k_DTM]
                assert q != numpy.inf  # for now
                if p != numpy.inf:
                    qp = float(q) / p
                else:
                    qp = q
                if qp != 1:
                    dd = dd ** qp
                weights = dd.sum(-1)

                # Back to the CPU. Not sure this is necessary, or the right way to do it.
                weights = numpy.array(weights)
                self.neighbors_ = numpy.array(nn)
            else:  # CPU
                from scipy.spatial import cKDTree

                kdtree = cKDTree(self.points_)
                qargs = {k: v for k, v in self.params_.items() if k in {"eps", "n_jobs"}}
                dd, self.neighbors_ = kdtree.query(self.points_, k=kmax, p=p, **qargs)
                if k < kmax:
                    self.neighbors_ = self.neighbors_[:, 0:k]
                if k_DTM < kmax:
                    dd = dd[:, 0:k_DTM]
                # weights = numpy.linalg.norm(dd, axis=1, ord=q)
                weights = (dd ** q).sum(-1)

            if self.density_type_ == "DTM":
                # We ignore constant factors, which don't matter for
                # clustering, although they do change thresholds
                dim = len(self.points_[0])
                weights = weights ** (-dim / q)
            else:
                # We ignore exponents, which become constant factors with log
                weights = -numpy.log(weights)

        if input_type == "points" and self.graph_type_ == "knn" and self.density_type_ not in {"DTM", "logDTM"}:
            self.points_ = X
            p = self.params_.get("p", 2)
            k = self.params_.get("k", 10)
            if self.params_.get("gpu"):
                import torch
                from pykeops.torch import LazyTensor

                # 'float64' is slow except on super expensive GPUs. Allow it with some param?
                XX = torch.tensor(self.points_, dtype=torch.float32)
                if p == numpy.inf:
                    # Requires a version of pykeops strictly more recent than 1.3
                    nn = (LazyTensor(XX[:, None, :]) - LazyTensor(XX[None, :, :])).abs().max(-1).argKmin(k, dim=1)
                elif p == 2:  # Any even integer?
                    nn = ((LazyTensor(XX[:, None, :]) - LazyTensor(XX[None, :, :])) ** p).sum(-1).argKmin(k, dim=1)
                else:
                    nn = (
                        ((LazyTensor(XX[:, None, :]) - LazyTensor(XX[None, :, :])).abs() ** p).sum(-1).argKmin(k, dim=1)
                    )
                # Back to the CPU. Not sure this is necessary, or the right way to do it.
                self.neighbors_ = numpy.array(nn)
            else:  # CPU
                from scipy.spatial import cKDTree

                kdtree = cKDTree(self.points_)
                # FIXME 'p'
                qargs = {k: v for k, v in self.params_.items() if k in {"eps", "n_jobs"}}
                _, self.neighbors_ = kdtree.query(self.points_, k=k, p=p, **qargs)

        if input_type == "points" and self.graph_type_ != "knn" and self.density_type_ in {"DTM", "logDTM"}:
            self.points_ = X
            q = self.params_.get("p_DTM", len(X[0]))
            p = self.params_.get("p", 2)
            k = self.params_.get("k", 10)
            k_DTM = self.params_.get("k_DTM", k)
            if self.params_.get("gpu"):
                import torch
                from pykeops.torch import LazyTensor

                # 'float64' is slow except on super expensive GPUs. Allow it with some param?
                XX = torch.tensor(self.points_, dtype=torch.float32)
                if p == numpy.inf:
                    assert False  # Not supported???
                    dd = (LazyTensor(XX[:, None, :]) - LazyTensor(XX[None, :, :])).abs().max(-1).Kmin(k_DTM, dim=1)
                elif p == 2:  # Any even integer?
                    dd = ((LazyTensor(XX[:, None, :]) - LazyTensor(XX[None, :, :])) ** p).sum(-1).Kmin(k_DTM, dim=1)
                else:
                    dd = (
                        ((LazyTensor(XX[:, None, :]) - LazyTensor(XX[None, :, :])).abs() ** p)
                        .sum(-1)
                        .Kmin(k_DTM, dim=1)
                    )
                assert q != numpy.inf  # for now
                if p != numpy.inf:
                    qp = float(q) / p
                else:
                    qp = q
                if qp != 1:
                    dd = dd ** qp
                weights = dd.sum(-1)
                # **1/q is a waste of time, whether we take another **-.25 or a log

                # Back to the CPU. Not sure this is necessary, or the right way to do it.
                weights = numpy.array(weights)
            else:  # CPU
                from scipy.spatial import cKDTree

                kdtree = cKDTree(self.points_)
                qargs = {k: v for k, v in self.params_.items() if k in {"eps", "n_jobs"}}
                dd, _ = kdtree.query(self.points_, k=k_DTM, p=p, **qargs)
                # weights = numpy.linalg.norm(dd, axis=1, ord=q)
                weights = (dd ** q).sum(-1)

            if self.density_type_ == "DTM":
                dim = len(self.points_[0])
                weights = weights ** (-dim / q)
            else:
                weights = -numpy.log(weights)

        if input_type == "distance_matrix" and self.density_type_ in {"DTM", "logDTM"}:
            q = self.params_.get("p_DTM", 2)
            X = numpy.array(X)
            k = self.params_.get("k_DTM")
            if not k:
                k = self.params_["k"]
            weights = (numpy.partition(X, k - 1)[:, 0:k] ** q).sum(-1)
            if self.density_type_ == "DTM":
                dim = len(self.points_[0])
                weights = weights ** (-dim / q)
            else:
                weights = -numpy.log(weights)

        if self.density_type_ == "kde":
            # FIXME: replace most assert with raise ValueError("blabla")
            assert input_type == "points"
            kde_params = self.params_.get("kde_params", dict())
            from sklearn.neighbors import KernelDensity

            weights = KernelDensity(**kde_params).fit(X).score_samples(X)

        if self.params_.get("symmetrize_graph"):
            self.neighbors_ = [set(line) for line in self.neighbors_]
            for i, line in enumerate(self.neighbors_):
                line.discard(i)
                for j in line:
                    self.neighbors_[j].add(i)

        self.weights_ = weights  # TODO remove
        self.leaf_labels_, self.children_, self.diagram_, self.max_density_per_cc_ = doit(
            list(self.neighbors_), weights
        )
        self.n_leaves_ = len(self.max_density_per_cc_) + len(self.children_)
        assert self.leaf_labels_.max() + 1 == len(self.max_density_per_cc_) + len(self.children_)
        if self.__merge_threshold:
            assert not self.__n_clusters
            self.__n_clusters = numpy.count_nonzero(
                self.diagram_[:, 1] - self.diagram_[:, 0] > self.__merge_threshold
            ) + len(self.max_density_per_cc_)
        if self.__n_clusters:
            renaming = merge(self.children_, self.n_leaves_, self.__n_clusters)
            self.labels_ = renaming[self.leaf_labels_]
        else:
            self.labels_ = self.leaf_labels_
            self.__n_clusters = self.n_leaves_
        return self

    def fit_predict(self, X, y=None, weights=None):
        """
        Equivalent to fit(), and returns the `labels_`.
        """
        return self.fit(X, y, weights).labels_

    # TODO: add argument k or threshold? Have a version where you can click and it shows the line and the corresponding k?
    def plot_diagram(self):
        """
        """
        import matplotlib.pyplot as plt

        plt.plot(self.diagram_[:, 0], self.diagram_[:, 1], "ro")
        l = self.diagram_[:, 1].min()
        r = max(self.diagram_[:, 0].max(), self.max_density_per_cc_.max())
        plt.plot([l, r], [l, r])
        plt.plot(
            self.max_density_per_cc_, numpy.full(self.max_density_per_cc_.shape, 1.1 * l - 0.1 * r), "ro", color="green"
        )
        plt.show()

    #    def predict(self, X):
    #        # X had better be the same as in fit()
    #        return self.labels_

    # Use set_params instead?
    @property
    def n_clusters_(self):
        return self.__n_clusters

    @n_clusters_.setter
    def n_clusters_(self, n_clusters):
        if n_clusters == self.__n_clusters:
            return
        self.__n_clusters = n_clusters
        self.__merge_threshold = None
        if hasattr(self, "leaf_labels_"):
            renaming = merge(self.children_, self.n_leaves_, self.__n_clusters)
            self.labels_ = renaming[self.leaf_labels_]

    @property
    def merge_threshold_(self):
        return self.__merge_threshold

    @merge_threshold_.setter
    def merge_threshold_(self, merge_threshold):
        if merge_threshold == self.__merge_threshold:
            return
        if hasattr(self, "leaf_labels_"):
            self.n_clusters_ = numpy.count_nonzero(self.diagram_[:, 1] - self.diagram_[:, 0] > merge_threshold) + len(
                self.max_density_per_cc_
            )
        else:
            self.__n_clusters = None
        self.__merge_threshold = merge_threshold


if __name__ == "__main__":
    import sys

    a = [(1, 2), (1.1, 1.9), (0.9, 1.8), (10, 0), (10.1, 0.05), (10.2, -0.1), (5.4, 0)]
    a = numpy.random.rand(500, 2)
    t = Tomato(
        input_type="points",
        metric="Euclidean",
        graph_type="knn",
        density_type="DTM",
        n_clusters=2,
        k=4,
        n_jobs=-1,
        eps=0.05,
    )
    t.fit(a)
    # print("neighbors\n",t.neighbors_)
    # print()
    # print("weights\n",t.weights_)
    # print()
    # print("diagram\n",t.diagram_)
    # print()
    print("max\n", t.max_density_per_cc_, file=sys.stderr)
    # print()
    print("leaf labels\n", t.leaf_labels_)
    # print()
    print("labels\n", t.labels_)
    # print()
    print("children\n", t.children_)
    # print()
    t.n_clusters_ = 2
    print("labels\n", t.labels_)
    t.plot_diagram()
