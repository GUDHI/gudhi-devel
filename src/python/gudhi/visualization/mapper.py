"""
@author: Mathieu Carriere
All rights reserved
"""

import numpy as np
import itertools
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import colorsys

from networkx                import cycle_basis
from scipy.sparse.csgraph    import dijkstra, shortest_path, connected_components
from sklearn.base            import BaseEstimator, TransformerMixin
from sklearn.preprocessing   import LabelEncoder
from sklearn.cluster         import DBSCAN, AgglomerativeClustering
from sklearn.metrics         import pairwise_distances
from scipy.spatial.distance  import directed_hausdorff
from scipy.sparse            import csgraph
from sklearn.neighbors       import KernelDensity, kneighbors_graph, radius_neighbors_graph, NearestNeighbors
from scipy.stats             import ks_2samp
from ..simplex_tree          import SimplexTree
from ..nerve_gic             import GraphInducedComplex


def estimate_scale(X, N=100, inp="point cloud", beta=0., C=10.):
    """
    Compute estimated scale of a point cloud or a distance matrix.

    Parameters:
        X (numpy array of shape (num_points) x (num_coordinates) if point cloud and (num_points) x (num_points) if distance matrix): input point cloud or distance matrix.
        N (int): subsampling iterations (default 100). See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
        inp (string): either "point cloud" or "distance matrix". Type of input data (default "point cloud").
        beta (double): exponent parameter (default 0.). See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
        C (double): constant parameter (default 10.). See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.

    Returns:
        delta (double): estimated scale that can be used with eg agglomerative clustering.
    """
    num_pts = X.shape[0]
    delta, m = 0., int(  num_pts / np.exp((1+beta) * np.log(np.log(num_pts)/np.log(C)))  )
    for _ in range(N):
        subpop = np.random.choice(num_pts, size=m, replace=False)
        if inp == "point cloud":
            d, _, _ = directed_hausdorff(X, X[subpop,:])
        if inp == "distance matrix":
            d = np.max(np.min(X[:,subpop], axis=1), axis=0)
        delta += d/N
    return delta


        

# implement union-find
def find(i, parents):
    if parents[i] == i:
        return i
    else:
        return find(parents[i], parents)

def union(i, j, parents, f):
    if f[i] <= f[j]:
        parents[j] = i
    else:
        parents[i] = j


class MapperComplex(BaseEstimator, TransformerMixin):
    """
    This is a class for computing Mapper simplicial complexes on point clouds or distance matrices. 
    """
    def __init__(self, filters, filter_bnds, colors, resolutions, gains, inp="point cloud", clustering=DBSCAN(), mask=0, N=100, beta=0., C=10.):
        """
        Constructor for the MapperComplex class.

        Attributes:
            inp (string): either "point cloud" or "distance matrix". Specifies the type of input data.
            filters (numpy array of shape (num_points) x (num_filters)): filters (sometimes called lenses) used to compute the Mapper. Each column of the numpy array defines a scalar function defined on the input points.
            filter_bnds (numpy array of shape (num_filters) x 2): limits of each filter, of the form [[f_1^min, f_1^max], ..., [f_n^min, f_n^max]]. If one of the values is numpy.nan, it can be computed from the points with the fit() method.
            colors (numpy array of shape (num_points) x (num_colors)): functions used to color the nodes of the output Mapper simplicial complex. More specifically, coloring is done by computing the means of these functions on the subpopulations corresponding to each node. It can be the same as filters.
            resolutions (numpy array of shape num_filters containing integers): resolution of each filter, ie number of intervals required to cover each filter image.
            gains (numpy array of shape num_filters containing doubles in [0,1]): gain of each filter, ie overlap percentage of the intervals covering each filter image.
            clustering (class): clustering class (default sklearn.cluster.DBSCAN()). Common clustering classes can be found in the scikit-learn library (such as AgglomerativeClustering for instance).
            mask (int): threshold on the size of the Mapper nodes (default 0). Any node associated to a subpopulation with less than **mask** points will be removed.
            N (int): subsampling iterations (default 100). See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
            beta (double): exponent parameter (default 0.). See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
            C (double): constant parameter (default 10.). See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.

            mapper_ (gudhi SimplexTree): Mapper simplicial complex computed after calling the fit() method
            node_info_ (dictionary): various information associated to the nodes of the Mapper. 
        """
        self.filters, self.filter_bnds, self.resolutions, self.gains, self.colors, self.clustering = filters, filter_bnds, resolutions, gains, colors, clustering
        self.input, self.mask, self.N, self.beta, self.C = inp, mask, N, beta, C

    def get_optimal_parameters_for_agglomerative_clustering(self, X, beta=0., C=10., N=100):
        """
        Compute optimal scale and resolutions for a point cloud or a distance matrix.

        Parameters:
            X (numpy array of shape (num_points) x (num_coordinates) if point cloud and (num_points) x (num_points) if distance matrix): input point cloud or distance matrix.
            beta (double): exponent parameter (default 0.). See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
            C (double): constant parameter (default 10.). See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
            N (int): subsampling iterations (default 100). See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.

        Returns:
            delta (double): optimal scale that can be used with agglomerative clustering.
            resolutions (numpy array of shape (num_filters): optimal resolutions associated to each filter.
        """
        num_pts, num_filt, delta = X.shape[0], self.filters.shape[1], 0
        delta = estimate_scale(X=X, N=N, inp=self.input, C=C, beta=beta)

        pairwise = pairwise_distances(X, metric="euclidean") if self.input == "point cloud" else X
        pairs = np.argwhere(pairwise <= delta)
        num_pairs = pairs.shape[0]
        res = []
        for f in range(num_filt):
            F = self.filters[:,f]
            minf, maxf = np.min(F), np.max(F)
            resf = 0
            for p in range(num_pairs):
                resf = max(resf, abs(F[pairs[p,0]] - F[pairs[p,1]]))
            res.append(int((maxf-minf)/resf))

        return delta, np.array(res)


    def fit(self, X, y=None):
        """
        Fit the MapperComplex class on a point cloud or a distance matrix: compute the Mapper and store it in a simplex tree called mapper_

        Parameters:
            X (numpy array of shape (num_points) x (num_coordinates) if point cloud and (num_points) x (num_points) if distance matrix): input point cloud or distance matrix.
            y (n x 1 array): point labels (unused).
        """
        num_pts, num_filters, num_colors = self.filters.shape[0], self.filters.shape[1], self.colors.shape[1]

        # If some resolutions are not specified, automatically compute them
        if self.resolutions is None or self.clustering is None:
            delta, resolutions = self.get_optimal_parameters_for_agglomerative_clustering(X=X, beta=self.beta, C=self.C, N=self.N)
            if self.clustering is None:
                if self.input == "point cloud":
                    self.clustering = AgglomerativeClustering(n_clusters=None, linkage="single", distance_threshold=delta, affinity="euclidean")  
                else:
                    self.clustering = AgglomerativeClustering(n_clusters=None, linkage="single", distance_threshold=delta, affinity="precomputed")
            if self.resolutions is None:
                self.resolutions = resolutions
                self.resolutions = np.array([int(r) for r in self.resolutions])

        if self.gains is None:
            self.gains = .33 * np.ones(num_filters)
 
        # If some filter limits are unspecified, automatically compute them
        if self.filter_bnds is None: 
            self.filter_bnds = np.hstack([np.min(self.filters, axis=0)[:,np.newaxis], np.max(self.filters, axis=0)[:,np.newaxis]])

        # Initialize attributes
        self.mapper_, self.node_info_ = SimplexTree(), {}

        if np.all(self.gains < .5):
            
            # Compute which points fall in which patch or patch intersections
            interval_inds, intersec_inds = np.empty(self.filters.shape), np.empty(self.filters.shape)
            for i in range(num_filters):
                f, r, g = self.filters[:,i], self.resolutions[i], self.gains[i]
                min_f, max_f = self.filter_bnds[i,0], np.nextafter(self.filter_bnds[i,1], np.inf)
                interval_endpoints, l = np.linspace(min_f, max_f, num=r+1, retstep=True)
                intersec_endpoints = []
                for j in range(1, len(interval_endpoints)-1):
                    intersec_endpoints.append(interval_endpoints[j] - g*l / (2 - 2*g))
                    intersec_endpoints.append(interval_endpoints[j] + g*l / (2 - 2*g))
                interval_inds[:,i] = np.digitize(f, interval_endpoints)
                intersec_inds[:,i] = 0.5 * (np.digitize(f, intersec_endpoints) + 1)

            # Build the binned_data map that takes a patch or a patch intersection and outputs the indices of the points contained in it
            binned_data = {}
            for i in range(num_pts):
                list_preimage = []
                for j in range(num_filters):
                    a, b = interval_inds[i,j], intersec_inds[i,j]
                    list_preimage.append([a])
                    if b == a:
                        list_preimage[j].append(a+1)
                    if b == a-1:
                        list_preimage[j].append(a-1)
                list_preimage = list(itertools.product(*list_preimage))
                for pre_idx in list_preimage:
                    try:
                        binned_data[pre_idx].append(i)
                    except KeyError:
                        binned_data[pre_idx] = [i]

        else:

            # Compute interval endpoints for each filter
            l_int, r_int = [], []
            for i in range(num_filters):
                L, R = [], []
                f, r, g = self.filters[:,i], self.resolutions[i], self.gains[i]
                min_f, max_f = self.filter_bnds[i,0], np.nextafter(self.filter_bnds[i,1], np.inf)
                interval_endpoints, l = np.linspace(min_f, max_f, num=r+1, retstep=True)
                for j in range(len(interval_endpoints)-1):
                    L.append(interval_endpoints[j]   - g*l / (2 - 2*g))
                    R.append(interval_endpoints[j+1] + g*l / (2 - 2*g))
                l_int.append(L)
                r_int.append(R)

            # Build the binned_data map that takes a patch or a patch intersection and outputs the indices of the points contained in it
            binned_data = {}
            for i in range(num_pts):
                list_preimage = []
                for j in range(num_filters):
                    fval = self.filters[i,j]
                    start, end = int(min(np.argwhere(np.array(r_int[j]) >= fval))), int(max(np.argwhere(np.array(l_int[j]) <= fval)))
                    list_preimage.append(list(range(start, end+1)))
                list_preimage = list(itertools.product(*list_preimage))
                for pre_idx in list_preimage:
                    try:
                        binned_data[pre_idx].append(i)
                    except KeyError:
                        binned_data[pre_idx] = [i]

        # Initialize the cover map, that takes a point and outputs the clusters to which it belongs
        cover, clus_base = [[] for _ in range(num_pts)], 0

        # For each patch
        for preimage in binned_data:

            # Apply clustering on the corresponding subpopulation
            idxs = np.array(binned_data[preimage])
            if len(idxs) > 1:
                clusters = self.clustering.fit_predict(X[idxs,:]) if self.input == "point cloud" else self.clustering.fit_predict(X[idxs,:][:,idxs])
            elif len(idxs) == 1:
                clusters = np.array([0])
            else:
                continue

            # Collect various information on each cluster
            num_clus_pre = np.max(clusters) + 1
            for clus_i in range(num_clus_pre):
                node_name = clus_base + clus_i
                subpopulation = idxs[clusters == clus_i]
                self.node_info_[node_name] = {}
                self.node_info_[node_name]["indices"] = subpopulation
                self.node_info_[node_name]["size"] = len(subpopulation)
                self.node_info_[node_name]["colors"] = np.mean(self.colors[subpopulation,:], axis=0)
                self.node_info_[node_name]["patch"] = preimage

            # Update the cover map
            for pt in range(clusters.shape[0]):
                node_name = clus_base + clusters[pt]
                if clusters[pt] != -1 and self.node_info_[node_name]["size"] >= self.mask:
                    cover[idxs[pt]].append(node_name)

            clus_base += np.max(clusters) + 1

        # Insert the simplices of the Mapper complex 
        for i in range(num_pts):
            self.mapper_.insert(cover[i], filtration=-3)
        self.mapper_.initialize_filtration()

        return self
