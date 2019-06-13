"""
@author: Mathieu Carriere
All rights reserved
"""

import numpy as np
import itertools

from metrics                 import BottleneckDistance
from sklearn.base            import BaseEstimator, TransformerMixin
from sklearn.cluster         import DBSCAN, AgglomerativeClustering
from sklearn.metrics         import pairwise_distances
from sklearn.neighbors       import radius_neighbors_graph, kneighbors_graph
from scipy.spatial.distance  import directed_hausdorff
from scipy.sparse            import csgraph

try:
    import gudhi as gd
    USE_GUDHI = True

except ImportError:
    USE_GUDHI = False
    print("Gudhi not found: MapperComplex not available")

#############################################
# Clustering ################################
#############################################

class MapperComplex(BaseEstimator, TransformerMixin):
    """
    This is a class for computing Mapper simplicial complexes on point clouds or distance matrices. 
    """
    def __init__(self, filters, filter_bnds, colors, resolutions, gains, inp="point cloud", clustering=DBSCAN(), mask=0):
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

            mapper_ (gudhi SimplexTree): Mapper simplicial complex computed after calling the fit() method
            node_info_ (dictionary): various information associated to the nodes of the Mapper. 
        """
        self.filters, self.filter_bnds, self.resolutions, self.gains, self.colors, self.clustering = filters, filter_bnds, resolutions, gains, colors, clustering
        self.input, self.mask = inp, mask

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
        m = int(  num_pts / np.exp((1+beta) * np.log(np.log(num_pts)/np.log(C)))  )
        for _ in range(N):
            subpop = np.random.choice(num_pts, size=m, replace=False)
            if self.input == "point cloud":
                d, _, _ = directed_hausdorff(X, X[subpop,:])
            if self.input == "distance matrix":
                d = np.max(np.min(X[:,subpop], axis=1), axis=0)
            delta += d/N

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
        if np.any(np.isnan(self.resolutions)):
            delta, resolutions = self.get_optimal_parameters_for_agglomerative_clustering(X=X, beta=0., C=10, N=100)
            #self.clustering = NNClustering(radius=delta, inp=self.input)  
            if self.input == "point cloud":
                self.clustering = AgglomerativeClustering(n_clusters=None, linkage="single", distance_threshold=delta, affinity="euclidean")  
            else:
                self.clustering = AgglomerativeClustering(n_clusters=None, linkage="single", distance_threshold=delta, affinity="precomputed")
            self.resolutions = np.where(np.isnan(self.resolutions), resolutions, self.resolutions)

        # If some filter limits are unspecified, automatically compute them
        self.filter_bnds = np.where(np.isnan(self.filter_bnds), np.hstack([np.min(self.filters, axis=0)[:,np.newaxis], np.max(self.filters, axis=0)[:,np.newaxis]]), self.filter_bnds)

        # Initialize attributes
        self.mapper_, self.node_info_ = gd.SimplexTree(), {}

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
                if len(subpopulation) >= self.mask:
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

    def compute_persistence_diagrams(self):
        """
        Compute the extended persistence diagrams of the Mapper simplicial complex associated to each color function.

        Returns:
            list_dgm (list of gudhi persistence diagrams): output extended persistence diagrams. There is one per color function.
        """
        num_cols, list_dgm = self.colors.shape[1], []

        # Compute an extended persistence diagram for each color
        for c in range(num_cols):

            # Retrieve all color values
            col_vals = {node_name: self.node_info_[node_name]["colors"][c] for node_name in self.node_info_.keys()}
            
            # Create a new simplicial complex by coning the Mapper with an extra point with name -2
            st = gd.SimplexTree()
            list_simplices, list_vertices = self.mapper_.get_skeleton(1), self.mapper_.get_skeleton(0)
            for (simplex, f) in list_simplices:
                st.insert(simplex + [-2], filtration=-3)

            # Assign ascending filtration values on the original simplices and descending filtration values on the coned simplices 
            min_val, max_val = min(col_vals), max(col_vals)
            for (vertex, f) in list_vertices:
                if st.find(vertex):
                    st.assign_filtration(vertex,        filtration = -2 + (col_vals[vertex[0]]-min_val)/(max_val-min_val))
                    st.assign_filtration(vertex + [-2], filtration =  2 - (col_vals[vertex[0]]-min_val)/(max_val-min_val))

            # Compute persistence
            st.make_filtration_non_decreasing()
            dgm = st.persistence()

            # Output extended persistence diagrams
            for point in range(len(dgm)):
                b,d = dgm[point][1][0], dgm[point][1][1]
                b,d = min_val+(2-abs(b))*(max_val-min_val), min_val+(2-abs(d))*(max_val-min_val)
                dgm[point] = tuple([dgm[point][0], tuple([b,d])])
            list_dgm.append(dgm)

        return list_dgm

    def compute_distribution(self, X, N=100):
        """
        Compute a bootstrap distribution of bottleneck distances. More specifically, subsample the input point cloud or distance matrix, compute the Mapper with the same parameters on this subsample, and compare its extended persistence diagrams with the original ones.

        Parameters:
            X (numpy array of shape (num_points) x (num_coordinates) if point cloud and (num_points) x (num_points) if distance matrix): input point cloud or distance matrix.
            N (int): bootstrap iterations (default 100).

        Returns:
            distribution: list of bottleneck distance values.
        """
        num_pts, distribution = len(X), []
        dgm = self.compute_persistence_diagrams()

        for bootstrap_id in range(N):

            print(str(bootstrap_id) + "th iteration")

            # Randomly select points
            idxs = np.random.choice(num_pts, size=num_pts, replace=True)
            Xboot = X[idxs,:] if self.input == "point cloud" else X[idxs,:][:,idxs]
            f_boot, c_boot = self.filters[idxs,:], self.colors[idxs,:]
            Mboot = self.__class__(filters=f_boot, filter_bnds=self.filter_bnds, colors=c_boot, resolutions=self.resolutions, gains=self.gains, inp=self.input, clustering=self.clustering).fit(Xboot)

            # Compute the corresponding persistence diagrams
            dgm_boot = Mboot.compute_persistence_diagrams()

            # Compute the bottleneck distances between them and keep the maximum
            df = 0.
            for i in range(len(dgm)):
                npts, npts_boot = len(dgm[i]), len(dgm_boot[i])
                D1 = np.array([[dgm[i][pt][1][0], dgm[i][pt][1][1]] for pt in range(npts) if dgm[i][pt][0] <= 1]) 
                D2 = np.array([[dgm_boot[i][pt][1][0], dgm_boot[i][pt][1][1]] for pt in range(npts_boot) if dgm_boot[i][pt][0] <= 1])
                bottle = BottleneckDistance().fit([D1])
                df = max(df, float(np.squeeze(bottle.transform([D2]))))
            distribution.append(df)

        return np.sort(distribution)
