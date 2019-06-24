"""
@author: Mathieu Carriere
All rights reserved
"""

import numpy as np
import itertools

#from metrics                 import BottleneckDistance
from sklearn.base            import BaseEstimator, TransformerMixin
from sklearn.preprocessing   import LabelEncoder
from sklearn.cluster         import DBSCAN, AgglomerativeClustering
from sklearn.metrics         import pairwise_distances
from scipy.spatial.distance  import directed_hausdorff
from scipy.sparse            import csgraph
from sklearn.neighbors       import KernelDensity, kneighbors_graph, radius_neighbors_graph, NearestNeighbors
import matplotlib.pyplot as plt

try:
    import gudhi as gd
    USE_GUDHI = True

except ImportError:
    USE_GUDHI = False
    print("Gudhi not found: MapperComplex not available")

#############################################
# Clustering ################################
#############################################

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
        if np.any(np.isnan(self.resolutions)):
            delta, resolutions = self.get_optimal_parameters_for_agglomerative_clustering(X=X, beta=0., C=10, N=100)
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
                bottle = gd.bottleneck_distance(D1, D2)
                df = max(df, bottle)
            distribution.append(df)

        return np.sort(distribution)




class DistanceToMeasure(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the distance-to-measure density estimator. 
    """
    def __init__(self, n_neighbors=30, inp="point cloud"):
        """
        Constructor for the DistanceToMeasure class.

        Attributes:
            num_neighbors (int): number of nearest neighbors used to estimate density.
            inp (string): either "point cloud" or "distance matrix". Type of input data (default "point cloud").
        """
        self.input = inp
        self.neighb = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean") if self.input == "point cloud" else NearestNeighbors(n_neighbors=n_neighbors, metric="precomputed")

    def fit(self, X, y=None):
        """
        Fit the DistanceToMeasure class on a point cloud or a distance matrix: compute the nearest neighbors of each point.

        Parameters:
            X (numpy array of shape (num_points) x (num_coordinates) if point cloud and (num_points) x (num_points) if distance matrix): input point cloud or distance matrix.
            y (n x 1 array): point labels (unused).
        """
        self.neighb.fit(X)

    def score_samples(self, X, y=None):
        """
        Compute estimated density values of a new point cloud or distance matrix.
     
        Parameters:
            X (numpy array of shape (num_points) x (num_coordinates) if point cloud and (num_points) x (num_points_indexed) if distance matrix): input point cloud or distance matrix.
            y (n x 1 array): point labels (unused).
        """
        dist, idxs = self.neighb.kneighbors(X) 
        return 1/np.sqrt(np.mean(np.square(dist), axis=1))
        
    
class ToMATo(BaseEstimator, TransformerMixin):
    """
    This is a class for computing ToMATo clustering.
    """
    def __init__(self, n_clusters=None, tau=None, density_estimator=DistanceToMeasure(), n_neighbors=None, radius=None, verbose=False):
        """
        Constructor for the ToMATo class.

        Attributes:
            tau (double): merging parameter (default None). If None, n_clusters is used instead.
            n_clusters (int): number of clusters (default None). If None, it is estimated on data.
            density_estimator (class): density estimator class (default DistanceToMeasure()). Common density estimator classes can be found in the scikit-learn library (such as KernelDensity for instance).
            n_neighbors (int): number of neighbors used to build nearest neighbor graph (default None). If None, delta-neighborhood graph is used.
            radius (double): threshold parameter of delta-neighborhood graph (default None). If None, it is estimated on data.
            verbose (bool): whether to print info during computation (default False).

            labels_ (numpy array of shape (num_points)): clustering labels computed after calling fit() method. 
        """
        self.tau, self.density_estimator, self.n_clusters = tau, density_estimator, n_clusters
        self.n_neighbors, self.radius = n_neighbors, radius
        self.verbose = verbose

    def find(self, i, parents):
        """
        Find function for Union-Find data structure.

        Parameters:
            i (int): ID of point for which parent is required.
            parents (numpy array of shape (num_points)): array storing parents of each point.
        """
        if parents[i] == i:
            return i
        else:
            return self.find(parents[i], parents)

    def union(self, i, j, parents, f):
        """
        Union function for Union-Find data structure. Peak of smaller function value is attached to peak of larger function value.

        Parameters:
            i (int): ID of first point to be merged.
            j (int): ID of second point to be merged.
            parents (numpy array of shape (num_points)): array storing parents of each point.
            f (numpy array of shape (num_points)): array storing function values of each point.
        """
        if f[i] > f[j]:
            parents[j] = i
        else:
            parents[i] = j

    def fit(self, X, y=None):
        """
        Fit the ToMATo class on a point cloud: compute the ToMATo clusters and store the corresponding labels in a numpy array called labels_

        Parameters:
            X (numpy array of shape (num_points) x (num_coordinates)): input point cloud.
            y (n x 1 array): point labels (unused).
        """
        num_pts = X.shape[0]

        if self.verbose:
            print("Computing density estimator")
        self.density_estimator.fit(X)
        self.density_values = self.density_estimator.score_samples(X)
        if self.verbose:
            plt.scatter(X[:,0], X[:,1], s=5., c=self.density_values)
            plt.show()

        if self.verbose:
            print("Computing underlying graph")
        if self.n_neighbors is not None: 
            A = kneighbors_graph(X, self.n_neighbors).toarray()
            A = np.minimum(A + A.T, np.ones(A.shape))
        elif self.radius is not None:
            A = radius_neighbors_graph(X, self.radius).toarray()
        else:
            radius = estimate_scale(X, N=100, inp="point cloud", C=10., beta=0.)
            if self.verbose:
                print("radius = " + str(radius))
            A = radius_neighbors_graph(X, radius).toarray()

        if self.verbose:
            print("Sorting points by density")
        sorted_idxs = np.flip(np.argsort(self.density_values))
        inv_sorted_idxs = np.arange(num_pts)
        for i in range(num_pts):
            inv_sorted_idxs[sorted_idxs[i]] = i

        if self.verbose:
            print("Computing tau")
        if self.tau is not None:
            tau = self.tau
        else:
            st = gd.SimplexTree()
            for i in range(num_pts):
                st.insert([i], filtration=-self.density_values[i])
            for i in range(num_pts):
                for j in range(i+1,num_pts):
                    if A[i,j] == 1.:
                        st.insert([i,j], filtration=max(-self.density_values[i],-self.density_values[j]))
            d = st.persistence()
            plot = gd.plot_persistence_diagram(d)
            plot.show()
            dgm = st.persistence_intervals_in_dimension(0)
            persistences = np.sort([abs(y-x) for (x,y) in dgm])
            if self.n_clusters is not None:
                tau = (persistences[-self.n_clusters-1] + persistences[-self.n_clusters]) / 2
            else:
                n_clusters = np.argmax(np.flip(persistences[1:-1] - persistences[:-2])) + 2
                tau = (persistences[-n_clusters-1] + persistences[-n_clusters]) / 2
        if self.verbose:
            print("tau = " + str(tau))

        if self.verbose:
            print("Applying UF sequentially")
        diag, parents = {}, -np.ones(num_pts, dtype=np.int32)
        for i in range(num_pts):

            current_pt = sorted_idxs[i]
            neighbors = np.squeeze(np.argwhere(A[current_pt,:] == 1.))
            higher_neighbors = [n for n in neighbors if inv_sorted_idxs[n] <= i] if len(neighbors.shape) > 0 else []

            if higher_neighbors == []:

                parents[current_pt] = current_pt
                diag[current_pt] = -np.inf

            else:

                g = higher_neighbors[np.argmax(self.density_values[np.array(higher_neighbors)])]
                pg = self.find(g, parents)
                parents[current_pt] = pg

                for neighbor in higher_neighbors:

                    pn = self.find(neighbor, parents)
                    val = min(self.density_values[pg], self.density_values[pn])

                    if pg != pn and val < self.density_values[current_pt] + tau and val > tau:
                        self.union(pg, pn, parents, self.density_values)
                        pp = pg if self.density_values[pg] < self.density_values[pn] else pn
                        diag[pp] = current_pt

        self.labels_ = np.array([self.find(n, parents) for n in range(num_pts)])
        self.labels_ = LabelEncoder().fit_transform(np.where(self.density_values[self.labels_] > tau, self.labels_, -np.ones(self.labels_.shape)))
