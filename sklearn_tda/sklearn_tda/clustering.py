"""
@author: Mathieu Carriere
All rights reserved
"""

import numpy as np
import itertools

from .metrics                import BottleneckDistance
from sklearn.base            import BaseEstimator, TransformerMixin
from sklearn.cluster         import DBSCAN
from sklearn.metrics         import pairwise_distances
from sklearn.neighbors       import radius_neighbors_graph, kneighbors_graph
from scipy.spatial.distance  import directed_hausdorff
from scipy.sparse            import csgraph

try:
    import gudhi as gd
    USE_GUDHI = True

except ImportError:
    USE_GUDHI = False
    print("Gudhi not found: GraphInducedComplex and MapperComplex not available")

#############################################
# Clustering ################################
#############################################

class NNClustering(BaseEstimator, TransformerMixin):

    def __init__(self, radius, metric="euclidean", inp="point cloud"):
        self.radius, self.metric, self.input = radius, metric, inp

    def fit_predict(self, X):
        if type(self.radius) is int:
            if self.input == "point cloud":
                adj = kneighbors_graph(X, n_neighbors=self.radius, metric=self.metric)
            if self.input == "distance matrix":
                adj = np.zeros(X.shape)
                idxs = np.argpartition(X, self.radius, axis=1)[:, :self.radius]
                for i in range(len(X)):
                    adj[i,idxs[i,:]] = np.ones(len(idxs[i]))                    
        else:
            if self.input == "point cloud":
                adj = radius_neighbors_graph(X, radius=self.radius, metric=self.metric)
            if self.input == "distance matrix":
                adj = np.where(X <= self.radius, np.ones(X.shape), np.zeros(X.shape))
        _, clusters = csgraph.connected_components(adj)
        return clusters


class MapperComplex(BaseEstimator, TransformerMixin):

    def __init__(self, filters=np.array([[0]]), filter_bnds="auto", colors=np.array([[0]]), resolutions=-1, gains=.3, clustering=DBSCAN(), 
                       mask=0, beta=0., C=10, N=100, inp="point cloud", verbose=False):
        self.filters, self.filter_bnds, self.resolutions, self.gains, self.colors, self.clustering = filters, filter_bnds, resolutions, gains, colors, clustering
        self.mask, self.verbose = mask, verbose
        self.input, self.beta, self.C, self.N = inp, beta, C, N

    def get_optimal_parameters_for_hierarchical_clustering(self, X):

        if self.filters.shape[0] == 1 and self.input == "point cloud":
            filters = X[:,self.filters.flatten()]
        else:
            filters = self.filters

        num_pts, num_filt, delta = X.shape[0], filters.shape[1], 0
        m = int(  num_pts / np.exp((1+self.beta) * np.log(np.log(num_pts)/np.log(self.C)))  )
        for _ in range(self.N):
            subpop = np.random.choice(num_pts, size=m, replace=False)
            if self.input == "point cloud":
                d, _, _ = directed_hausdorff(X, X[subpop,:])
            if self.input == "distance matrix":
                d = np.max(np.min(X[:,subpop], axis=1), axis=0)
            delta += d/self.N

        if self.input == "point cloud":
            pairwise = pairwise_distances(X, metric="euclidean")
        if self.input == "distance matrix":
            pairwise = X
        pairs = np.argwhere(pairwise <= delta)
        num_pairs = pairs.shape[0]
        res = []
        for f in range(num_filt):
            F = filters[:,f]
            minf, maxf = np.min(F), np.max(F)
            resf = 0
            for p in range(num_pairs):
                resf = max(resf, abs(F[pairs[p,0]] - F[pairs[p,1]]))
            res.append(int((maxf-minf)/resf))

        return delta, res


    def fit(self, X, y=None):

        if self.filters.shape[0] == 1:
            if self.input == "point cloud":
                filters = X[:, self.filters.flatten()]
            else:
                print("Cannot set filters as coordinates when input is a distance matrix---using eccentricity instead")
                filters = np.max(X, axis=1)[:,np.newaxis]
        else:
            filters = self.filters

        if self.colors.shape[0] == 1:
            if self.input == "point cloud":
                colors = X[:, self.colors.flatten()]
            else:
                print("Cannot set colors as coordinates when input is a distance matrix---using null function instead")
                colors = np.zeros([X.shape[0],1])
        else:
            colors = self.colors

        if isinstance(self.gains, float):
            gains = self.gains * np.ones([filters.shape[1]])
        else:
            gains = self.gains

        if self.resolutions == -1:
            delta, resolutions = self.get_optimal_parameters_for_hierarchical_clustering(X)
            clustering = NNClustering(radius=delta, inp=self.input)
        else:
            resolutions = self.resolutions
            clustering  = self.clustering

        self.st_, self.graph_ = gd.SimplexTree(), []
        self.clus_colors_, self.clus_size_, self.clus_name_, self.clus_subpop_ = dict(), dict(), dict(), dict()

        num_filters, num_colors = filters.shape[1], colors.shape[1]
        interval_inds, intersec_inds = np.empty(filters.shape), np.empty(filters.shape)
        for i in range(num_filters):
            f, r, g = filters[:,i], resolutions[i], gains[i]
            if self.filter_bnds == "auto":
                min_f, max_f = np.min(f), np.max(f)
                epsilon = pow(10, np.log10(abs(max_f)) - 5)
                interval_endpoints, l = np.linspace(min_f - epsilon, max_f + epsilon, num=r+1, retstep=True)
            else:
                min_f, max_f = self.filter_bnds[i,0], self.filter_bnds[i,1]
                interval_endpoints, l = np.linspace(min_f, max_f, num=r+1, retstep=True)
            intersec_endpoints = []
            for j in range(1, len(interval_endpoints)-1):
                intersec_endpoints.append(interval_endpoints[j] - g*l / (2 - 2*g))
                intersec_endpoints.append(interval_endpoints[j] + g*l / (2 - 2*g))
            interval_inds[:,i] = np.digitize(f, interval_endpoints)
            intersec_inds[:,i] = 0.5 * (np.digitize(f, intersec_endpoints) + 1)
            if self.verbose:
                print(interval_inds[:,i])
                print(intersec_inds[:,i])

        num_pts = filters.shape[0]
        binned_data = dict()
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
                if pre_idx in binned_data:
                    binned_data[pre_idx].append(i)
                else:
                    binned_data[pre_idx] = [i]

        if self.verbose:
            print(binned_data)

        cover = []
        for i in range(num_pts):
            cover.append([])

        clus_base = 0
        for preimage in binned_data:

            idxs = np.array(binned_data[preimage])
            if self.input == "point cloud":
                clusters = clustering.fit_predict(X[idxs,:])
            if self.input == "distance matrix":
                clusters = clustering.fit_predict(X[idxs,:][:,idxs])

            if self.verbose:
                print("clusters in preimage " + str(preimage) + " = " + str(clusters))

            num_clus_pre = np.max(clusters) + 1
            for i in range(num_clus_pre):
                subpopulation = idxs[clusters == i]
                color_vals = np.mean(colors[subpopulation,:], axis=0)
                self.clus_colors_[clus_base + i] = color_vals
                self.clus_size_  [clus_base + i] = len(subpopulation)
                self.clus_name_  [clus_base + i] = preimage
                self.clus_subpop_[clus_base + i] = subpopulation

            for i in range(clusters.shape[0]):
                if clusters[i] != -1:
                    cover[idxs[i]].append(clus_base + clusters[i])

            clus_base += np.max(clusters) + 1

        for i in range(num_pts):
            self.st_.insert(cover[i], filtration=-3)

        for simplex in self.st_.get_skeleton(2):
            if len(simplex[0]) > 1:
                idx1, idx2 = simplex[0][0], simplex[0][1]
                if self.mask <= self.clus_size_[idx1] and self.mask <= self.clus_size_[idx2]:
                    self.graph_.append([simplex[0]])
            else:
                clus_idx = simplex[0][0]
                if self.mask <= self.clus_size_[clus_idx]:
                    self.graph_.append([simplex[0], self.clus_colors_[clus_idx], self.clus_size_[clus_idx], self.clus_name_[clus_idx], self.clus_subpop_[clus_idx]])

        return self

    def persistence_diagram(self):
        list_dgm = []
        num_cols = self.clus_colors_[list(self.clus_colors_.keys())[0]].shape[0]
        for c in range(num_cols):
            col_vals = []
            for key, elem in self.clus_colors_.items():
                col_vals.append(elem[c])
            st = gd.SimplexTree()
            list_simplices, list_vertices = self.st_.get_skeleton(1), self.st_.get_skeleton(0)
            for simplex in list_simplices:
                st.insert(simplex[0] + [-2], filtration = -3)
            min_val, max_val = min(col_vals), max(col_vals)
            for vertex in list_vertices:
                if st.find(vertex[0]):
                    st.assign_filtration(vertex[0],        filtration = -2 + (col_vals[vertex[0][0]]-min_val)/(max_val-min_val))
                    st.assign_filtration(vertex[0] + [-2], filtration =  2 - (col_vals[vertex[0][0]]-min_val)/(max_val-min_val))
            st.make_filtration_non_decreasing()
            dgm = st.persistence()
            for point in range(len(dgm)):
                b,d = dgm[point][1][0], dgm[point][1][1]
                b,d = min_val+(2-abs(b))*(max_val-min_val), min_val+(2-abs(d))*(max_val-min_val)
                dgm[point] = tuple([dgm[point][0], tuple([b,d])])
            list_dgm.append(dgm)
        return list_dgm

    def compute_distribution(self, X, N=100):
        num_pts, distribution = len(X), []
        for bootstrap_id in range(N):
            if self.verbose:
                print(str(bootstrap_id) + "th iteration")
            idxs = np.random.choice(num_pts, size=num_pts, replace=True)
            if self.input == "point cloud":
                Xboot = X[idxs,:]
            if self.input == "distance matrix":
                Xboot = X[idxs,:][:,idxs]
            filters_boot = self.filters[idxs,:] if self.filters.shape[0] > 1 else self.filters
            colors_boot  = self.colors[idxs,:]  if self.colors.shape[0]  > 1 else self.colors
            resolutions_boot, gains_boot, clustering_boot = self.resolutions, self.gains, self.clustering 
            Mboot = self.__class__(filters=filters_boot, colors=colors_boot, resolutions=resolutions_boot, gains=gains_boot, clustering=clustering_boot).fit(Xboot)
            dgm1, dgm2 = self.persistence_diagram(), Mboot.persistence_diagram()
            ndg, df = len(dgm1), 0
            for nd in range(ndg):
                npts1, npts2 = len(dgm1[nd]), len(dgm2[nd])
                D1, D2 = [], []
                for pt in range(npts1):
                    if dgm1[nd][pt][0] <= 1:
                        D1.append([dgm1[nd][pt][1][0], dgm1[nd][pt][1][1]])
                for pt in range(npts2):
                    if dgm2[nd][pt][0] <= 1:
                        D2.append([dgm2[nd][pt][1][0], dgm2[nd][pt][1][1]])
                D1, D2 = np.array(D1), np.array(D2)
                bottle = BottleneckDistance().fit([D1])
                df = max(df, bottle.transform([D2])[0][0])
            distribution.append(df)
        return distribution
            
        

class GraphInducedComplex(BaseEstimator, TransformerMixin):

    def __init__(self, graph=-1, graph_subsampling=100, graph_subsampling_power=0.001, graph_subsampling_constant=10,
                       cover_type="functional", filter=0, resolution=-1, gain=0.33, Voronoi_subsampling=1000,
                       mask=0, color=0, verbose=False, input="point cloud"):

        if USE_GUDHI == False:
            raise ImportError("Error: Gudhi not imported")

        self.cc_ = gd.CoverComplex()
        self.cc_.set_type("GIC")
        self.cc_.set_mask(mask)
        self.cc_.set_verbose(verbose)
        self.graph_, self.graph_subsampling_, self.graph_subsampling_constant_, self.graph_subsampling_power_ = graph, graph_subsampling, graph_subsampling_constant, graph_subsampling_power
        self.cover_type_, self.filter_, self.resolution_, self.gain_, self.Voronoi_subsampling_ = cover_type, filter, resolution, gain, Voronoi_subsampling
        self.color_, self.input = color, input

    def fit(self, X, y=None):

        # Read input
        if self.input == "point cloud":
            self.cc_.set_point_cloud_from_range(X)
        elif self.input == "distance matrix":
            self.cc_.set_distances_from_range(X)

        # Set color function
        if type(self.color_) is int:
            self.cc_.set_color_from_coordinate(self.color_)
        if type(self.color_) is np.ndarray:
            self.cc_.set_color_from_range(self.color_)

        # Set underlying neighborhood graph for connected components
        if self.graph_ == -1:
            self.cc_.set_subsampling(self.graph_subsampling_constant_, self.graph_subsampling_power_)
            self.cc_.set_graph_from_automatic_rips(self.graph_subsampling_)
        else:
            self.cc_.set_graph_from_rips(self.graph_)

        # Set cover of point cloud
        if self.cover_type_ == "functional":
            ###### Function values
            if type(self.filter_) is int:
                self.cc_.set_function_from_coordinate(self.filter_)
            if type(self.filter_) is np.ndarray:
                self.cc_.set_function_from_range(self.filter_)
            ###### Gain
            self.cc_.set_gain(self.gain_)
            ###### Resolution
            if self.resolution_ == -1:
                self.cc_.set_automatic_resolution()
            else:
                if type(self.resolution_) is int:
                    self.cc_.set_resolution_with_interval_number(self.resolution_)
                else:
                    self.cc_.set_resolution_with_interval_length(self.resolution_)
            ###### Cover computation
            self.cc_.set_cover_from_function()
        if self.cover_type_ == "Voronoi":
            self.cc_.set_cover_from_Voronoi(self.Voronoi_subsampling_)

        # Compute simplices
        self.cc_.find_simplices()
        self.cc_.create_simplex_tree()

        return self

    def print_result(self, output_type="txt"):
        if output_type == "txt":
            self.cc_.write_info()
        if output_type == "dot":
            self.cc_.plot_dot()
        if output_type == "off":
            self.cc_.plot_off()

    def compute_p_value(self, bootstrap=10):
        self.cc_.compute_distribution(bootstrap)
        return self.cc_.compute_p_value()

    def compute_confidence_level_from_distance(self, bootstrap=10, distance=1.0):
        self.cc_.compute_distribution(bootstrap)
        return self.cc_.compute_confidence_level_from_distance(distance)

    def compute_distance_from_confidence_level(self, bootstrap=10, alpha=0.1):
        self.cc_.compute_distribution(bootstrap)
        return self.cc_.compute_distance_from_confidence_level(alpha)

    def subpopulation(self, node_index=0):
        return self.cc_.subpopulation(node_index)

    def persistence_diagram(self):
        return self.cc_.compute_PD()
