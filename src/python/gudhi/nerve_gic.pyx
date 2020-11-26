# This file is part of the Gudhi Library - https://gudhi.inria.fr/ -
# which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full
# license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2018 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from __future__ import print_function
from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp cimport bool
import errno
import os
from libc.stdint cimport intptr_t

from gudhi.simplex_tree cimport *
from gudhi.simplex_tree import SimplexTree

import numpy as np
import itertools
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import colorsys

from networkx                import cycle_basis
from scipy.sparse.csgraph    import dijkstra, shortest_path, connected_components
from scipy.sparse            import csr_matrix
from sklearn.base            import BaseEstimator, TransformerMixin
from sklearn.cluster         import DBSCAN, AgglomerativeClustering
from sklearn.metrics         import pairwise_distances
from scipy.spatial.distance  import directed_hausdorff
from scipy.stats             import ks_2samp

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2018 Inria"
__license__ = "GPL v3"

cdef extern from "Nerve_gic_interface.h" namespace "Gudhi":
    cdef cppclass Nerve_gic_interface "Gudhi::cover_complex::Nerve_gic_interface":
        Nerve_gic_interface()
        double compute_confidence_level_from_distance(double distance)
        double compute_distance_from_confidence_level(double alpha)
        void compute_distribution(int N)
        double compute_p_value()
        vector[pair[double, double]] compute_PD()
        void find_simplices()
        void create_simplex_tree(Simplex_tree_interface_full_featured* simplex_tree)
        bool read_point_cloud(string off_file_name)
        double set_automatic_resolution()
        void set_color_from_coordinate(int k)
        void set_color_from_file(string color_file_name)
        void set_color_from_range(vector[double] color)
        void set_cover_from_file(string cover_file_name)
        void set_cover_from_function()
        void set_cover_from_range(vector[vector[int]] assignments)
        void set_cover_from_Euclidean_Voronoi(int m)
        void set_function_from_coordinate(int k)
        void set_function_from_file(string func_file_name)
        void set_function_from_range(vector[double] function)
        void set_gain(double g)
        double set_graph_from_automatic_euclidean_rips(int N)
        void set_graph_from_file(string graph_file_name)
        void set_graph_from_OFF()
        void set_graph_from_euclidean_rips(double threshold)
        void set_mask(int nodemask)
        void set_resolution_with_interval_length(double resolution)
        void set_resolution_with_interval_number(int resolution)
        void set_subsampling(double constant, double power)
        void set_type(string type)
        void set_verbose(bool verbose)
        vector[int] subpopulation(int c)
        double subcolor(int c)
        void write_info(string data_name, string cover_name, string color_name)
        void plot_DOT(string data_name)
        void plot_OFF(string data_name)
        void set_point_cloud_from_range(vector[vector[double]] cloud)
        void set_distances_from_range(vector[vector[double]] distance_matrix)

# NGIComplex python interface
cdef class NGIComplex:
    """Cover complex data structure for Nerve and Graph Induced complexes.

    The data structure is a simplicial complex, representing a Graph Induced
    simplicial Complex (GIC) or a Nerve, and whose simplices are computed with
    a cover C of a point cloud P. To compute a GIC, one also needs a
    graph G built on top of P, whose cliques with vertices belonging to
    different elements of C correspond to the simplices of the GIC.
    """

    cdef Nerve_gic_interface * thisptr

    # Fake constructor that does nothing but documenting the constructor
    def __init__(self):
        """NGIComplex constructor.
        """

    # The real cython constructor
    def __cinit__(self):
        self.thisptr = new Nerve_gic_interface()

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr

    def __is_defined(self):
        """Returns true if NGIComplex pointer is not NULL.
         """
        return self.thisptr != NULL

    def set_point_cloud_from_range(self, cloud):
        """ Reads and stores the input point cloud from a vector stored in
        memory.

        :param cloud: Input vector containing the point cloud.
        :type cloud: vector[vector[double]]
        """
        return self.thisptr.set_point_cloud_from_range(cloud)

    def set_distances_from_range(self, distance_matrix):
        """ Reads and stores the input distance matrix from a vector stored in
        memory.

        :param distance_matrix: Input vector containing the distance matrix.
        :type distance_matrix: vector[vector[double]]
        """
        return self.thisptr.set_distances_from_range(distance_matrix)

    def compute_confidence_level_from_distance(self, distance):
        """Computes the confidence level of a specific bottleneck distance
        threshold.

        :param distance: Bottleneck distance.
        :type distance: double
        :rtype: double
        :returns: Confidence level.
        """
        return self.thisptr.compute_confidence_level_from_distance(distance)

    def compute_distance_from_confidence_level(self, alpha):
        """Computes the bottleneck distance threshold corresponding to a
        specific confidence level.

        :param alpha: Confidence level.
        :type alpha: double
        :rtype: double
        :returns: Bottleneck distance.
        """
        return self.thisptr.compute_distance_from_confidence_level(alpha)

    def compute_distribution(self, N=100):
        """Computes bootstrapped distances distribution.

        :param N: Loop number (default value is 100).
        :type alpha: int
        """
        self.thisptr.compute_distribution(N)

    def compute_p_value(self):
        """Computes the p-value, i.e. the opposite of the confidence level of
        the largest bottleneck distance preserving the points in the
        persistence diagram of the output simplicial complex.

        :rtype: double
        :returns: p-value.
        """
        return self.thisptr.compute_p_value()

    def compute_PD(self):
        """Computes the extended persistence diagram of the complex.
        """
        return self.thisptr.compute_PD()

    def create_simplex_tree(self):
        """
        :returns: A simplex tree created from the Cover complex.
        :rtype: SimplexTree
        """
        stree = SimplexTree()
        cdef intptr_t stree_int_ptr=stree.thisptr
        self.thisptr.create_simplex_tree(
            <Simplex_tree_interface_full_featured*>stree_int_ptr)
        return stree

    def find_simplices(self):
        """Computes the simplices of the simplicial complex.
        """
        self.thisptr.find_simplices()

    def read_point_cloud(self, off_file):
        """Reads and stores the input point cloud from .(n)OFF file.

        :param off_file: Name of the input .OFF or .nOFF file.
        :type off_file: string
        :rtype: bool
        :returns: Read file status.
        """
        if os.path.isfile(off_file):
            return self.thisptr.read_point_cloud(off_file.encode('utf-8'))
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                    off_file)

    def set_automatic_resolution(self):
        """Computes the optimal length of intervals (i.e. the smallest interval
        length avoiding discretization artifacts - see :cite:`Carriere17c`) for a
        functional cover.

        :rtype: double
        :returns: reso interval length used to compute the cover.
        """
        return self.thisptr.set_automatic_resolution()

    def set_color_from_coordinate(self, k=0):
        """Computes the function used to color the nodes of the simplicial
        complex from the k-th coordinate.

        :param k: Coordinate to use (start at 0). Default value is 0.
        :type k: int
        """
        return self.thisptr.set_color_from_coordinate(k)

    def set_color_from_file(self, color_file_name):
        """Computes the function used to color the nodes of the simplicial
        complex from a file containing the function values.

        :param color_file_name: Name of the input color file.
        :type color_file_name: string
        """
        if os.path.isfile(color_file_name):
            self.thisptr.set_color_from_file(color_file_name.encode('utf-8'))
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                    color_file_name)

    def set_color_from_range(self, color):
        """Computes the function used to color the nodes of the simplicial
        complex from a vector stored in memory.

        :param color: Input vector of values.
        :type color: vector[double]
        """
        self.thisptr.set_color_from_range(color)

    def set_cover_from_file(self, cover_file_name):
        """Creates the cover C from a file containing the cover elements of
        each point (the order has to be the same as in the input file!).

        :param cover_file_name: Name of the input cover file.
        :type cover_file_name: string
        """
        if os.path.isfile(cover_file_name):
            self.thisptr.set_cover_from_file(cover_file_name.encode('utf-8'))
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                    cover_file_name)

    def set_cover_from_function(self):
        """Creates a cover C from the preimages of the function f.
        """
        self.thisptr.set_cover_from_function()

    def set_cover_from_range(self, assignments):
        """Creates a cover C from a vector stored in memory.
        """
        self.thisptr.set_cover_from_range(assignments)

    def set_cover_from_Voronoi(self, m=100):
        """Creates the cover C from the Voronoï cells of a subsampling of the
        point cloud.

        :param m: Number of points in the subsample. Default value is 100.
        :type m: int
        """
        self.thisptr.set_cover_from_Euclidean_Voronoi(m)

    def set_function_from_coordinate(self, k):
        """Creates the function f from the k-th coordinate of the point cloud.

        :param k: Coordinate to use (start at 0).
        :type k: int
        """
        self.thisptr.set_function_from_coordinate(k)

    def set_function_from_file(self, func_file_name):
        """Creates the function f from a file containing the function values.

        :param func_file_name: Name of the input function file.
        :type func_file_name: string
        """
        if os.path.isfile(func_file_name):
            self.thisptr.set_function_from_file(func_file_name.encode('utf-8'))
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                    func_file_name)

    def set_function_from_range(self, function):
        """Creates the function f from a vector stored in memory.

        :param function: Input vector of values.
        :type function: vector[double]
        """
        self.thisptr.set_function_from_range(function)

    def set_gain(self, g = 0.3):
        """Sets a gain from a value stored in memory.

        :param g: Gain (default value is 0.3).
        :type g: double
        """
        self.thisptr.set_gain(g)

    def set_graph_from_automatic_rips(self, N=100):
        """Creates a graph G from a Rips complex whose threshold value is
        automatically tuned with subsampling - see :cite:`Carriere17c`.

        :param N: Number of subsampling iteration (the default reasonable value
            is 100, but there is no guarantee on how to choose it).
        :type N: int
        :rtype: double
        :returns: Delta threshold used for computing the Rips complex.
        """
        return self.thisptr.set_graph_from_automatic_euclidean_rips(N)

    def set_graph_from_file(self, graph_file_name):
        """Creates a graph G from a file containing the edges.

        :param graph_file_name: Name of the input graph file. The graph file
            contains one edge per line, each edge being represented by the IDs
            of its two nodes.
        :type graph_file_name: string
        """
        if os.path.isfile(graph_file_name):
            self.thisptr.set_graph_from_file(graph_file_name.encode('utf-8'))
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                    graph_file_name)

    def set_graph_from_OFF(self):
        """Creates a graph G from the triangulation given by the input OFF
        file.
        """
        self.thisptr.set_graph_from_OFF()

    def set_graph_from_rips(self, threshold):
        """Creates a graph G from a Rips complex.

        :param threshold: Threshold value for the Rips complex.
        :type threshold: double
        """
        self.thisptr.set_graph_from_euclidean_rips(threshold)

    def set_mask(self, nodemask):
        """Sets the mask, which is a threshold integer such that nodes in the
        complex that contain a number of data points which is less than or
        equal to this threshold are not displayed.

        :param nodemask: Threshold.
        :type nodemask: int
        """
        self.thisptr.set_mask(nodemask)

    def set_resolution_with_interval_length(self, resolution):
        """Sets a length of intervals from a value stored in memory.

        :param resolution: Length of intervals.
        :type resolution: double
        """
        self.thisptr.set_resolution_with_interval_length(resolution)

    def set_resolution_with_interval_number(self, resolution):
        """Sets a number of intervals from a value stored in memory.

        :param resolution: Number of intervals.
        :type resolution: int
        """
        self.thisptr.set_resolution_with_interval_number(resolution)

    def set_subsampling(self, constant, power):
        """Sets the constants used to subsample the data set. These constants
        are explained in :cite:`Carriere17c`.

        :param constant: Constant.
        :type constant: double
        :param power: Power.
        :type resolution: double
        """
        self.thisptr.set_subsampling(constant, power)

    def set_type(self, type):	
        """Specifies whether the type of the output simplicial complex.	
        :param type: either "GIC" or "Nerve".	
        :type type: string	
        """	
        self.thisptr.set_type(type.encode('utf-8'))

    def set_verbose(self, verbose):
        """Specifies whether the program should display information or not.

        :param verbose: true = display info, false = do not display info.
        :type verbose: boolean
        """
        self.thisptr.set_verbose(verbose)

    def subpopulation(self, c):
        """Returns the data subset corresponding to a specific node of the
        created complex.

        :param c: ID of the node.
        :type c: int
        :rtype: vector[int]
        :returns: Vector of IDs of data points.
        """
        return self.thisptr.subpopulation(c)

    def subcolor(self, c):
        """Returns the mean color value corresponding to a specific node of the
        created complex.

        :param c: ID of the node.
        :type c: int
        :rtype: double
        :returns: Mean color value of data points.
        """
        return self.thisptr.subcolor(c)

    def write_info(self, dataname, covername, colorname):
        """Creates a .txt file called SC.txt describing the 1-skeleton, which can
        then be plotted with e.g. KeplerMapper.
        """
        return self.thisptr.write_info(dataname, covername, colorname)

    def plot_dot(self, name):
        """Creates a .dot file called SC.dot for neato (part of the graphviz
        package) once the simplicial complex is computed to get a visualization of
        its 1-skeleton in a .pdf file.
        """
        return self.thisptr.plot_DOT(name)

    def plot_off(self, name):
        """Creates a .off file called SC.off for 3D visualization, which contains
        the 2-skeleton of the GIC. This function assumes that the cover has been
        computed with Voronoi. If data points are in 1D or 2D, the remaining
        coordinates of the points embedded in 3D are set to 0.
        """
        return self.thisptr.plot_OFF(name)

















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


class _MapperComplex(BaseEstimator, TransformerMixin):
    """
    This is a class for computing Mapper simplicial complexes on point clouds or distance matrices. Used internally.
    """
    def __init__(self, filters, filter_bnds, colors, resolutions, gains, inp="point cloud", clustering=DBSCAN(), mask=0, N=100, beta=0., C=10.):
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

        return self


class CoverComplex(BaseEstimator, TransformerMixin):
    """
    Constructor for the CoverComplex class. This class wraps Mapper, Nerve and Graph Induced complexes in a single interface. Graph Induced and Nerve complexes can still be called from the class NGIComplex (with a few more functionalities, such as defining datasets or graphs through files). Key differences between Mapper, Nerve and Graph Induced complexes (GIC) are: Mapper nodes are defined with given input clustering method while GIC nodes are defined with given input graph and Nerve nodes are defined with cover elements, GIC accepts partitions instead of covers while Mapper and Nerve require cover elements to overlap. Also, note that when the cover is functional (i.e., preimages of filter functions), GIC only accepts one scalar-valued filter with gain < 0.5, meaning that the arguments "resolutions" and "gains" should have length 1. If you have more than one scalar filter, or if the gain is more than 0.5, the cover should be computed beforehand and fed to the class with the "assignments" argument. On the other hand, Mapper and Nerve complexes accept "resolutions" and "gains" with any length. 

    Attributes:
        complex_type (string): type of cover complex. Either "mapper", "gic" or "nerve". 
        input_type (string): type of input data. Either "point cloud" or "distance matrix".
        cover (string): specifies the cover. Either "functional" (preimages of filter function), "voronoi" or "precomputed".
        colors (numpy array of shape (num_points) x (num_colors)): functions used to color the nodes of the cover complex. More specifically, coloring is done by computing the means of these functions on the subpopulations corresponding to each node. If None, first coordinate is used if input is point cloud, and eccentricity is used if input is distance matrix.
        mask (int): threshold on the size of the cover complex nodes (default 0). Any node associated to a subpopulation with less than **mask** points will be removed.
        voronoi_samples (int): number of Voronoi germs used for partitioning the input dataset. Used only if complex_type = "gic" and cover = "voronoi".
        assignments (list of length (num_points) of lists of integers): cover assignment for each point. Used only if complex_type = "gic" or "nerve" and cover = "precomputed".
        filters (numpy array of shape (num_points) x (num_filters)): filter functions (sometimes called lenses) used to compute the cover. Each column of the numpy array defines a scalar function defined on the input points. Used only if cover = "functional".
        filter_bnds (numpy array of shape (num_filters) x 2): limits of each filter, of the form [[f_1^min, f_1^max], ..., [f_n^min, f_n^max]]. If one of the values is numpy.nan, it can be computed from the dataset with the fit() method. Used only if cover = "functional".
        resolutions (numpy array of shape num_filters containing integers): resolution of each filter function, ie number of intervals required to cover each filter image. Must be of length 1 if complex_type = "gic". Used only if cover = "functional". If None, it is estimated from data.
        gains (numpy array of shape num_filters containing doubles in [0,1]): gain of each filter function, ie overlap percentage of the intervals covering each filter image. Must be of length 1 if complex_type = "gic". Used only if cover = "functional".
        N (int): subsampling iterations (default 100) for estimating scale and resolutions. Used only if cover = "functional" and clustering or resolutions = None. See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
        beta (double): exponent parameter (default 0.) for estimating scale and resolutions. Used only if cover = "functional" and clustering or resolutions = None. See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
        C (double): constant parameter (default 10.) for estimating scale and resolutions. Used only if cover = "functional" and clustering or resolutions = None. See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
        clustering (class): clustering class (default sklearn.cluster.DBSCAN()). Common clustering classes can be found in the scikit-learn library (such as AgglomerativeClustering for instance). Used only if complex_type = "mapper". If None, it is set to hierarchical clustering, with scale estimated from data.
        graph (string): type of graph to use for GIC. Used only if complex_type = "gic". Currently accepts "rips" only.
        rips_threshold (float): Rips parameter. Used only if complex_type = "gic" and graph = "rips".
        input_name (string): name of dataset. Used when generating plots.
        cover_name (string): name of cover. Used when generating plots.
        color_name (string): name of color function. Used when generating plots.
        verbose (bool): whether to display info while computing.
 
        simplex_tree (gudhi SimplexTree): simplicial complex representing the cover complex computed after calling the fit() method.
        node_info_ (dictionary): various information associated to the nodes of the cover complex. 
    """
    def __init__(self, complex_type="mapper", input_type="point cloud", cover="functional", colors=None, mask=0,
                       voronoi_samples=100, assignments=None, filters=None, filter_bnds=None, resolutions=None, gains=None, N=100, beta=0., C=10.,
                       clustering=None,
                       graph="rips", rips_threshold=None, 
                       input_name="data", cover_name="cover", color_name="color", verbose=False):

        self.complex_type, self.input_type, self.cover, self.colors, self.mask = complex_type, input_type, cover, colors, mask
        self.voronoi_samples, self.assignments, self.filters, self.filter_bnds, self.resolutions, self.gains, self.clustering = voronoi_samples, assignments, filters, filter_bnds, resolutions, gains, clustering
        self.graph, self.rips_threshold, self.N, self.beta, self.C = graph, rips_threshold, N, beta, C
        self.input_name, self.cover_name, self.color_name, self.verbose = input_name, cover_name, color_name, verbose

    def fit(self, X, y=None):
        """
        Fit the CoverComplex class on a point cloud or a distance matrix: compute the cover complex and store it in a simplex tree called simplex_tree.

        Parameters:
            X (numpy array of shape (num_points) x (num_coordinates) if point cloud and (num_points) x (num_points) if distance matrix): input point cloud or distance matrix.
            y (n x 1 array): point labels (unused).
        """
        self.data = X

        if self.colors is None:
            if self.input_type == "point cloud":
                self.colors = X[:,0] if self.complex_type == "gic" else X[:,0:1]
            elif self.input_type == "distance matrix":
                self.colors = X.max(axis=0) if self.complex_type == "gic" else X.max(axis=0)[:,np.newaxis]

        if self.filters is None:
            if self.input_type == "point cloud":
                self.filters = X[:,0] if self.complex_type == "gic" else X[:,0:1]
            elif self.input_type == "distance matrix":
                self.filters = X.max(axis=0) if self.complex_type == "gic" else X.max(axis=0)[:,np.newaxis]

        if self.complex_type == "gic" or self.complex_type == "nerve":

            if self.complex_type is not "nerve" or self.cover is not "functional": # Nerve + functional cover is a special case where it's better to use Mapper than C++ code

                ct = "GIC" if self.complex_type == "gic" else "Nerve"
                self.complex = NGIComplex()
                self.complex.set_type(ct)
                self.complex.set_verbose(self.verbose)

                if self.input_type == "point cloud":
                    self.complex.set_point_cloud_from_range(X)
                elif self.input_type == "distance matrix":
                    self.complex.set_distances_from_range(X)
            
                self.complex.set_color_from_range(self.colors)
            
                if self.complex_type == "gic":
                    if self.graph == "rips":
                        if self.rips_threshold is not None:
                            self.complex.set_graph_from_rips(self.rips_threshold)
                        else:
                            self.complex.set_subsampling(self.C, self.beta)
                            self.complex.set_graph_from_automatic_rips(self.N)

            if self.cover == "voronoi":
                assert self.complex_type is not "nerve"
                self.complex.set_cover_from_Voronoi(self.voronoi_samples)

            elif self.cover == "functional":

                if self.complex_type == "gic":

                    self.complex.set_function_from_range(self.filters)

                    if self.resolutions is None:
                        self.complex.set_automatic_resolution()
                    else:
                        self.complex.set_resolution_with_interval_number(self.resolutions[0])

                    if self.gains is None:
                        self.complex.set_gain(.33)
                    else:
                        self.complex.set_gain(self.gains[0])

                    self.complex.set_cover_from_function()

                elif self.complex_type == "nerve": # it's actually better to use Mapper with constant clustering
                    self.complex = _MapperComplex(filters=self.filters, filter_bnds=self.filter_bnds, colors=self.colors, 
                                                 resolutions=self.resolutions, gains=self.gains, inp=self.input_type, 
                                                 clustering=self._constant_clustering, mask=self.mask, N=self.N, beta=self.beta, C=self.C)
                    self.complex.fit(X)
                    self.simplex_tree = self.complex.mapper_
                    self.node_info = self.complex.node_info_

            elif self.cover == "precomputed":
                self.complex.set_cover_from_range(self.assignments)

            if self.complex_type is not "nerve" or self.cover is not "functional":

                self.complex.set_mask(self.mask)

                self.complex.find_simplices()
                simplex_tree = self.complex.create_simplex_tree()
            
                self.simplex_tree = SimplexTree()
                idv, names = 0, {}
                for v,_ in simplex_tree.get_skeleton(0):
                    if len(self.complex.subpopulation(v[0])) > self.mask:
                        names[v[0]] = idv
                        self.simplex_tree.insert([idv])
                        idv += 1
                for s,_ in simplex_tree.get_filtration():
                    if len(s) >= 2 and np.all([len(self.complex.subpopulation(v)) > self.mask for v in s]):
                        self.simplex_tree.insert([names[v] for v in s])
                self.node_info = {}
                for v,_ in simplex_tree.get_skeleton(0):
                    if len(self.complex.subpopulation(v[0])) > self.mask:
                        node = names[v[0]]
                        pop = self.complex.subpopulation(v[0])
                        self.node_info[node] = {"indices": pop, "size": len(pop), "colors": [self.complex.subcolor(v[0])]}
    

        elif self.complex_type == "mapper":

            assert self.cover is not "voronoi"
            self.complex = _MapperComplex(filters=self.filters, filter_bnds=self.filter_bnds, colors=self.colors, 
                                         resolutions=self.resolutions, gains=self.gains, inp=self.input_type, 
                                         clustering=self.clustering, mask=self.mask, N=self.N, beta=self.beta, C=self.C)
            self.complex.fit(X)
            self.simplex_tree = self.complex.mapper_
            self.node_info = self.complex.node_info_

        return self

    def get_networkx(self, get_attrs=False):
        """
        Turn the 1-skeleton of the cover complex computed after calling fit() method into a networkx graph.

        Parameters:
            get_attrs (bool): if True, the color functions will be used as attributes for the networkx graph.

        Returns:
            G (networkx graph): graph representing the 1-skeleton of the cover complex.
        """
        st = self.simplex_tree
        G = nx.Graph()
        for (splx,_) in st.get_skeleton(1):	
            if len(splx) == 1:
                G.add_node(splx[0])
            if len(splx) == 2:
                G.add_edge(splx[0], splx[1])
        if get_attrs:
            attrs = {k: {"attr_name": self.node_info[k]["colors"]} for k in G.nodes()}
            nx.set_node_attributes(G, attrs)
        return G

    class _constant_clustering():
        def fit_predict(X):
            return np.zeros([len(X)], dtype=np.int32)

    def compute_topological_features(self, threshold=0.):
        """
        Compute the topological features (connected components, up/down branches, loops) of the 1-skeleton of the cover complex. Connected components and loops are computed with scipy functions, and branches are detected with Union-Find and 0-dimensional persistence of the 1-skeleton.

        Parameters:
            threshold (float): any topological feature whose size is less than this parameter (relative to the first color function) will be discarded.

        Returns:
            dgm (list of (dim,(a,b)) tuples): list of feature characteristics. dim is the topological dimension of the feature (0 for CCs and branches, 1 for loops), a,b are the min and max of the first color function along the feature.
            bnds (list of lists): list of feature points. Each element of this list is the list of point IDs forming the corresponding feature. 
        """
        st = self.simplex_tree
        function = [self.node_info[v[0]]["colors"][0] for v,_ in st.get_skeleton(0)]
        num_nodes = st.num_vertices()
        dgm, bnd = [], []

        # connected_components
        A = np.zeros([num_nodes, num_nodes])
        for (splx,_) in st.get_skeleton(1):
            if len(splx) == 2:	
                A[splx[0], splx[1]] = 1
                A[splx[1], splx[0]] = 1
        _, ccs = connected_components(A, directed=False)
        for ccID in np.unique(ccs):
            pts = np.argwhere(ccs == ccID).flatten()
            vals = [function[p] for p in pts]
            if np.abs(min(vals) - max(vals)) >= threshold:
                dgm.append((0, (min(vals), max(vals))))
                bnd.append(pts)

        # loops
        G = self.get_networkx()
        bndall = cycle_basis(G)
        for pts in bndall:
            vals = [function[p] for p in pts]
            if np.abs(min(vals) - max(vals)) >= threshold:	
                dgm.append((1,(min(vals), max(vals))))
                bnd.append(pts)
        
        # branches
        for topo_type in ["downbranch", "upbranch"]:

            # upranch is downbranch of opposite function
            if topo_type == "upbranch":
                function = [-f for f in function]

            # sort vertices according to function values and compute inverse function 
            sorted_idxs = np.argsort(np.array(function))
            inv_sorted_idxs = np.zeros(num_nodes)
            for i in range(num_nodes):
                inv_sorted_idxs[sorted_idxs[i]] = i

            # go through all vertices in ascending function order
            persistence_diag, persistence_set, parents, visited = {}, {}, -np.ones(num_nodes, dtype=np.int32), {}
            for i in range(num_nodes):

                current_pt = sorted_idxs[i]
                neighbors = np.ravel(np.argwhere(A[current_pt,:] == 1))
                lower_neighbors = [n for n in neighbors if inv_sorted_idxs[n] <= i] if len(neighbors) > 0 else []

                # no lower neighbors: current point is a local minimum
                if lower_neighbors == []:
                    parents[current_pt] = current_pt

                # some lower neighbors exist
                else:

                    # find parent pg of lower neighbors with lowest function value
                    neigh_parents = [find(n, parents) for n in lower_neighbors]
                    pg = neigh_parents[np.argmin([function[n] for n in neigh_parents])]

                    # set parent of current point to pg
                    parents[current_pt] = pg

                    # for each lower neighbor, we will create a persistence diagram point and corresponding set of nodes
                    for neighbor in lower_neighbors:

                        # get parent pn
                        pn = find(neighbor, parents)
                        val = function[pn]
                        persistence_set[pn] = []

                        # we will create persistence set only if parent pn is not local minimum pg
                        if pn != pg:
                            # go through all strictly lower nodes with parent pn
                            for v in sorted_idxs[:i]:
                                if find(v, parents) == pn:
                                    # if it is already part of another persistence set, continue
                                    try:
                                        visited[v]
                                    # else, mark visited and include it in current persistence set
                                    except KeyError:
                                        visited[v] = True
                                        persistence_set[pn].append(v)

                            # add current point to persistence set
                            persistence_set[pn].append(current_pt)

                            # do union and create persistence point corresponding to persistence set if persistence is sufficiently large
                            if np.abs(function[pn]-function[current_pt]) >= threshold:
                                persistence_diag[pn] = current_pt
                                union(pg, pn, parents, function)

            for key, val in iter(persistence_diag.items()):
                if topo_type == "downbranch":
                    dgm.append((0, (function[key],  function[val])))
                elif topo_type == "upbranch":
                    dgm.append((0, (-function[val], -function[key])))
                bnd.append(persistence_set[key])

        bnd = [list(b) for b in bnd]
        self.persistence_diagram, self.persistence_sets = dgm, bnd 
        return dgm, bnd

    def bootstrap_topological_features(self, N):
        """
        Use bootstrap to empirically assess stability of the features. This function computes a distribution of bottleneck distances, that can used afterwards to run tests on each topological feature.

        Parameters:
            N (int): number of bootstrap iterations.
        """
        if self.complex_type == "mapper":

            dgm = self.persistence_diagram
            num_pts, distribution = len(self.data), []
            for bootstrap_id in range(N):

                print(str(bootstrap_id) + "th iteration")

                # Randomly select points
                idxs = np.random.choice(num_pts, size=num_pts, replace=True)
                Xboot = self.data[idxs,:] if self.input_type == "point cloud" else self.data[idxs,:][:,idxs]
                f_boot, c_boot = self.filters[idxs,:], self.colors[idxs,:]
                Mboot = self.__class__(complex_type="mapper", filters=f_boot, filter_bnds=self.filter_bnds, colors=c_boot, resolutions=self.resolutions, gains=self.gains, 
                                      input_type=self.input_type, clustering=self.clustering).fit(Xboot)

                # Compute the corresponding persistence diagrams
                dgm_boot, _ = Mboot.compute_topological_features()

            # Compute the bottleneck distance
            npts, npts_boot = len(dgm), len(dgm_boot)
            D1 = np.array([[dgm[pt][1][0], dgm[pt][1][1]] for pt in range(npts)]) 
            D2 = np.array([[dgm_boot[pt][1][0], dgm_boot[pt][1][1]] for pt in range(npts_boot)])
            try:
                from gudhi.bottleneck import bottleneck_distance
                bottle = bottleneck_distance(D1, D2)
            except ImportError:
                print("Gudhi built without CGAL")
                raise
            distribution.append(bottle)
            self.distribution = np.sort(distribution)

        elif self.complex_type == "gic":

            self.complex.compute_PD()
            self.complex.compute_distribution(N)

    def get_distance_from_confidence_level(self, alpha=.95):
        """
        Compute the bottleneck distance threshold corresponding to a specific confidence level.

        Parameters:
            alpha (float): confidence level.

        Returns:
            distance value (float); each feature whose size is above this distance is sure at confidence level alpha.
        """
        if self.complex_type == "gic":
            return self.complex.compute_distance_from_confidence_level(alpha)
        elif self.complex_type == "mapper":
            return self.distribution[int(alpha*len(self.distribution))]

    def get_confidence_level_from_distance(self, distance):
        """
        Compute the confidence level of a specific bottleneck distance threshold.

        Parameters:
            distance (float): bottleneck distance threshold.

        Returns:
            confidence level (float); each feature whose size is above the distance threshold is sure at this confidence level.
        """
        if self.complex_type == "gic":
            return self.complex.compute_confidence_level_from_distance(distance)
        elif self.complex_type == "mapper":
            return len(np.argwhere(self.distribution <= distance))/len(self.distribution)

    def get_pvalue(self):
        """
        Compute the p-value, i.e. the opposite of the confidence level of the largest bottleneck distance preserving the topological features.

        Returns:
            p-value (float)
        """
        if self.complex_type == "gic":
            return self.complex.compute_p_value()
        elif self.complex_type == "mapper":
            distancemin = min([np.abs(pt[1][0]-pt[1][1]) for pt in self.persistence_diagram])
            return 1.-self.compute_confidence_from_distance(distancemin)

    def compute_differential_coordinates(self, nodes=None, features=None, sparse=False):
        """
        Compute the coordinates that best explain a set of nodes VS the rest of the nodes (in the 1-skeleton of the cover complex) with a Kolmogorov-Smirnov test. Only works if input_type is "point cloud".

        Parameters:
            nodes (list of integers): list of nodes to try. For instance, one can take the list of nodes obtained after calling "compute_topological_features"
            features (list of integers): the coordinates to try. All coordinates are tested if None.
            sparse (bool): set to True if your data is sparse and there will be speedup, otherwise use False.

        Returns:
            features (list of integers): the list of coordinates, ranked from smallest to largest p-values.
            p-values (list of float): the corresponding p-values. 
        """
        if self.input_type == "distance matrix":
            print("Need coordinates for running differential coordinates!")
            raise

        node_info = self.node_info
        X = self.data
        nodes = [s[0] for s,_ in self.simplex_tree.get_skeleton(0)] if nodes is None else nodes

        if features is None:
            features = np.arange(X.shape[1])

        list_idxs1 = list(np.unique(np.concatenate([node_info[node_name]["indices"] for node_name in nodes])))
        list_idxs2 = list(set(np.arange(X.shape[0]))-set(list_idxs1))
        pvals = []
        for f in features:
            if sparse:
                Xsp = csr_matrix(X)
                group1, group2 = np.squeeze(np.array(Xsp[list_idxs1,f].todense())), np.squeeze(np.array(Xsp[list_idxs2,f].todense()))
            else:
                group1, group2 = X[list_idxs1,f], X[list_idxs2,f]
            _,pval = ks_2samp(group1, group2)
            pvals.append(pval)
        pvals = np.array(pvals)
        F, P = features[np.argsort(pvals)], np.sort(pvals) 
        return F, P

    def print_to_dot(self, epsv=.2, epss=.4):
        """
        Write the cover complex in a DOT file, that can be processed with, e.g., neato.

        Parameters:
            epsv (float): scale the node colors between [epsv, 1-epsv]
            epss (float): scale the node sizes between [epss, 1-epss]
        """
        st = self.simplex_tree 
        node_info = self.node_info

        threshold = 0.
        maxv, minv = max([node_info[k]["colors"][0] for k in node_info.keys()]), min([node_info[k]["colors"][0] for k in node_info.keys()])
        maxs, mins = max([node_info[k]["size"]      for k in node_info.keys()]), min([node_info[k]["size"]      for k in node_info.keys()])  

        f = open(self.input_name + ".dot", "w")
        f.write("graph MAP{")
        cols = []
        for (simplex,_) in st.get_skeleton(0):
            cnode = (1.-2*epsv) * (node_info[simplex[0]]["colors"][0] - minv)/(maxv-minv) + epsv if maxv != minv else 0
            snode = (1.-2*epss) * (node_info[simplex[0]]["size"]-mins)/(maxs-mins) + epss if maxs != mins else 1
            f.write(  str(simplex[0]) + "[shape=circle width=" + str(snode) + " fontcolor=black color=black label=\""  + "\" style=filled fillcolor=\"" + str(cnode) + ", 1, 1\"]")
            cols.append(cnode)
        for (simplex,_) in st.get_filtration():
            if len(simplex) == 2:
                f.write("  " + str(simplex[0]) + " -- " + str(simplex[1]) + " [weight=15];")
        f.write("}")
        f.close()

        L = np.linspace(epsv, 1.-epsv, 100)
        colsrgb = []
        for c in L:	
            colsrgb.append(colorsys.hsv_to_rgb(c,1,1))
        fig, ax = plt.subplots(figsize=(6, 1))
        fig.subplots_adjust(bottom=0.5)
        my_cmap = matplotlib.colors.ListedColormap(colsrgb, name=self.color_name)
        cb = matplotlib.colorbar.ColorbarBase(ax, cmap=my_cmap, norm=matplotlib.colors.Normalize(vmin=minv, vmax=maxv), orientation="horizontal")
        cb.set_label(self.color_name)
        fig.savefig("colorbar_" + self.color_name + ".pdf", format="pdf")
        plt.close()

    def print_to_txt(self):
        """
        Write the cover complex to a TXT file, that can be processed with KeplerMapper.
        """
        st = self.simplex_tree
        if self.complex_type == "gic":
            self.complex.write_info(self.input_name.encode("utf-8"), self.cover_name.encode("utf-8"), self.color_name.encode("utf-8"))
        elif self.complex_type == "mapper":
            f = open(self.input_name + ".txt", "w")
            f.write(self.input_name + "\n")
            f.write(self.cover_name + "\n")
            f.write(self.color_name + "\n")
            f.write(str(self.complex.resolutions[0]) + " " + str(self.complex.gains[0]) + "\n")
            f.write(str(st.num_vertices()) + " " + str(len(list(st.get_skeleton(1)))-st.num_vertices()) + "\n")
            name2id = {}
            idv = 0
            for s,_ in st.get_skeleton(0):
                f.write(str(idv) + " " + str(self.node_info[s[0]]["colors"][0]) + " " + str(self.node_info[s[0]]["size"]) + "\n")
                name2id[s[0]] = idv
                idv += 1
            for s,_ in st.get_skeleton(1):
                if len(s) == 2:
                    f.write(str(name2id[s[0]]) + " " + str(name2id[s[1]]) + "\n")
            f.close()
