# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Mathieu Carrière
#
# Copyright (C) 2021 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

import numpy as np
import itertools
import matplotlib
import matplotlib.pyplot as plt

from sklearn.base            import BaseEstimator, TransformerMixin
from sklearn.cluster         import DBSCAN, AgglomerativeClustering
from sklearn.metrics         import pairwise_distances
from scipy.spatial.distance  import directed_hausdorff

from .. import SimplexTree, CoverComplex

class CoverComplexPy(BaseEstimator, TransformerMixin):
    """
    This is a mother class for MapperComplex, GraphInducedComplex and NerveComplex. 

    Attributes:
        simplex_tree (gudhi SimplexTree): simplicial complex representing the cover complex computed after calling the fit() method.
        node_info (dictionary): various information associated to the nodes of the cover complex. 
    """
    def __init__(self, input_name="data", cover_name="cover", color_name="color", verbose=False):
        """
        Constructor for the CoverComplexPy class.

        Parameters:
        """
        self.input_name, self.cover_name, self.color_name, self.verbose = input_name, cover_name, color_name, verbose

    def get_networkx(self, get_attrs=False):
        """
        Turn the 1-skeleton of the cover complex computed after calling fit() method into a networkx graph.
        This function requires networkx (https://networkx.org/documentation/stable/install.html).

        Parameters:
            get_attrs (bool): if True, the color functions will be used as attributes for the networkx graph.

        Returns:
            G (networkx graph): graph representing the 1-skeleton of the cover complex.
        """
        try:
            import networkx as nx
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
        except ImportError:
            print("Networkx not found, nx graph not computed")

    def print_to_dot(self, epsv=.2, epss=.4):
        """
        Write the cover complex in a DOT file called "{self.input_name}.dot", that can be processed with, e.g., neato.

        Parameters:
            epsv (float): scale the node colors between [epsv, 1-epsv]
            epss (float): scale the node sizes between [epss, 1-epss]
        """
        st = self.simplex_tree 
        node_info = self.node_info

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
        try:
            import colorsys
            for c in L:
                colsrgb.append(colorsys.hsv_to_rgb(c,1,1))
            fig, ax = plt.subplots(figsize=(6, 1))
            fig.subplots_adjust(bottom=0.5)
            my_cmap = matplotlib.colors.ListedColormap(colsrgb, name=self.color_name)
            cb = matplotlib.colorbar.ColorbarBase(ax, cmap=my_cmap, norm=matplotlib.colors.Normalize(vmin=minv, vmax=maxv), orientation="horizontal")
            cb.set_label(self.color_name)
            fig.savefig("colorbar_" + self.color_name + ".pdf", format="pdf")
            plt.close()
        except ImportError:
            print("colorsys not found, colorbar not printed")

    def print_to_txt(self):
        """
        Write the cover complex to a TXT file called "{self.input_name}.txt", that can be processed with KeplerMapper.
        """
        st = self.simplex_tree
        f = open(self.input_name + ".txt", "w")
        f.write(self.input_name + "\n")
        f.write(self.cover_name + "\n")
        f.write(self.color_name + "\n")
        f.write("0 0\n")
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

    class _constant_clustering():
        def fit_predict(X):
            return np.zeros([len(X)], dtype=np.int32)


class MapperComplex(CoverComplexPy, BaseEstimator, TransformerMixin):
    """
    This is a class for computing Mapper simplicial complexes on point clouds or distance matrices.
    """
    def __init__(self, input_type="point cloud", colors=None, mask=0,
                       filters=None, filter_bnds=None, resolutions=None, gains=None, clustering=DBSCAN(), N=100, beta=0., C=10.,
                       input_name="data", cover_name="cover", color_name="color", verbose=False):
        """
        Constructor for the MapperComplex class.

        Parameters:
            input_type (string): type of input data. Either "point cloud" or "distance matrix".
            colors (numpy array of shape (num_points) x (num_colors)): functions used to color the nodes of the cover complex. More specifically, coloring is done by computing the means of these functions on the subpopulations corresponding to each node. If None, first coordinate is used if input is point cloud, and eccentricity is used if input is distance matrix.
            mask (int): threshold on the size of the cover complex nodes (default 0). Any node associated to a subpopulation with less than **mask** points will be removed.
            filters (numpy array of shape (num_points) x (num_filters)): filter functions (sometimes called lenses) used to compute the cover. Each column of the numpy array defines a scalar function defined on the input points.
            filter_bnds (numpy array of shape (num_filters) x 2): limits of each filter, of the form [[f_1^min, f_1^max], ..., [f_n^min, f_n^max]]. If one of the values is numpy.nan, it can be computed from the dataset with the fit() method.
            resolutions (numpy array of shape num_filters containing integers): resolution of each filter function, ie number of intervals required to cover each filter image. If None, it is estimated from data.
            gains (numpy array of shape num_filters containing doubles in [0,1]): gain of each filter function, ie overlap percentage of the intervals covering each filter image.
            clustering (class): clustering class (default sklearn.cluster.DBSCAN()). Common clustering classes can be found in the scikit-learn library (such as AgglomerativeClustering for instance). If None, it is set to hierarchical clustering, with scale estimated from data.
            N (int): subsampling iterations (default 100) for estimating scale and resolutions. Used only if clustering or resolutions = None. See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
            beta (double): exponent parameter (default 0.) for estimating scale and resolutions. Used only if clustering or resolutions = None. See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
            C (double): constant parameter (default 10.) for estimating scale and resolutions. Used only if clustering or resolutions = None. See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
            input_name (string): name of dataset. Used when generating plots.
            cover_name (string): name of cover. Used when generating plots.
            color_name (string): name of color function. Used when generating plots.
            verbose (bool): whether to display info while computing.
        """
        self.filters, self.filter_bnds, self.resolutions, self.gains, self.colors, self.clustering = filters, filter_bnds, resolutions, gains, colors, clustering
        self.input_type, self.mask, self.N, self.beta, self.C = input_type, mask, N, beta, C
        CoverComplexPy.__init__(self, input_name, cover_name, color_name, verbose)

    def estimate_scale(self, X, N=100, inp="point cloud", beta=0., C=10.):
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
        num_filt, delta = self.filters.shape[1], 0
        delta = self.estimate_scale(X=X, N=N, inp=self.input_type, C=C, beta=beta)

        pairwise = pairwise_distances(X, metric="euclidean") if self.input_type == "point cloud" else X
        pairs = np.argwhere(pairwise <= delta)
        num_pairs = pairs.shape[0]
        res = []
        for f in range(num_filt):
            F = self.filters[:,f]
            resf = 0
            for p in range(num_pairs):
                resf = max(resf, abs(F[pairs[p,0]] - F[pairs[p,1]]))
            res.append(resf)

        return delta, np.array(res)

    def fit(self, X, y=None):
        """
        Fit the MapperComplex class on a point cloud or a distance matrix: compute the Mapper complex and store it in a simplex tree called simplex_tree.

        Parameters:
            X (numpy array of shape (num_points) x (num_coordinates) if point cloud and (num_points) x (num_points) if distance matrix): input point cloud or distance matrix.
            y (n x 1 array): point labels (unused).
        """
        if self.filters is None:
            if self.input_type == "point cloud":
                self.filters = X[:,0:1]
            elif self.input_type == "distance matrix":
                self.filters = X.max(axis=0)[:,None]
        if self.colors is None:
            if self.input_type == "point cloud":
                self.colors = X[:,0:1]
            elif self.input_type == "distance matrix":
                self.colors = X.max(axis=0)[:,None]

        num_pts, num_filters = self.filters.shape[0], self.filters.shape[1]

        # If some filter limits are unspecified, automatically compute them
        if self.filter_bnds is None: 
            self.filter_bnds = np.hstack([np.min(self.filters, axis=0)[:,np.newaxis], np.max(self.filters, axis=0)[:,np.newaxis]])

        # If some resolutions are not specified, automatically compute them
        if self.gains is None:
            self.gains = .33 * np.ones(num_filters)
        if self.resolutions is None or self.clustering is None:
            delta, resolutions = self.get_optimal_parameters_for_agglomerative_clustering(X=X, beta=self.beta, C=self.C, N=self.N)
            if self.clustering is None:
                if self.input_type == "point cloud":
                    self.clustering = AgglomerativeClustering(n_clusters=None, linkage="single", distance_threshold=delta, affinity="euclidean")  
                else:
                    self.clustering = AgglomerativeClustering(n_clusters=None, linkage="single", distance_threshold=delta, affinity="precomputed")
            if self.resolutions is None:
                self.resolutions = np.multiply(resolutions, 1./self.gains)
                self.resolutions = np.array([int( (self.filter_bnds[ir,1]-self.filter_bnds[ir,0])/r) for ir, r in enumerate(self.resolutions)])

        # Initialize attributes
        self.simplex_tree, self.node_info = SimplexTree(), {}

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
                clusters = self.clustering.fit_predict(X[idxs,:]) if self.input_type == "point cloud" else self.clustering.fit_predict(X[idxs,:][:,idxs])
            elif len(idxs) == 1:
                clusters = np.array([0])
            else:
                continue

            # Collect various information on each cluster
            num_clus_pre = np.max(clusters) + 1
            for clus_i in range(num_clus_pre):
                node_name = clus_base + clus_i
                subpopulation = idxs[clusters == clus_i]
                self.node_info[node_name] = {}
                self.node_info[node_name]["indices"] = subpopulation
                self.node_info[node_name]["size"] = len(subpopulation)
                self.node_info[node_name]["colors"] = np.mean(self.colors[subpopulation,:], axis=0)
                self.node_info[node_name]["patch"] = preimage

            # Update the cover map
            for pt in range(clusters.shape[0]):
                node_name = clus_base + clusters[pt]
                if clusters[pt] != -1 and self.node_info[node_name]["size"] >= self.mask:
                    cover[idxs[pt]].append(node_name)

            clus_base += np.max(clusters) + 1

        # Insert the simplices of the Mapper complex 
        for i in range(num_pts):
            self.simplex_tree.insert(cover[i])

        return self

class GraphInducedComplex(CoverComplexPy, BaseEstimator, TransformerMixin):
    """
    This is a class for computing graph induced simplicial complexes on point clouds or distance matrices.
    """
    def __init__(self, input_type="point cloud", cover="functional", color=None, mask=0,
                       voronoi_samples=100, assignments=None, filter=None, filter_bnds=None, resolution=None, gain=None, N=100, beta=0., C=10.,
                       graph="rips", rips_threshold=None,
                       input_name="data", cover_name="cover", color_name="color", verbose=False):
        """
        Constructor for the GraphInducedComplex class.

        Parameters:
            input_type (string): type of input data. Either "point cloud" or "distance matrix".
            cover (string): specifies the cover. Either "functional" (preimages of filter function), "voronoi" or "precomputed".
            color (numpy array of shape (num_points)): function used to color the nodes of the cover complex. More specifically, coloring is done by computing the means of this function on the subpopulations corresponding to each node. If None, first coordinate is used if input is point cloud, and eccentricity is used if input is distance matrix.
            mask (int): threshold on the size of the cover complex nodes (default 0). Any node associated to a subpopulation with less than **mask** points will be removed.
            voronoi_samples (int): number of Voronoi germs used for partitioning the input dataset. Used only if cover = "voronoi".
            assignments (list of length (num_points) of lists of integers): cover assignment for each point. Used only if cover = "precomputed".
            filter (numpy array of shape (num_points)): filter function (sometimes called lens) used to compute the cover. Used only if cover = "functional".
            filter_bnds (numpy array of shape 2): limits of the filter function, of the form [f^min, f^max]. If one of the values is numpy.nan, it can be computed from the dataset with the fit() method. Used only if cover = "functional".
            resolution (int): resolution of the filter function, ie number of intervals required to cover each filter image. Used only if cover = "functional". If None, it is estimated from data.
            gain (double in [0,1]): gain of the filter function, ie overlap percentage of the intervals covering each filter image. Used only if cover = "functional".
            N (int): subsampling iterations (default 100) for estimating scale and resolutions. Used only if cover = "functional". See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
            beta (double): exponent parameter (default 0.) for estimating scale and resolutions. Used only if cover = "functional". See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
            C (double): constant parameter (default 10.) for estimating scale and resolutions. Used only if cover = "functional". See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
            graph (string): type of graph to use for GIC. Currently accepts "rips" only.
            rips_threshold (float): Rips parameter. Used only if graph = "rips".
            input_name (string): name of dataset. Used when generating plots.
            cover_name (string): name of cover. Used when generating plots.
            color_name (string): name of color function. Used when generating plots.
            verbose (bool): whether to display info while computing.
        """

        self.input_type, self.cover, self.color, self.mask = input_type, cover, color, mask
        self.voronoi_samples, self.assignments, self.filter, self.filter_bnds, self.resolution, self.gain = voronoi_samples, assignments, filter, filter_bnds, resolution, gain
        self.graph, self.rips_threshold, self.N, self.beta, self.C = graph, rips_threshold, N, beta, C
        CoverComplexPy.__init__(self, input_name, cover_name, color_name, verbose)

    def fit(self, X, y=None):
        """
        Fit the GraphInducedComplex class on a point cloud or a distance matrix: compute the graph induced complex and store it in a simplex tree called simplex_tree.

        Parameters:
            X (numpy array of shape (num_points) x (num_coordinates) if point cloud and (num_points) x (num_points) if distance matrix): input point cloud or distance matrix.
            y (n x 1 array): point labels (unused).
        """
        self.data = X
        self.complex = CoverComplex()
        self.complex.set_type("GIC")
        self.complex.set_verbose(self.verbose)

        if self.input_type == "point cloud":
            self.complex.set_point_cloud_from_range(X)
        elif self.input_type == "distance matrix":
            self.complex.set_distances_from_range(X)
            
        # Set vertex color
        if self.color is None:
            if self.input_type == "point cloud":
                self.color = X[:,0]
            elif self.input_type == "distance matrix":
                self.color = X.max(axis=0)
        self.complex.set_color_from_range(self.color)
            
        # Set underlying graph
        if self.graph == "rips":
            if self.rips_threshold is not None:
                self.complex.set_graph_from_rips(self.rips_threshold)
            else:
                self.complex.set_subsampling(self.C, self.beta)
                self.complex.set_graph_from_automatic_rips(self.N)

        # Set vertex cover
        if self.cover == "voronoi":
            self.complex.set_cover_from_Voronoi(self.voronoi_samples)

        elif self.cover == "functional":

            if self.filter is None:
                if self.input_type == "point cloud":
                    self.filter = X[:,0]
                elif self.input_type == "distance matrix":
                    self.filter = X.max(axis=0)
            self.complex.set_function_from_range(self.filter)

            if self.resolution is None:
                self.complex.set_automatic_resolution()
            else:
                self.complex.set_resolution_with_interval_number(self.resolution)

            if self.gain is None:
                self.complex.set_gain(.33)
            else:
                self.complex.set_gain(self.gain)

            self.complex.set_cover_from_function()

        elif self.cover == "precomputed":
            self.complex.set_cover_from_range(self.assignments)

        # Compute simplex tree
        self.complex.set_mask(self.mask)
        self.complex.find_simplices()
        simplex_tree = self.complex.create_simplex_tree()
            
        # Normalize vertex names of simplex tree
        self.simplex_tree = SimplexTree()
        idv, names = 0, {}
        for v,_ in simplex_tree.get_skeleton(0):
            if len(self.complex.subpopulation(v[0])) > self.mask:
                names[v[0]] = idv
                self.simplex_tree.insert([idv])
                idv += 1
        for s,_ in simplex_tree.get_simplices():
            if len(s) >= 2 and np.all([len(self.complex.subpopulation(v)) > self.mask for v in s]):
                self.simplex_tree.insert([names[v] for v in s])

        # Store vertex info
        self.node_info = {}
        for v,_ in simplex_tree.get_skeleton(0):
            if len(self.complex.subpopulation(v[0])) > self.mask:
                node = names[v[0]]
                pop = self.complex.subpopulation(v[0])
                self.node_info[node] = {"indices": pop, "size": len(pop), "colors": [self.complex.subcolor(v[0])]}
    
        return self

class NerveComplex(CoverComplexPy, BaseEstimator, TransformerMixin):
    """
    This is a class for computing nerve simplicial complexes on point clouds or distance matrices.
    """
    def __init__(self, input_type="point cloud", color=None, mask=0,
                       assignments=None,
                       input_name="data", cover_name="cover", color_name="color", verbose=False):
        """
        Constructor for the NerveComplex class.

        Parameters:
            input_type (string): type of input data. Either "point cloud" or "distance matrix".
            color (numpy array of shape (num_points) x (num_colors)): functions used to color the nodes of the cover complex. More specifically, coloring is done by computing the means of these functions on the subpopulations corresponding to each node. If None, first coordinate is used if input is point cloud, and eccentricity is used if input is distance matrix.
            mask (int): threshold on the size of the cover complex nodes (default 0). Any node associated to a subpopulation with less than **mask** points will be removed.
            assignments (list of length (num_points) of lists of integers): cover assignment for each point.
            input_name (string): name of dataset. Used when generating plots.
            cover_name (string): name of cover. Used when generating plots.
            color_name (string): name of color function. Used when generating plots.
            verbose (bool): whether to display info while computing.
        """

        self.input_type, self.color, self.mask = input_type, color, mask
        self.assignments = assignments
        CoverComplexPy.__init__(self, input_name, cover_name, color_name, verbose)

    def fit(self, X, y=None):
        """
        Fit the NerveComplex class on a point cloud or a distance matrix: compute the nerve complex and store it in a simplex tree called simplex_tree.

        Parameters:
            X (numpy array of shape (num_points) x (num_coordinates) if point cloud and (num_points) x (num_points) if distance matrix): input point cloud or distance matrix.
            y (n x 1 array): point labels (unused).
        """
        self.data = X
        self.complex = CoverComplex()
        self.complex.set_type("Nerve")
        self.complex.set_verbose(self.verbose)

        if self.input_type == "point cloud":
            self.complex.set_point_cloud_from_range(X)
        elif self.input_type == "distance matrix":
            self.complex.set_distances_from_range(X)
            
        # Set vertex color
        if self.color is None:
            if self.input_type == "point cloud":
                self.color = X[:,0]
            elif self.input_type == "distance matrix":
                self.color = X.max(axis=0)
        self.complex.set_color_from_range(self.color)
            
        # Set vertex cover
        self.complex.set_cover_from_range(self.assignments)

        # Compute simplex tree
        self.complex.set_mask(self.mask)
        self.complex.find_simplices()
        simplex_tree = self.complex.create_simplex_tree()
            
        # Normalize vertex names of simplex tree
        self.simplex_tree = SimplexTree()
        idv, names = 0, {}
        for v,_ in simplex_tree.get_skeleton(0):
            if len(self.complex.subpopulation(v[0])) > self.mask:
                names[v[0]] = idv
                self.simplex_tree.insert([idv])
                idv += 1
        for s,_ in simplex_tree.get_simplices():
            if len(s) >= 2 and np.all([len(self.complex.subpopulation(v)) > self.mask for v in s]):
                self.simplex_tree.insert([names[v] for v in s])

        # Store vertex info
        self.node_info = {}
        for v,_ in simplex_tree.get_skeleton(0):
            if len(self.complex.subpopulation(v[0])) > self.mask:
                node = names[v[0]]
                pop = self.complex.subpopulation(v[0])
                self.node_info[node] = {"indices": pop, "size": len(pop), "colors": [self.complex.subcolor(v[0])]}
    
        return self

