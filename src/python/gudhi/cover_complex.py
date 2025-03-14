# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Mathieu Carrière
#
# Copyright (C) 2021 Inria
#
# Modification(s):
#   - 2025/03 Ziyad Oulhaj: Added parallel processing for the MapperComplex fit method 

import numpy as np
import itertools
import matplotlib
import matplotlib.pyplot as plt
from joblib import Parallel, delayed

from sklearn.base            import BaseEstimator
from sklearn.cluster         import DBSCAN, AgglomerativeClustering
from sklearn.metrics         import pairwise_distances
from scipy.spatial.distance  import directed_hausdorff

from sklearn import __version__ as sklearn_version
from sklearn.utils.fixes import parse_version
agglomerative_clustering_metric = parse_version(sklearn_version) >= parse_version('1.2.0')

from . import SimplexTree, CoverComplex

# Takes the Binned data map and the index of a patch and computes the list of clusters in that patch. It is used for parallel processing in the MapperComplex fit method.
def cluster_patch(Binned_data, p_ind, clustering, X, input_type):
    data_bin=Binned_data[p_ind]
    if len(data_bin) > 1:
        clusters = clustering.fit_predict(X[data_bin,:]) if input_type == "point cloud" \
        else clustering.fit_predict(X[data_bin,:][:,data_bin])
    else:
        clusters = np.array([])

    return clusters

def _save_to_html(dat, lens, color, param, points, edges, html_output_filename):

    from ._kepler_mapper import KeplerMapper
    network = {}
    mapper = KeplerMapper(verbose=0)
    data = np.zeros((3,3))
    projected_data = mapper.fit_transform( data, projection="sum", scaler=None )

    from collections import defaultdict
    nodes = defaultdict(list)
    links = defaultdict(list)
    custom = defaultdict(list)

    for point in points:
        nodes[  str(int(point[0]))  ] = [  int(point[0]), point[1], int(point[2])  ]
        links[  str(int(point[0]))  ] = []
        custom[  int(point[0])  ] = point[1]

    for edge in edges:
        links[  str(edge[0])  ].append(  str(edge[1])  )
        links[  str(edge[1])  ].append(  str(edge[0])  )

    custom_val = custom.values()
    m = min(custom_val)
    M = max(custom_val)
    
    network["nodes"] = nodes
    network["links"] = links
    network["meta"] = lens


    mapper.visualize(network, color_function=color, path_html=html_output_filename, title=dat,
    graph_link_distance=30, graph_gravity=0.1, graph_charge=-120, custom_tooltips=custom, width_html=0,
    height_html=0, show_tooltips=True, show_title=True, show_meta=True, res=param[0], gain=param[1], minimum=m, maximum=M)
    message = repr(html_output_filename) + " is generated. You can now use your favorite web browser to visualize it."
    print(message)


class CoverComplexPy(BaseEstimator):
    """
    This is a mother class for MapperComplex, GraphInducedComplex and NerveComplex.

    Attributes:
        simplex_tree_ (gudhi SimplexTree): simplicial complex representing the cover complex computed after calling the fit() method.
        node_info_ (dictionary): various information associated to the nodes of the cover complex.
    """
    def __init__(self, verbose=False):
        """
        Constructor for the CoverComplexPy class.

        Parameters
        ----------
        """
        self.verbose = verbose

    def get_networkx(self, set_attributes_from_colors=False):
        """
        Turn the 1-skeleton of the cover complex computed after calling fit() method into a networkx graph.
        This function requires networkx (https://networkx.org/documentation/stable/install.html).

        Parameters
        ----------
        set_attributes_from_colors : bool
            if True, the color functions will be used as attributes for the networkx graph.

        Returns
        -------
        G : networkx graph
            graph representing the 1-skeleton of the cover complex.
        """
        import networkx as nx
        st = self.simplex_tree_
        G = nx.Graph()
        for (splx,_) in st.get_skeleton(1):	
            if len(splx) == 1:
                G.add_node(splx[0])
            if len(splx) == 2:
                G.add_edge(splx[0], splx[1])
        if set_attributes_from_colors:
            attrs = {k: {"attr_name": self.node_info_[k]["colors"]} for k in G.nodes()}
            nx.set_node_attributes(G, attrs)
        return G

    def save_to_dot(self, file_name="cover_complex", color_name="color", eps_color=.1, eps_size=.1):
        """
        Write the 0-skeleton of the cover complex in a DOT file called "{file_name}.dot", that can be processed with, e.g., neato. The vertices of the cover complex are colored with the first color function, ie, the first column of self.colors.  This function also produces an extra pdf file "colorbar_{color_name}.pdf" containing a colorbar corresponding to the node colors in the DOT file.

        Parameters
        ----------
        file_name : string
            name for the output .dot file, default "cover_complex" 
        color_name : string
            name for the output .pdf showing the colorbar of the color used for the Mapper nodes, default "color"
        eps_color : float
            scale the node colors between [eps_color, 1-eps_color]. Should be between 0 and 1/2. When close to 0., the color varies a lot across the nodes, if close to 1/2, the color tends to be more uniform.
        eps_size : float
            scale the node sizes between [eps_size, 1-eps_size]. Should be between 0 and 1/2. When close to 0., the size varies a lot across the nodes, if close to 1/2, the nodes tend to have the same size.
        """
        st = self.simplex_tree_
        node_info_ = self.node_info_

        maxv, minv = max(node_info_[k]["colors"][0] for k in node_info_.keys()), min(node_info_[k]["colors"][0] for k in node_info_.keys())
        maxs, mins = max(node_info_[k]["size"]      for k in node_info_.keys()), min(node_info_[k]["size"]      for k in node_info_.keys())

        if not file_name.lower().endswith(".dot"):
            file_name += ".dot"

        with open(file_name, "w") as f:
            f.write("graph MAP{")
            cols = []
            for (simplex,_) in st.get_skeleton(0):
                cnode = (1.-2*eps_color) * (node_info_[simplex[0]]["colors"][0] - minv)/(maxv-minv) + eps_color if maxv != minv else 0
                snode = (1.-2*eps_size) * (node_info_[simplex[0]]["size"]-mins)/(maxs-mins) + eps_size if maxs != mins else 1
                f.write(  str(simplex[0]) + "[shape=circle width=" + str(snode) + " fontcolor=black color=black label=\""  + "\" style=filled fillcolor=\"" + str(cnode) + ", 1, 1\"]")
                cols.append(cnode)
            for (simplex,_) in st.get_simplices():
                if len(simplex) == 2:
                    f.write("  " + str(simplex[0]) + " -- " + str(simplex[1]) + " [weight=15];")
            f.write("}")
        
        L = np.linspace(eps_color, 1.-eps_color, 100)
        colsrgb = []
        import colorsys
        for c in L:
            colsrgb.append(colorsys.hsv_to_rgb(c,1,1))
        fig, ax = plt.subplots(figsize=(6, 1))
        fig.subplots_adjust(bottom=0.5)
        my_cmap = matplotlib.colors.ListedColormap(colsrgb, name=color_name)
        cb = matplotlib.colorbar.ColorbarBase(ax, cmap=my_cmap, norm=matplotlib.colors.Normalize(vmin=minv, vmax=maxv), orientation="horizontal")
        cb.set_label(color_name)
        fig.savefig("colorbar_" + color_name + ".pdf", format="pdf")
        plt.close()
        
    def save_to_txt(self, file_name="cover_complex", data_name="data", cover_name="cover", color_name="color"):
        """
        Write the cover complex to a TXT file called "{file_name}.txt", that can be processed with the KeplerMapper Python script "KeplerMapperVisuFromTxtFile.py" available under "src/Nerve_GIC/utilities/".

        Parameters
        ----------
        file_name : string
            name for the output .txt file, default "cover_complex" 
        data_name : string
            name to use for the data on which the cover complex was computed, default "data". It will be used when generating an html visualization with KeplerMapperVisuFromTxtFile.py 
        cover_name : string
            name to use for the cover used to compute the cover complex, default "cover". It will be used when generating an html visualization with KeplerMapperVisuFromTxtFile.py
        color_name : string
            name to use for the color used to color the cover complex nodes, default "color". It will be used when generating an html visualization with KeplerMapperVisuFromTxtFile.py
        """
        st = self.simplex_tree_

        if not file_name.lower().endswith(".txt"):
            file_name += ".txt"

        with open(file_name, "w") as f:
            f.write(data_name + "\n")
            f.write(cover_name + "\n")
            f.write(color_name + "\n")
            f.write(str(self.resolutions[0]) + " " + str(self.gains[0]) + "\n")
            f.write(str(st.num_vertices()) + " " + str(len(list(st.get_skeleton(1)))-st.num_vertices()) + "\n")
            name2id = {}
            idv = 0
            for s,_ in st.get_skeleton(0):
                f.write(str(idv) + " " + str(self.node_info_[s[0]]["colors"][0]) + " " + str(self.node_info_[s[0]]["size"]) + "\n")
                name2id[s[0]] = idv
                idv += 1
            for s,_ in st.get_skeleton(1):
                if len(s) == 2:
                    f.write(str(name2id[s[0]]) + " " + str(name2id[s[1]]) + "\n")
    
    def save_to_html(self, file_name="cover_complex", data_name="data", cover_name="cover", color_name="color"):
        """
        Write the cover complex to an HTML file called "{file_name}.html", that can be visualized in a browser. This function is based on a fork of https://github.com/MLWave/kepler-mapper

        Parameters
        ----------
        file_name : string
            name for the output .html file, default "cover_complex" 
        data_name : string
            name to use for the data on which the cover complex was computed, default "data".
        cover_name : string
            name to use for the cover used to compute the cover complex, default "cover".
        color_name : string
            name to use for the color used to color the cover complex nodes, default "color".
        """

        st = self.simplex_tree_

        if not file_name.lower().endswith(".html"):
            file_name += ".html"

        points, edges, name2id, idv = [], [], {}, 0
        for s,_ in st.get_skeleton(0):
            points.append([ idv, self.node_info_[s[0]]["colors"][0], self.node_info_[s[0]]["size"] ])
            name2id[s[0]] = idv
            idv += 1
        for s,_ in st.get_skeleton(1):
            if len(s) == 2:
                edges.append([ name2id[s[0]] , name2id[s[1]] ])

        _save_to_html(data_name, cover_name, color_name, [self.resolutions[0], self.gains[0]], points, edges, file_name)


    class _constant_clustering():
        def fit_predict(X):
            return np.zeros([len(X)], dtype=np.int32)


class MapperComplex(CoverComplexPy):
    """
    This is a class for computing Mapper simplicial complexes on point clouds or distance matrices.
    """
    def __init__(self, *, input_type="point cloud", colors=None, min_points_per_node=0, filter_bnds=None, resolutions=None, gains=None, clustering=DBSCAN(), N=100, beta=0., C=10., verbose=False):
        """
        Constructor for the MapperComplex class.

        Parameters
        ----------
        input_type : string
            type of input data. Either "point cloud" or "distance matrix".
        min_points_per_node : int
            threshold on the size of the cover complex nodes (default 0). Any node associated to a subpopulation with less than **min_points_per_node** points will be removed.
        filter_bnds : list of lists or numpy array of shape (num_filters) x 2)
            limits of each filter, of the form [[f_1^min, f_1^max], ..., [f_n^min, f_n^max]]. If one of the values is numpy.nan, it can be computed from the dataset with the fit() method.
        resolutions : list or numpy array of shape num_filters containing integers
            resolution of each filter function, ie number of intervals required to cover each filter image. If None, it is estimated from data.
        gains : list or numpy array of shape num_filters containing doubles in [0,1]
            gain of each filter function, ie overlap percentage of the intervals covering each filter image. If None, it is set as 1/3 for all filters, since in the automatic parameter selection method in http://www.jmlr.org/papers/volume19/17-291/17-291.pdf, any arbitrary value between 1/3 and 1/2 works, so we go with the minimal one (ensuring that the complex is a graph if only given one filter).
        clustering : class
            clustering class (default sklearn.cluster.DBSCAN()). Common clustering classes can be found in the scikit-learn library (such as AgglomerativeClustering for instance). If None, it is set to hierarchical clustering, with scale estimated from data.
        N : int
            subsampling iterations (default 100) for estimating scale and resolutions. Used only if clustering or resolutions = None. See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
        beta : float
            exponent parameter (default 0.) for estimating scale and resolutions. Used only if clustering or resolutions = None. See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
        C : float 
            constant parameter (default 10.) for estimating scale and resolutions. Used only if clustering or resolutions = None. See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
        verbose : bool
            whether to display info while computing.
        """

        self.filter_bnds, self.resolutions, self.gains, self.clustering = filter_bnds, resolutions, gains, clustering
        self.input_type, self.min_points_per_node, self.N, self.beta, self.C = input_type, min_points_per_node, N, beta, C
        CoverComplexPy.__init__(self, verbose)

    def estimate_scale(self, X, N=100, beta=0., C=10.):
        """
        Compute estimated scale of a point cloud or a distance matrix.

        Parameters
        ----------
            X : numpy array of shape (num_points) x (num_coordinates) if point cloud and (num_points) x (num_points) if distance matrix
                input point cloud or distance matrix.
            N : int
                subsampling iterations (default 100). See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
            beta : float 
                exponent parameter (default 0.). See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
            C : float
                constant parameter (default 10.). See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.

        Returns
        -------
        delta : float
            estimated scale that can be used with, e.g., agglomerative clustering.
        """
        num_pts = X.shape[0]
        delta, m = 0., int(  num_pts / np.exp((1+beta) * np.log(np.log(num_pts)/np.log(C)))  )
        for _ in range(N):
            subpop = np.random.choice(num_pts, size=m, replace=False)
            if self.input_type == "point cloud":
                d, _, _ = directed_hausdorff(X, X[subpop,:])
            if self.input_type == "distance matrix":
                d = np.max(np.min(X[:,subpop], axis=1), axis=0)
            delta += d/N
        return delta

    def get_optimal_parameters_for_agglomerative_clustering(self, X, beta=0., C=10., N=100):
        """
        Compute optimal scale and resolutions for a point cloud or a distance matrix.

        Parameters
        ----------
            X : numpy array of shape (num_points) x (num_coordinates) if point cloud and (num_points) x (num_points) if distance matrix
                input point cloud or distance matrix.
            beta : float
                exponent parameter (default 0.). See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
            C : float 
                constant parameter (default 10.). See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
            N : int
                subsampling iterations (default 100). See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.

        Returns
        -------
        delta : float
            optimal scale that can be used with agglomerative clustering.
        resolutions : numpy array of shape (num_filters)
            optimal resolutions associated to each filter.
        """
        num_filt, delta = self.filters.shape[1], 0
        delta = self.estimate_scale(X=X, N=N, C=C, beta=beta)

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

    def fit(self, X, y=None, filters=None, colors=None, n_jobs=1, backend='loky'):
        """
        Fit the MapperComplex class on a point cloud or a distance matrix: compute the Mapper complex and store it in a simplex tree called `simplex_tree_`.

        Parameters
        ----------
            X : numpy array of shape (num_points) x (num_coordinates) if point cloud and (num_points) x (num_points) if distance matrix
                input point cloud or distance matrix.
            y : n x 1 array
                point labels (unused).
            filters : list of lists or numpy array of shape (num_points) x (num_filters) 
                filter functions (sometimes called lenses) used to compute the cover. Each column of the numpy array defines a scalar function defined on the input points.
            colors : list of lists or numpy array of shape (num_points) x (num_colors)
                functions used to color the nodes of the cover complex. More specifically, coloring is done by computing the means of these functions on the subpopulations corresponding to each node. If None, first coordinate is used if input is point cloud, and eccentricity is used if input is distance matrix.
            n_jobs :  The maximum number of concurrently running jobs (default 1). See https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html for details.
            backend : Specify the parallelization backend implementation (default 'loky'). See https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html for details.
        """

        if self.resolutions is not None:
            self.resolutions = np.array(self.resolutions)
            if len(self.resolutions.shape) == 0:
                self.resolutions = np.broadcast_to(self.resolutions, 1)
        if self.gains is not None:
            self.gains = np.array(self.gains)
            if len(self.gains.shape) == 0:
                self.gains = np.broadcast_to(self.gains, 1)
        if self.filter_bnds is not None:
            self.filter_bnds = np.array(self.filter_bnds)

        self.filters, self.colors = filters, colors

        if self.filters is None:
            if self.input_type == "point cloud":
                self.filters = X[:,0:1]
            elif self.input_type == "distance matrix":
                self.filters = X.max(axis=0)[:,None]
        else:
            if isinstance(self.filters, np.ndarray) == False:
                self.filters = np.array(self.filters).T

        if self.colors is None:
            if self.input_type == "point cloud":
                self.colors = X[:,0:1]
            elif self.input_type == "distance matrix":
                self.colors = X.max(axis=0)[:,None]
        else:
            if isinstance(self.colors, np.ndarray) == False:
                self.colors = np.array(self.colors).T

        if len(self.filters.shape) == 1: # if self.filters is a 1D filter, convert it to an array of shape [n,1]
            self.filters = np.reshape(self.filters, [len(X),1])
        if len(self.colors.shape) == 1: # if self.colors is a 1D filter, convert it to an array of shape [n,1]
            self.colors = np.reshape(self.colors, [len(X),1])
      
        num_pts, num_filters = self.filters.shape[0], self.filters.shape[1]

        # If some filter limits are unspecified, automatically compute them
        if self.filter_bnds is None:
            self.filter_bnds = np.hstack([np.min(self.filters, axis=0)[:,np.newaxis], np.max(self.filters, axis=0)[:,np.newaxis]])

        # If some resolutions are not specified, automatically compute them
        if self.gains is None:
            self.gains = np.broadcast_to(1/3, num_filters)
        if self.resolutions is None or self.clustering is None:
            delta, resolutions = self.get_optimal_parameters_for_agglomerative_clustering(X=X, beta=self.beta, C=self.C, N=self.N)
            if self.clustering is None:
                self.clustering = AgglomerativeClustering(n_clusters=None, linkage="single", distance_threshold=delta, **{
        "metric" if agglomerative_clustering_metric else "affinity":
        "euclidean" if self.input_type == "point cloud" else "precomputed"
                                                        })
            if self.resolutions is None:
                self.resolutions = np.multiply(resolutions, 1./self.gains)
                self.resolutions = np.array([int( (self.filter_bnds[ir,1]-self.filter_bnds[ir,0])/r) for ir, r in enumerate(self.resolutions)])

        # Initialize attributes
        self.simplex_tree_, self.node_info_ = SimplexTree(), {}

        # Compute the endpoints of the cover intervals for all filters 
        column_indices=np.indices((num_filters,self.resolutions.max()))[1]

        steps=np.repeat(((self.filter_bnds[:,1]-self.filter_bnds[:,0])/self.resolutions).reshape((num_filters,1)),self.resolutions.max(),axis=1)

        mins=np.repeat(self.filter_bnds[:,0].reshape((num_filters,1)),self.resolutions.max(),axis=1)

        gains=np.repeat(self.gains.reshape((num_filters,1)),self.resolutions.max(),axis=1)

        epsilons=gains/(2-2*gains)*steps

        left_endpoints=mins+steps*column_indices-epsilons

        right_endpoints=mins+steps*(column_indices+1)+epsilons
        
        # Find for each point the intervals it is associated to
        filter_values=np.repeat(self.filters.T.reshape(num_filters,1,num_pts),self.resolutions.max(),axis=1)

        comparison=(filter_values>=np.repeat(left_endpoints.reshape((num_filters,self.resolutions.max(),1)),num_pts,axis=2))*(filter_values<np.repeat(right_endpoints.reshape((num_filters,self.resolutions.max(),1)),num_pts,axis=2))

        # Compute the list of possible patches
        patch_list=[list(range(self.resolutions[f])) for f in range(num_filters)]
        patches=list(itertools.product(*patch_list))
        
        # Initialize the Binned data map that associates each patch to the points that belong to it
        # Initialize the cover map that associates each point to the clusters it belongs to
        Binned_data=[np.array([]) for p in patches]
        cover_map=[[] for pt in range(num_pts)]
        
        # Fill the Binned data map
        for p_ind in range(len(patches)):

            p=patches[p_ind]
            point_list=[set(np.where(comparison[f,p[f],:])[0]) for f in range(num_filters)]

            Binned_data_current=point_list[0]
            for f in range(1,num_filters):
                Binned_data_current=Binned_data_current.intersection(point_list[f])
            Binned_data[p_ind]=np.array(list(Binned_data_current))

        
        # Compute the clustering in each patch in parallel
        clusters_list = Parallel(n_jobs=n_jobs,backend=backend)(delayed(cluster_patch)(Binned_data, p_ind, self.clustering, X, self.input_type) for p_ind in range(len(patches)))
        
        # Go through the list of clusters
        current_max=0
        for clusters_ind in range(len(clusters_list)):
            # Change the name of the clusters to avoid confusion
            clusters_list[clusters_ind]=clusters_list[clusters_ind]+current_max
            clusters=clusters_list[clusters_ind]
            current_max+=len(np.unique(clusters))
            
            # Get information about each individual cluster
            for clus in np.unique(clusters):
                subpopulation = Binned_data[clusters_ind][clusters == clus]
                self.node_info_[clus] = {}
                self.node_info_[clus]["indices"] = subpopulation
                self.node_info_[clus]["size"] = len(subpopulation)
                self.node_info_[clus]["colors"] = np.mean(self.colors[subpopulation,:], axis=0)
                self.node_info_[clus]["patch"] = patches[clusters_ind]
            # Fill the cover map
            [cover_map[pt].append(clus) for pt,clus in zip(Binned_data[clusters_ind], clusters)]

        # Insert the simplices into the Mapper
        for splx in cover_map:
            self.simplex_tree_.insert(splx)

        return self


class GraphInducedComplex(CoverComplexPy):
    """
    This is a class for computing graph induced simplicial complexes on point clouds or distance matrices.
    """
    def __init__(self, *, input_type="point cloud", cover="functional", min_points_per_node=0,
                          voronoi_samples=100, assignments=None,  filter_bnds=None, resolution=None, gain=None, N=100, beta=0., C=10.,
                          graph="rips", rips_threshold=None, verbose=False):
        """
        Constructor for the GraphInducedComplex class.

        Parameters:
            input_type (string): type of input data. Either "point cloud" or "distance matrix".
            cover (string): specifies the cover. Either "functional" (preimages of filter function), "voronoi" or "precomputed".
            min_points_per_node (int): threshold on the size of the cover complex nodes (default 0). Any node associated to a subpopulation with less than **min_points_per_node** points will be removed.
            voronoi_samples (int): number of Voronoi germs used for partitioning the input dataset. Used only if cover = "voronoi".
            assignments (list of length (num_points) of lists of integers): cover assignment for each point. Used only if cover = "precomputed".
            filter_bnds (list or numpy array of shape 2): limits of the filter function, of the form [f^min, f^max]. If one of the values is numpy.nan, it can be computed from the dataset with the fit() method. Used only if cover = "functional".
            resolution (int): resolution of the filter function, ie number of intervals required to cover each filter image. Used only if cover = "functional". If None, it is estimated from data.
            gain (double in [0,1]): gain of the filter function, ie overlap percentage of the intervals covering each filter image. Used only if cover = "functional".
            N (int): subsampling iterations (default 100) for estimating scale and resolutions. Used only if cover = "functional". See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
            beta (double): exponent parameter (default 0.) for estimating scale and resolutions. Used only if cover = "functional". See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
            C (double): constant parameter (default 10.) for estimating scale and resolutions. Used only if cover = "functional". See http://www.jmlr.org/papers/volume19/17-291/17-291.pdf for details.
            graph (string): type of graph to use for GIC. Currently accepts "rips" only.
            rips_threshold (float): Rips parameter. Used only if graph = "rips".
            verbose (bool): whether to display info while computing.
        """

        self.input_type, self.cover, self.min_points_per_node = input_type, cover, min_points_per_node
        self.voronoi_samples, self.assignments, self.filter_bnds, self.resolution, self.gain = voronoi_samples, assignments, filter_bnds, resolution, gain
        self.graph, self.rips_threshold, self.N, self.beta, self.C = graph, rips_threshold, N, beta, C
        CoverComplexPy.__init__(self, verbose)

    def fit(self, X, y=None, filter=None, color=None):
        """
        Fit the GraphInducedComplex class on a point cloud or a distance matrix: compute the graph induced complex and store it in a simplex tree called `simplex_tree_`.

        Parameters:
            X (numpy array of shape (num_points) x (num_coordinates) if point cloud and (num_points) x (num_points) if distance matrix): input point cloud or distance matrix.
            y (n x 1 array): point labels (unused).
            filter (list or numpy array of shape (num_points)): filter function (sometimes called lens) used to compute the cover. Used only if cover = "functional".
            color (list or numpy array of shape (num_points)): function used to color the nodes of the cover complex. More specifically, coloring is done by computing the means of this function on the subpopulations corresponding to each node. If None, first coordinate is used if input is point cloud, and eccentricity is used if input is distance matrix.
        """
        self.data, self.filter, self.color = X, filter, color
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
        else:
            self.color = np.array(self.color)

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
            else:
                self.filter = np.array(self.filter)

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
        self.complex.set_mask(self.min_points_per_node)
        self.complex.find_simplices()
        simplex_tree_ = self.complex.create_simplex_tree()

        # Normalize vertex names of simplex tree
        self.simplex_tree_ = SimplexTree()
        idv, names = 0, {}
        for v,_ in simplex_tree_.get_skeleton(0):
            if len(self.complex.subpopulation(v[0])) > self.min_points_per_node:
                names[v[0]] = idv
                self.simplex_tree_.insert([idv])
                idv += 1
        for s,_ in simplex_tree_.get_simplices():
            if len(s) >= 2 and np.all([len(self.complex.subpopulation(v)) > self.min_points_per_node for v in s]):
                self.simplex_tree_.insert([names[v] for v in s])

        # Store vertex info
        self.node_info_ = {}
        for v,_ in simplex_tree_.get_skeleton(0):
            if len(self.complex.subpopulation(v[0])) > self.min_points_per_node:
                node = names[v[0]]
                pop = self.complex.subpopulation(v[0])
                self.node_info_[node] = {"indices": pop, "size": len(pop), "colors": [self.complex.subcolor(v[0])]}

        return self

class NerveComplex(CoverComplexPy):
    """
    This is a class for computing nerve simplicial complexes on point clouds or distance matrices.
    """
    def __init__(self, *, input_type="point cloud", min_points_per_node=0, verbose=False):
        """
        Constructor for the NerveComplex class.

        Parameters:
            input_type (string): type of input data. Either "point cloud" or "distance matrix".
            min_points_per_node (int): threshold on the size of the cover complex nodes (default 0). Any node associated to a subpopulation with less than **min_points_per_node** points will be removed.
            verbose (bool): whether to display info while computing.
        """

        self.input_type, self.min_points_per_node = input_type, min_points_per_node
        CoverComplexPy.__init__(self, verbose)

    def fit(self, X, y=None, assignments=None, color=None):
        """
        Fit the NerveComplex class on a point cloud or a distance matrix: compute the nerve complex and store it in a simplex tree called `simplex_tree_`.

        Parameters:
            X (numpy array of shape (num_points) x (num_coordinates) if point cloud and (num_points) x (num_points) if distance matrix): input point cloud or distance matrix.
            y (n x 1 array): point labels (unused).
            assignments (list of length (num_points) of lists of integers): cover assignment for each point.
            color (numpy array of shape (num_points) x (num_colors)): functions used to color the nodes of the cover complex. More specifically, coloring is done by computing the means of these functions on the subpopulations corresponding to each node. If None, first coordinate is used if input is point cloud, and eccentricity is used if input is distance matrix.
        """
        self.data, self.color = X, color
        self.assignments = assignments
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
        else:
            self.color = np.array(self.color)
        self.complex.set_color_from_range(self.color)

        # Set vertex cover
        self.complex.set_cover_from_range(self.assignments)

        # Compute simplex tree
        self.complex.set_mask(self.min_points_per_node)
        self.complex.find_simplices()
        simplex_tree_ = self.complex.create_simplex_tree()

        # Normalize vertex names of simplex tree
        self.simplex_tree_ = SimplexTree()
        idv, names = 0, {}
        for v,_ in simplex_tree_.get_skeleton(0):
            if len(self.complex.subpopulation(v[0])) > self.min_points_per_node:
                names[v[0]] = idv
                self.simplex_tree_.insert([idv])
                idv += 1
        for s,_ in simplex_tree_.get_simplices():
            if len(s) >= 2 and np.all([len(self.complex.subpopulation(v)) > self.min_points_per_node for v in s]):
                self.simplex_tree_.insert([names[v] for v in s])

        # Store vertex info
        self.node_info_ = {}
        for v,_ in simplex_tree_.get_skeleton(0):
            if len(self.complex.subpopulation(v[0])) > self.min_points_per_node:
                node = names[v[0]]
                pop = self.complex.subpopulation(v[0])
                self.node_info_[node] = {"indices": pop, "size": len(pop), "colors": [self.complex.subcolor(v[0])]}

        return self
