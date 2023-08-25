""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2018 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi import CoverComplex
from gudhi.cover_complex import MapperComplex, GraphInducedComplex, NerveComplex
import pytest
import numpy as np
from sklearn.cluster import AgglomerativeClustering

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2018 Inria"
__license__ = "MIT"


def test_empty_constructor():
    # Try to create an empty CoverComplex
    cover = CoverComplex()
    assert cover._is_defined() == True

def test_non_existing_file_read():
    # Try to open a non existing file
    cover = CoverComplex()
    with pytest.raises(FileNotFoundError):
        cover.read_point_cloud("pouetpouettralala.toubiloubabdou")

def test_files_creation():
    # Create test file
    cloud_file = open("cloud", "w")
    cloud_file.write("nOFF\n3\n3 0 0\n0 0 0\n2 1 0\n4 0 0")
    cloud_file.close()
    cover_file = open("cover", "w")
    cover_file.write("1\n2\n3")
    cover_file.close()
    graph_file = open("graph", "w")
    graph_file.write("0 1\n0 2\n1 2")
    graph_file.close()

def test_nerve():
    nerve = CoverComplex()
    nerve.set_type("Nerve")
    assert nerve.read_point_cloud("cloud") == True
    nerve.set_color_from_coordinate()
    nerve.set_graph_from_file("graph")
    nerve.set_cover_from_file("cover")
    nerve.find_simplices()
    stree = nerve.create_simplex_tree()

    assert stree.num_vertices() == 3
    assert (stree.num_simplices() - stree.num_vertices()) == 0
    assert stree.dimension() == 0

def test_graph_induced_complex():
    gic = CoverComplex()
    gic.set_type("GIC")
    assert gic.read_point_cloud("cloud") == True
    gic.set_color_from_coordinate()
    gic.set_graph_from_file("graph")
    gic.set_cover_from_file("cover")
    gic.find_simplices()
    stree = gic.create_simplex_tree()

    assert stree.num_vertices() == 3
    assert (stree.num_simplices() - stree.num_vertices()) == 4
    assert stree.dimension() == 2

def test_voronoi_graph_induced_complex():
    gic = CoverComplex()
    gic.set_type("GIC")
    assert gic.read_point_cloud("cloud") == True
    gic.set_color_from_coordinate()
    gic.set_graph_from_file("graph")
    gic.set_cover_from_Voronoi(2)
    gic.find_simplices()
    stree = gic.create_simplex_tree()

    assert stree.num_vertices() == 2
    assert (stree.num_simplices() - stree.num_vertices()) == 1
    assert stree.dimension() == 1

def test_cover_complex():

    #              x
    #              |
    #              x
    #              |
    #      x---x---x
    #      |       |
    #      x       x
    #      |       |
    #      x---x---x
    #      |
    #      x
    #      |
    #      x

    X = np.array([[1,1],[1,1.5],[1,2],[1,2.5],[1,3],[1.5,2],[1.5,3],[2,2],[2,2.5],[2,3],[2,3.5],[2,4]])
    F = np.array([[1,1,1,1,1,1.5,1.5,2,2,2,2,2],[1,1.5,2,2.5,3,2,3,2,2.5,3,3.5,4]]).T

    M = GraphInducedComplex(input_type="point cloud", cover="functional", min_points_per_node=0, filter_bnds=np.array([.5,4.5]), 
         resolution=4, gain=.3, graph='rips', rips_threshold=.6, verbose=True).fit(X, filter=F[:,1], color=None)

    assert list(M.simplex_tree_.get_filtration()) == [([0], 0.0), ([1], 0.0), ([0, 1], 0.0), ([2], 0.0), ([1, 2], 0.0), ([3], 0.0), ([2, 3], 0.0)]

    M = GraphInducedComplex(input_type="point cloud", cover="voronoi", voronoi_samples=2, min_points_per_node=0, graph='rips', rips_threshold=.6, verbose=True).fit(X, color=None)

    assert list(M.simplex_tree_.get_filtration()) == [([0], 0.0), ([1], 0.0), ([0, 1], 0.0)]

    M = MapperComplex(input_type="point cloud", min_points_per_node=0, filter_bnds=np.array([[.5,2.5],[.5,4.5]]), 
         resolutions=np.array([2,4]), gains=np.array([.3,.3]), clustering=AgglomerativeClustering(n_clusters=None, linkage='single', distance_threshold=.6)).fit(X, filters=F, colors=None)

    assert list(M.simplex_tree_.get_filtration()) == [([0], 0.), ([1], 0.), ([0, 1], 0.), ([2], 0.), ([1, 2], 0.), ([3], 0.), ([1, 3], 0.), ([4], 0.), ([2, 4], 0.), ([3, 4], 0.), ([5], 0.), ([4, 5], 0.)]

    M = NerveComplex(input_type="point cloud", min_points_per_node=0).fit(X, assignments=[[0],[0,1],[1],[1,2],[2],[1,3],[2,4],[3],[3,4],[4,5],[5],[5]], color=None)

    assert list(M.simplex_tree_.get_filtration()) == [([0], 0.), ([1], 0.), ([0, 1], 0.), ([2], 0.), ([1, 2], 0.), ([3], 0.), ([1, 3], 0.), ([4], 0.), ([2, 4], 0.), ([3, 4], 0.), ([5], 0.), ([4, 5], 0.)]
