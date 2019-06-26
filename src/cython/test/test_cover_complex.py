from gudhi import CoverComplex

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2018 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2018 Inria"
__license__ = "MIT"


def test_empty_constructor():
    # Try to create an empty CoverComplex
    cover = CoverComplex()
    assert cover.__is_defined() == True

def test_non_existing_file_read():
    # Try to open a non existing file
    cover = CoverComplex()
    assert (cover.read_point_cloud('pouetpouettralala.toubiloubabdou') == False)

def test_files_creation():
    # Create test file
    cloud_file = open('cloud', 'w')
    cloud_file.write('nOFF\n3\n3 0 0\n0 0 0\n2 1 0\n4 0 0')
    cloud_file.close()
    cover_file = open('cover', 'w')
    cover_file.write('1\n2\n3')
    cover_file.close()
    graph_file = open('graph', 'w')
    graph_file.write('0 1\n0 2\n1 2')
    graph_file.close()

def test_nerve():
    nerve = CoverComplex()
    nerve.set_type('Nerve')
    assert (nerve.read_point_cloud('cloud') == True)
    nerve.set_color_from_coordinate()
    nerve.set_graph_from_file('graph')
    nerve.set_cover_from_file('cover')
    nerve.find_simplices()
    stree = nerve.create_simplex_tree()

    assert (stree.num_vertices() == 3)
    assert ((stree.num_simplices() - stree.num_vertices()) == 0)
    assert (stree.dimension() == 0)

def test_graph_induced_complex():
    gic = CoverComplex()
    gic.set_type('GIC')
    assert (gic.read_point_cloud('cloud') == True)
    gic.set_color_from_coordinate()
    gic.set_graph_from_file('graph')
    gic.set_cover_from_file('cover')
    gic.find_simplices()
    stree = gic.create_simplex_tree()

    assert (stree.num_vertices() == 3)
    assert ((stree.num_simplices() - stree.num_vertices()) == 4)
    assert (stree.dimension() == 2)

def test_voronoi_graph_induced_complex():
    gic = CoverComplex()
    gic.set_type('GIC')
    assert (gic.read_point_cloud('cloud') == True)
    gic.set_color_from_coordinate()
    gic.set_graph_from_file('graph')
    gic.set_cover_from_Voronoi(2)
    gic.find_simplices()
    stree = gic.create_simplex_tree()

    assert (stree.num_vertices() == 2)
    assert ((stree.num_simplices() - stree.num_vertices()) == 1)
    assert (stree.dimension() == 1)
