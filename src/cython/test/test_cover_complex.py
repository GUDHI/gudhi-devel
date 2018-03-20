from gudhi import CoverComplex

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Vincent Rouvreau

   Copyright (C) 2018 Inria

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2018 Inria"
__license__ = "GPL v3"


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
    assert (nerve.read_point_cloud('cloud') == True)
    nerve.set_type('Nerve')
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
    gic.set_graph_from_file('graph')
    gic.set_cover_from_Voronoi(2)
    gic.find_simplices()
    stree = gic.create_simplex_tree()

    assert (stree.num_vertices() == 2)
    assert ((stree.num_simplices() - stree.num_vertices()) == 1)
    assert (stree.dimension() == 1)
