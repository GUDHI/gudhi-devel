#!/usr/bin/env python

import gudhi
import argparse

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

parser = argparse.ArgumentParser(description='Coordinate GIC '
                                 'from points read in a OFF file.',
                                 epilog='Example: '
                                 'example/coordinate_graph_induced_complex.py '
                                 '-f ../data/points/KleinBottle5D.off 0 -v'
                                 '- Constructs the coordinate GIC with the '
                                 'points from the given OFF file.')
parser.add_argument("-f", "--file", type=str, required=True)
parser.add_argument("-c", "--coordinate", type=int, default=0)
parser.add_argument("-v", "--verbose", default=False, action='store_true' , help='Flag for program verbosity')

args = parser.parse_args()

nerve_complex = gudhi.CoverComplex()
nerve_complex.set_verbose(args.verbose)

if (nerve_complex.read_point_cloud(args.file)):
    nerve_complex.set_type('GIC')
    nerve_complex.set_color_from_coordinate(args.coordinate)
    nerve_complex.set_function_from_coordinate(args.coordinate)
    nerve_complex.set_graph_from_automatic_rips()
    nerve_complex.set_automatic_resolution()
    nerve_complex.set_gain()
    nerve_complex.set_cover_from_function()
    nerve_complex.find_simplices()
    nerve_complex.plot_dot()
    simplex_tree = nerve_complex.create_simplex_tree()
    nerve_complex.compute_PD()
    if (args.verbose):
        print('Iterator on coordinate GIC simplices')
        result_str = 'Coordinate GIC is of dimension ' + \
            repr(simplex_tree.dimension()) + ' - ' + \
            repr(simplex_tree.num_simplices()) + ' simplices - ' + \
            repr(simplex_tree.num_vertices()) + ' vertices.'
        print(result_str)
        for filtered_value in simplex_tree.get_filtration():
            print(filtered_value[0])
