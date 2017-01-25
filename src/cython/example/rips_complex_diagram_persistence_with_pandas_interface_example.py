#!/usr/bin/env python

import gudhi
import pandas
import argparse

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Vincent Rouvreau

   Copyright (C) 2016 INRIA

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
__copyright__ = "Copyright (C) 2016 INRIA"
__license__ = "GPL v3"

print("#####################################################################")
print("RipsComplex creation from points read in a file")

parser = argparse.ArgumentParser(description='RipsComplex creation from '
                                 'points read in a file.',
                                 epilog='Example: '
                                 'example/rips_complex_diagram_persistence_with_pandas_interface_example.py '
                                 '../data/2000_random_points_on_3D_Torus.csv '
                                 '- Constructs a rips complex with the '
                                 'points from the given file. File format '
                                 'is X1, X2, ..., Xn')
parser.add_argument("-f", "--file", type=str, required=True)
parser.add_argument("-e", "--max-edge-length", type=float, default=0.5)
parser.add_argument('--no-diagram', default=False, action='store_true' , help='Flag for not to display the diagrams')

args = parser.parse_args()

points = pandas.read_csv(args.file, header=None)

message = "RipsComplex with max_edge_length=" + repr(args.max_edge_length)
print(message)

rips_complex = gudhi.RipsComplex(points=points.values,
                                 max_edge_length=args.max_edge_length)

simplex_tree = rips_complex.create_simplex_tree(max_dimension=len(points.values[0]))

message = "Number of simplices=" + repr(simplex_tree.num_simplices())
print(message)

diag = simplex_tree.persistence()

print("betti_numbers()=")
print(simplex_tree.betti_numbers())

if args.no_diagram == False:
    gudhi.diagram_persistence(diag)
