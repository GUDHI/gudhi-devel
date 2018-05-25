#!/usr/bin/env python

import gudhi
import argparse

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Vincent Rouvreau

   Copyright (C) 2016 Inria

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
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "GPL v3"

parser = argparse.ArgumentParser(description='RipsComplex creation from '
                                 'a distance matrix read in a csv file.',
                                 epilog='Example: '
                                 'example/rips_complex_diagram_persistence_from_distance_matrix_file_example.py '
                                 '-f ../data/distance_matrix/lower_triangular_distance_matrix.csv -e 12.0 -d 3'
                                 '- Constructs a Rips complex with the '
                                 'distance matrix from the given csv file.')
parser.add_argument("-f", "--file", type=str, required=True)
parser.add_argument("-e", "--max_edge_length", type=float, default=0.5)
parser.add_argument("-d", "--max_dimension", type=int, default=1)
parser.add_argument("-b", "--band_boot", type=float, default=0.)
parser.add_argument('--no-diagram', default=False, action='store_true' , help='Flag for not to display the diagrams')

args = parser.parse_args()

print("#####################################################################")
print("RipsComplex creation from distance matrix read in a csv file")

message = "RipsComplex with max_edge_length=" + repr(args.max_edge_length)
print(message)

distance_matrix = gudhi.read_lower_triangular_matrix_from_csv_file(csv_file=args.file)
rips_complex = gudhi.RipsComplex(distance_matrix=distance_matrix, max_edge_length=args.max_edge_length)
simplex_tree = rips_complex.create_simplex_tree(max_dimension=args.max_dimension)

message = "Number of simplices=" + repr(simplex_tree.num_simplices())
print(message)

diag = simplex_tree.persistence()

print("betti_numbers()=")
print(simplex_tree.betti_numbers())

if args.no_diagram == False:
    pplot = gudhi.plot_persistence_diagram(diag, band_boot=args.band_boot)
    pplot.show()
