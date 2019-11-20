#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plot
import gudhi

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

parser = argparse.ArgumentParser(
    description="RipsComplex creation from " "a distance matrix read in a csv file.",
    epilog="Example: "
    "example/rips_complex_diagram_persistence_from_distance_matrix_file_example.py "
    "-f ../data/distance_matrix/lower_triangular_distance_matrix.csv -e 12.0 -d 3"
    "- Constructs a Rips complex with the "
    "distance matrix from the given csv file.",
)
parser.add_argument("-f", "--file", type=str, required=True)
parser.add_argument("-e", "--max_edge_length", type=float, default=0.5)
parser.add_argument("-d", "--max_dimension", type=int, default=1)
parser.add_argument("-b", "--band", type=float, default=0.0)
parser.add_argument(
    "--no-diagram",
    default=False,
    action="store_true",
    help="Flag for not to display the diagrams",
)

args = parser.parse_args()

print("#####################################################################")
print("RipsComplex creation from distance matrix read in a csv file")

message = "RipsComplex with max_edge_length=" + repr(args.max_edge_length)
print(message)

distance_matrix = gudhi.read_lower_triangular_matrix_from_csv_file(csv_file=args.file)
rips_complex = gudhi.RipsComplex(
    distance_matrix=distance_matrix, max_edge_length=args.max_edge_length
)
simplex_tree = rips_complex.create_simplex_tree(max_dimension=args.max_dimension)

message = "Number of simplices=" + repr(simplex_tree.num_simplices())
print(message)

diag = simplex_tree.persistence()

print("betti_numbers()=")
print(simplex_tree.betti_numbers())

if args.no_diagram == False:
    gudhi.plot_persistence_diagram(diag, band=args.band)
    plot.show()
