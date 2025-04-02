#!/usr/bin/env python

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2017 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__author__ = "Vincent Rouvreau"
__maintainer__ = ""
__copyright__ = "Copyright (C) 2017 Inria"
__license__ = "MIT"


import sys
import argparse
import gudhi


parser = argparse.ArgumentParser(
    description="RipsComplex creation from " "a correlation matrix read in a csv file.",
    epilog="Example: "
    "example/rips_complex_diagram_persistence_from_correlation_matrix_file_example.py "
    "-f ../data/correlation_matrix/lower_triangular_correlation_matrix.csv -e 12.0 -d 3"
    "- Constructs a Rips complex with the "
    "correlation matrix from the given csv file.",
)
parser.add_argument("-f", "--file", type=str, required=True)
parser.add_argument("-c", "--min_edge_correlation", type=float, default=0.5)
parser.add_argument("-d", "--max_dimension", type=int, default=1)
parser.add_argument("-b", "--band", type=float, default=0.0)
parser.add_argument(
    "--no-diagram",
    default=False,
    action="store_true",
    help="Flag for not to display the diagrams",
)

args = parser.parse_args()

if not (-1.0 < args.min_edge_correlation < 1.0):
    print("Wrong value of the threshold corelation (should be between -1 and 1).")
    sys.exit(1)

print("#####################################################################")
print("Caution: as persistence diagrams points will be under the diagonal,")
print("bottleneck distance and persistence graphical tool will not work")
print("properly, this is a known issue.")

print("#####################################################################")
print("RipsComplex creation from correlation matrix read in a csv file")

message = "RipsComplex with min_edge_correlation=" + repr(args.min_edge_correlation)
print(message)

correlation_matrix = gudhi.read_lower_triangular_matrix_from_csv_file(csv_file=args.file)
# Given a correlation matrix M, we compute component-wise M'[i,j] = 1-M[i,j] to get a distance matrix:
distance_matrix = [
    [1.0 - correlation_matrix[i][j] for j in range(len(correlation_matrix[i]))]
    for i in range(len(correlation_matrix))
]

rips_complex = gudhi.RipsComplex(
    distance_matrix=distance_matrix, max_edge_length=1.0 - args.min_edge_correlation
)
simplex_tree = rips_complex.create_simplex_tree(max_dimension=args.max_dimension)

message = "Number of simplices=" + repr(simplex_tree.num_simplices())
print(message)

diag = simplex_tree.persistence()

print("betti_numbers()=")
print(simplex_tree.betti_numbers())

# invert the persistence diagram
invert_diag = [
    (diag[pers][0], (1.0 - diag[pers][1][0], 1.0 - diag[pers][1][1]))
    for pers in range(len(diag))
]

if args.no_diagram == False:
    import matplotlib.pyplot as plot

    gudhi.plot_persistence_diagram(invert_diag, band=args.band)
    plot.show()
